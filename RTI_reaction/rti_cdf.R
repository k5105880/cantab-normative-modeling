# Deterministic new-subject marginal prediction for Student/log1p RTI models.
# For every posterior draw, b_i ~ N(0, z_i' Sigma_b z_i) is integrated with
# Gaussian quadrature. Student residuals are integrated analytically with pt().
# Draws are averaged on the probability/density scale using log-sum-exp.
# No finite predictive sample, empirical CDF, or probability clipping is used.

.rti_scalar_integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != as.integer(x) || x < lower)
    stop(name, " must be an integer >= ", lower, call. = FALSE)
  as.integer(x)
}

# Golub-Welsch nodes/weights for E[f(Z)], Z ~ N(0,1). The probabilists'
# Hermite Jacobi matrix has off-diagonal sqrt(1), sqrt(2), ... .
# Caching avoids repeated eigendecomposition while scoring data chunks.
rti_gauss_normal <- local({
  cache <- new.env(parent = emptyenv())
  function(nodes = 40L) {
    nodes <- .rti_scalar_integer(nodes, "nodes", 2L)
    key <- as.character(nodes)
    if (exists(key, envir = cache, inherits = FALSE)) return(cache[[key]])
    J <- matrix(0, nodes, nodes)
    J[cbind(seq_len(nodes - 1L), 2:nodes)] <- sqrt(seq_len(nodes - 1L))
    J <- J + t(J)
    e <- eigen(J, symmetric = TRUE)
    ord <- order(e$values)
    weights <- e$vectors[1L, ord]^2
    # Very remote nodes can have numerically zero weights; they contribute zero.
    keep <- weights > 0
    out <- list(nodes = e$values[ord][keep],
                weights = weights[keep] / sum(weights[keep]))
    cache[[key]] <- out
    out
  }
})

.rti_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

.rti_logmean_cols <- function(x) {
  m <- apply(x, 2L, max)
  shifted <- sweep(x, 2L, m, "-")
  out <- m + log(colMeans(exp(shifted)))
  out[is.infinite(m) & m < 0] <- -Inf
  # Sums of normalized weights may round a log probability above zero by ulps;
  # probability outputs are bounded below, after integration, only at exactly 1.
  out
}

.rti_z_from_logs <- function(log_lower, log_upper) {
  ifelse(log_lower <= log_upper,
         qnorm(log_lower, lower.tail = TRUE, log.p = TRUE),
         qnorm(log_upper, lower.tail = FALSE, log.p = TRUE))
}

# Pure numerical component. mu, sigma, and re_variance are draws x rows;
# response is one log1p observation per row, and nu one finite positive df/draw.
# Returns the posterior-averaged marginal density on the log1p response scale.
rti_cdf_components <- function(response, mu, sigma, nu, re_variance,
                               quad_nodes = 40L, keep_draws = FALSE) {
  if (!is.matrix(mu) || !is.numeric(mu) || any(dim(mu) < 1L))
    stop("mu must be a nonempty numeric draws-by-rows matrix", call. = FALSE)
  if (!is.numeric(response) || length(response) != ncol(mu) ||
      !is.numeric(nu) || length(nu) != nrow(mu) ||
      !is.matrix(sigma) || !is.numeric(sigma) ||
      !is.matrix(re_variance) || !is.numeric(re_variance) ||
      !identical(dim(sigma), dim(mu)) || !identical(dim(re_variance), dim(mu)))
    stop("response, mu, sigma, nu, and RE variance have incompatible dimensions", call. = FALSE)
  if (any(!is.finite(response)) || any(!is.finite(mu)) ||
      any(!is.finite(sigma)) || any(sigma <= 0) ||
      any(!is.finite(nu)) || any(nu <= 0) ||
      any(!is.finite(re_variance)) || any(re_variance < 0))
    stop("All parameters must be finite; sigma/nu positive and RE variance nonnegative", call. = FALSE)
  g <- rti_gauss_normal(quad_nodes)
  centered <- matrix(response, nrow(mu), ncol(mu), byrow = TRUE) - mu
  re_sd <- sqrt(re_variance)
  ll <- lu <- ld <- matrix(-Inf, nrow(mu), ncol(mu))
  for (k in seq_along(g$nodes)) {
    q <- (centered - re_sd * g$nodes[k]) / sigma
    lw <- log(g$weights[k])
    ll <- .rti_logadd(ll, matrix(pt(q, df = nu, lower.tail = TRUE, log.p = TRUE),
                                nrow(mu), ncol(mu)) + lw)
    lu <- .rti_logadd(lu, matrix(pt(q, df = nu, lower.tail = FALSE, log.p = TRUE),
                                nrow(mu), ncol(mu)) + lw)
    ld <- .rti_logadd(ld, matrix(dt(q, df = nu, log = TRUE), nrow(mu), ncol(mu)) -
                       log(sigma) + lw)
  }
  log_lower <- pmin(0, .rti_logmean_cols(ll))
  log_upper <- pmin(0, .rti_logmean_cols(lu))
  out <- data.frame(log_p_lower = log_lower, log_p_upper = log_upper,
             z = .rti_z_from_logs(log_lower, log_upper),
             log_density = .rti_logmean_cols(ld))
  if(keep_draws) {
    attr(out,'log_lower_draws')<-pmin(ll,0)
    attr(out,'log_upper_draws')<-pmin(lu,0)
  }
  out
}

# Rows must contain all iterations of chain 1, then chain 2, etc. This is the
# ordering of posterior::as_draws_matrix and brms predictions with all draw IDs.
# Scaling the small tail avoids underflow; MCSE is transformed to z only after
# averaging probabilities. This is simulation error, not a subject-level CI.
rti_probability_precision <- function(log_lower,log_upper,chains,mcse_limit=.03) {
  chains<-.rti_scalar_integer(chains,'chains')
  if(!is.matrix(log_lower)||!identical(dim(log_lower),dim(log_upper))||
     any(dim(log_lower)<1L)||any(!is.finite(log_lower))||any(!is.finite(log_upper))||
     any(log_lower>1e-10)||any(log_upper>1e-10)||nrow(log_lower)%%chains!=0L||
     nrow(log_lower)/chains<4L)stop('Complete ordered finite probability chains are required')
  if(length(mcse_limit)!=1L||!is.finite(mcse_limit)||mcse_limit<=0)
    stop('Positive finite posterior MCSE limit required')
  lower<-.rti_logmean_cols(log_lower);upper<-.rti_logmean_cols(log_upper)
  if(any(abs(exp(lower)+exp(upper)-1)>1e-8))stop('Posterior tail probabilities do not sum to one')
  z<-.rti_z_from_logs(pmin(lower,0),pmin(upper,0))
  result<-lapply(seq_len(ncol(log_lower)),function(j) {
    tail<-if(lower[j]<=upper[j])log_lower[,j]else log_upper[,j]
    shift<-max(tail)
    p<-matrix(exp(tail-shift),ncol=chains)
    mc<-as.numeric(posterior::mcse_mean(p))
    ess<-as.numeric(posterior::ess_mean(p))
    se<-if(is.finite(mc)&&mc>=0&&is.finite(z[j]))
      exp(log(mc)+shift-dnorm(z[j],log=TRUE))else NA_real_
    valid<-is.finite(se)&&is.finite(ess)&&ess>0
    data.frame(posterior_mcse_z=se,posterior_tail_ess=ess,
      posterior_mc_ok=valid&&se<=mcse_limit,
      posterior_mc_status=if(!valid)'posterior_MC_undefined'else
        if(se>mcse_limit)'posterior_MC_unresolved'else'ok',
      posterior_n_chains=chains,posterior_iterations_per_chain=nrow(p))
  })
  do.call(rbind,result)
}

rti_apply_marginal_precision <- function(p) {
  if(!all(c('z','numerical_ok','posterior_mc_ok','posterior_mcse_z','posterior_mc_status')%in%names(p))||
     anyNA(p$numerical_ok)||anyNA(p$posterior_mc_ok))stop('Invalid marginal precision result')
  mc_ok<-p$posterior_mc_ok & is.finite(p$posterior_mcse_z)
  p$z_estimate<-p$z
  p$score_ok<-p$numerical_ok & is.finite(p$z) & mc_ok
  p$score_status<-ifelse(!p$numerical_ok,'numerical_unresolved',
    ifelse(!mc_ok,p$posterior_mc_status,'ok'))
  p$score_status[p$score_status=='ok' & !p$score_ok]<-'posterior_MC_undefined'
  p$z<-ifelse(p$score_ok,p$z,NA_real_)
  p
}

# Check CDF and density convergence independently: symmetry may give exactly
# z=0 at every quadrature order despite a badly underestimated density.
# Only failing rows are retried. Differences always compare the final estimate
# with its immediately preceding quadrature order, and the used order is exposed.
rti_cdf_checked <- function(response, mu, sigma, nu, re_variance,
                            quad_nodes = 40L, check_nodes = 80L,
                            z_tolerance = 0.02, log_score_tolerance = 0.01,
                            max_nodes = max(320L, check_nodes), keep_draws = FALSE) {
  quad_nodes <- .rti_scalar_integer(quad_nodes, "quad_nodes", 2L)
  check_nodes <- .rti_scalar_integer(check_nodes, "check_nodes", 2L)
  max_nodes <- .rti_scalar_integer(max_nodes, "max_nodes", 2L)
  if (check_nodes <= quad_nodes || max_nodes < check_nodes)
    stop("Require quad_nodes < check_nodes <= max_nodes", call. = FALSE)
  if (length(z_tolerance) != 1L || !is.finite(z_tolerance) || z_tolerance <= 0 ||
      length(log_score_tolerance) != 1L || !is.finite(log_score_tolerance) ||
      log_score_tolerance <= 0)
    stop("z_tolerance and log_score_tolerance must be finite and positive", call. = FALSE)
  compare <- function(coarse, fine) {
    zd <- abs(fine$z - coarse$z)
    ld <- abs(fine$log_density - coarse$log_density)
    ok <- is.finite(fine$z) & is.finite(fine$log_density) &
      is.finite(zd) & zd <= z_tolerance & is.finite(ld) & ld <= log_score_tolerance &
      is.finite(fine$log_p_lower) & is.finite(fine$log_p_upper) &
      abs(exp(fine$log_p_lower) + exp(fine$log_p_upper) - 1) < 1e-10
    data.frame(quad_z_diff = zd, quad_log_score_diff = ld, numerical_ok = ok)
  }
  coarse <- rti_cdf_components(response, mu, sigma, nu, re_variance, quad_nodes)
  fine <- rti_cdf_components(response, mu, sigma, nu, re_variance, check_nodes,keep_draws)
  lower_draws<-attr(fine,'log_lower_draws');upper_draws<-attr(fine,'log_upper_draws')
  diagnostics <- compare(coarse, fine)
  used <- rep.int(check_nodes, nrow(fine))
  current_nodes <- check_nodes
  while (any(!diagnostics$numerical_ok) && current_nodes < max_nodes) {
    failed <- which(!diagnostics$numerical_ok)
    next_nodes <- as.integer(min(current_nodes * 2, max_nodes))
    next_fine <- rti_cdf_components(response[failed], mu[, failed, drop = FALSE],
                                    sigma[, failed, drop = FALSE], nu,
                                    re_variance[, failed, drop = FALSE], next_nodes,keep_draws)
    if(keep_draws) {
      lower_draws[,failed]<-attr(next_fine,'log_lower_draws')
      upper_draws[,failed]<-attr(next_fine,'log_upper_draws')
    }
    diagnostics[failed, ] <- compare(fine[failed, , drop = FALSE], next_fine)
    fine[failed, ] <- next_fine
    used[failed] <- next_nodes
    current_nodes <- next_nodes
  }
  out<-cbind(fine, diagnostics, quad_nodes_used = used)
  if(keep_draws) {
    attr(out,'log_lower_draws')<-lower_draws
    attr(out,'log_upper_draws')<-upper_draws
  }
  out
}

.rti_validate_covariance <- function(cov_draws) {
  d <- dim(cov_draws)
  if (!is.numeric(cov_draws) || length(d) != 3L || d[1L] < 1L ||
      d[2L] != d[3L] || any(!is.finite(cov_draws)))
    stop("Subject_ID covariance must be a finite draws-by-coefficient-by-coefficient array", call. = FALSE)
  cn <- dimnames(cov_draws)[[2L]]
  if (is.null(cn) || anyDuplicated(cn) ||
      !identical(cn, dimnames(cov_draws)[[3L]]) ||
      !(setequal(cn, c("Intercept", "Time_since_baseline_1y")) ||
        setequal(cn, c("PlatformRS6", "PlatformConnect", "Time_since_baseline_1y"))))
    stop("Unsupported Subject_ID random effects: expected intercept/time or RS6/Connect/time", call. = FALSE)
  tolerance <- 1e-10 * max(1, max(abs(cov_draws)))
  if (max(abs(cov_draws - aperm(cov_draws, c(1, 3, 2)))) > tolerance)
    stop("Subject_ID covariance draws are not symmetric", call. = FALSE)
  for (s in seq_len(d[1L])) {
    if (min(eigen(cov_draws[s, , ], symmetric = TRUE, only.values = TRUE)$values) < -tolerance)
      stop("Subject_ID covariance draw is not positive semidefinite", call. = FALSE)
  }
  invisible(cn)
}

.rti_re_design <- function(coef_names, newdata) {
  tt <- newdata$Time_since_baseline_1y
  if (!is.numeric(tt) || length(tt) != nrow(newdata) || any(!is.finite(tt)))
    stop("Time_since_baseline_1y must be present and finite", call. = FALSE)
  z <- matrix(0, nrow(newdata), length(coef_names),
              dimnames = list(NULL, coef_names))
  z[, "Time_since_baseline_1y"] <- tt
  if ("Intercept" %in% coef_names) {
    z[, "Intercept"] <- 1
  } else {
    platform <- as.character(newdata$Platform)
    if (length(platform) != nrow(newdata) || anyNA(platform) ||
        any(!platform %in% c("RS6", "Connect")))
      stop("Platform must contain only RS6 or Connect for platform random effects", call. = FALSE)
    z[, "PlatformRS6"] <- as.numeric(platform == "RS6")
    z[, "PlatformConnect"] <- as.numeric(platform == "Connect")
  }
  z
}

.rti_design_variance <- function(cov_draws, z) {
  out <- matrix(0, dim(cov_draws)[1L], nrow(z))
  for (j in seq_len(ncol(z))) {
    for (k in seq_len(ncol(z))) {
      out <- out + outer(cov_draws[, j, k], z[, j] * z[, k])
    }
  }
  if (any(!is.finite(out)) || any(out < -1e-8))
    stop("Computed new-subject RE variances are invalid", call. = FALSE)
  pmax(out, 0) # preserve matrix shape; remove only roundoff in valid PSD matrices
}

# Public pure helper, also useful when independently auditing the RE covariance.
rti_re_variance <- function(cov_draws, newdata) {
  coef_names <- .rti_validate_covariance(cov_draws)
  .rti_design_variance(cov_draws, .rti_re_design(coef_names, newdata))
}

.rti_validate_fit_structure <- function(fit) {
  if (!inherits(fit, "brmsfit")) stop("fit must be a brmsfit", call. = FALSE)
  if (!identical(fit$family$family, "student") ||
      !identical(fit$family$link, "identity"))
    stop("Only Student models with identity mean link are supported", call. = FALSE)
  tt <- brms::brmsterms(fit$formula)
  if (!identical(fit$formula$formula[[2L]], as.name("NM_response")) ||
      isTRUE(attr(fit$formula$formula, "nl")) || length(tt$nlpars))
    stop("Expected an ordinary univariate NM_response formula", call. = FALSE)
  if (length(fit$formula$pfix) || !all(names(tt$dpars) %in% c("mu", "sigma")))
    stop("Only scalar posterior nu and an ordinary fixed-predictor sigma formula are supported", call. = FALSE)
  mu <- tt$dpars$mu
  sigma <- tt$dpars$sigma
  allowed_mu <- c("formula", "fe", "re", "sm", "allvars", "family", "dpar",
                  "resp", "respform", "transform")
  allowed_sigma <- c("formula", "fe", "allvars", "family", "dpar", "resp",
                     "respform", "transform")
  if (!all(names(mu) %in% allowed_mu) ||
      (!is.null(sigma) && !all(names(sigma) %in% allowed_sigma)))
    stop("Unsupported mean/sigma structure (e.g. random sigma, autocorrelation, or special effects)", call. = FALSE)
  re <- mu$re
  if (is.null(re) || nrow(re) != 1L || !identical(as.character(re$group), "Subject_ID"))
    stop("Exactly one Subject_ID grouping term is required", call. = FALSE)
  group_call <- re$gcall[[1L]]
  if (!identical(group_call$dist, "gaussian") ||
      any(vapply(c("by", "pw", "cov", "type"),
                 function(x) !identical(group_call[[x]], ""), logical(1))))
    stop("Only ordinary independent Gaussian subject effects are supported", call. = FALSE)
  invisible(TRUE)
}

# Scores a NEW subject, even if Subject_ID matches a training subject. Spline
# terms remain in posterior_linpred(re_formula=NA); only subject effects are
# excluded there and integrated separately through their covariance draws.
# quad_z_diff and quad_log_score_diff compare the final two quadrature orders.
# Failed rows retry at doubled orders up to max_nodes; numerical_ok requires
# both z and log-density convergence. Failed rows must be withheld from z
# analysis and likelihood comparisons. Probability columns may round to 0 or 1,
# while z uses independently evaluated log tails and is never clipped.
# The raw log score is log p(log1p(raw)) - log1p(raw), i.e. correct Jacobian.
rti_predict_marginal <- function(fit, newdata, draw_ids = NULL,
                                 quad_nodes = 40L, check_nodes = 80L,
                                 chunk_size = 40L, z_tolerance = 0.02,
                                 log_score_tolerance = 0.01,
                                 max_nodes = max(320L, check_nodes),
                                 posterior_mcse_limit = .03) {
  if (!requireNamespace("brms", quietly = TRUE) ||
      !requireNamespace("posterior", quietly = TRUE))
    stop("brms and posterior are required", call. = FALSE)
  .rti_validate_fit_structure(fit)
  quad_nodes <- .rti_scalar_integer(quad_nodes, "quad_nodes", 2L)
  check_nodes <- .rti_scalar_integer(check_nodes, "check_nodes", 2L)
  max_nodes <- .rti_scalar_integer(max_nodes, "max_nodes", 2L)
  chunk_size <- .rti_scalar_integer(chunk_size, "chunk_size", 1L)
  if (check_nodes <= quad_nodes || max_nodes < check_nodes)
    stop("Require quad_nodes < check_nodes <= max_nodes", call. = FALSE)
  if (length(z_tolerance) != 1L || !is.finite(z_tolerance) || z_tolerance <= 0 ||
      length(log_score_tolerance) != 1L || !is.finite(log_score_tolerance) ||
      log_score_tolerance <= 0)
    stop("z_tolerance and log_score_tolerance must be finite and positive", call. = FALSE)
  if (!is.data.frame(newdata) || !nrow(newdata) ||
      !all(c("NM_raw", "row_id_global") %in% names(newdata)))
    stop("newdata requires rows, NM_raw, and row_id_global", call. = FALSE)
  if (!is.numeric(newdata$NM_raw) || any(!is.finite(newdata$NM_raw)) ||
      any(newdata$NM_raw < 0) || anyNA(newdata$row_id_global))
    stop("RTI raw responses must be finite/nonnegative and row IDs nonmissing", call. = FALSE)
  response <- log1p(newdata$NM_raw)
  if ("NM_response" %in% names(newdata) &&
      (!is.numeric(newdata$NM_response) || any(!is.finite(newdata$NM_response)) ||
       max(abs(newdata$NM_response - response)) > 1e-8))
    stop("NM_response must equal log1p(NM_raw)", call. = FALSE)
  total_draws <- posterior::ndraws(fit)
  if (is.null(draw_ids)) draw_ids <- seq_len(total_draws)
  if (!is.numeric(draw_ids) || !length(draw_ids) || any(!is.finite(draw_ids)) ||
      any(draw_ids != as.integer(draw_ids)) || any(draw_ids < 1L | draw_ids > total_draws) ||
      anyDuplicated(draw_ids))
    stop("draw_ids must be distinct valid posterior draw indices", call. = FALSE)
  # brms::posterior_linpred selects draws in posterior order, even when IDs
  # arrive unsorted. Canonicalize before also subsetting covariance and nu,
  # otherwise a random draw subset silently breaks within-draw dependence.
  draw_ids <- sort(as.integer(draw_ids))
  complete_chains<-identical(draw_ids,seq_len(total_draws))
  chains<-posterior::nchains(fit)
  vc <- brms::VarCorr(fit, summary = FALSE)
  if (!identical(names(vc), "Subject_ID") || is.null(vc$Subject_ID$cov))
    stop("Exactly one Subject_ID covariance block is required", call. = FALSE)
  cov_draws <- vc$Subject_ID$cov[draw_ids, , , drop = FALSE]
  coef_names <- .rti_validate_covariance(cov_draws)
  design <- .rti_re_design(coef_names, newdata)
  nu_draws <- posterior::as_draws_matrix(fit, variable = "nu")
  if (!identical(colnames(nu_draws), "nu"))
    stop("Expected one scalar posterior nu parameter per draw", call. = FALSE)
  nu <- as.numeric(nu_draws[draw_ids, "nu"])
  if (any(!is.finite(nu)) || any(nu <= 0)) stop("Invalid posterior nu", call. = FALSE)
  result <- vector("list", ceiling(nrow(newdata) / chunk_size))
  for (part in seq_along(result)) {
    idx <- seq.int((part - 1L) * chunk_size + 1L, min(part * chunk_size, nrow(newdata)))
    nd <- newdata[idx, , drop = FALSE]
    # Avoid mgcv's singleton smooth-matrix dimension drop, then remove padding.
    model_nd<-if(nrow(nd)==1L)nd[c(1L,1L),,drop=FALSE]else nd
    keep<-seq_len(nrow(nd))
    mu <- brms::posterior_linpred(fit, newdata = model_nd, re_formula = NA,
                                 draw_ids = draw_ids, allow_new_levels = TRUE)
    sigma <- brms::posterior_linpred(fit, newdata = model_nd, re_formula = NA,
                                    dpar = "sigma", transform = TRUE,
                                    draw_ids = draw_ids, allow_new_levels = TRUE)
    mu<-mu[,keep,drop=FALSE];sigma<-sigma[,keep,drop=FALSE]
    v <- .rti_design_variance(cov_draws, design[idx, , drop = FALSE])
    fine <- rti_cdf_checked(response[idx], mu, sigma, nu, v, quad_nodes, check_nodes,
                            z_tolerance, log_score_tolerance, max_nodes,keep_draws=complete_chains)
    precision<-if(complete_chains)
      rti_probability_precision(attr(fine,'log_lower_draws'),attr(fine,'log_upper_draws'),
                                chains,posterior_mcse_limit)else
      data.frame(posterior_mcse_z=rep(NA_real_,nrow(nd)),posterior_tail_ess=NA_real_,
        posterior_mc_ok=FALSE,posterior_mc_status='not_assessed_draw_subset',
        posterior_n_chains=NA_integer_,posterior_iterations_per_chain=NA_integer_)
    raw_score <- fine$log_density - response[idx]
    if (any(!is.finite(raw_score)))
      stop("Nonfinite raw marginal log score; do not silently drop scoring rows", call. = FALSE)
    lower <- exp(fine$log_p_lower)
    upper <- exp(fine$log_p_upper)
    result[[part]] <- cbind(data.frame(row_id_global = nd$row_id_global,
                                p_lower = lower, p_upper = upper, z = fine$z,
                                log_score_raw = raw_score,
                                quad_z_diff = fine$quad_z_diff,
                                numerical_ok = fine$numerical_ok,
                                quad_log_score_diff = fine$quad_log_score_diff,
                                quad_nodes_used = fine$quad_nodes_used),precision)
  }
  ans <- do.call(rbind, result)
  rownames(ans) <- NULL
  ans
}
