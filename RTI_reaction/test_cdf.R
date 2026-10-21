# Run: Rscript --vanilla test_cdf.R
# Optional RTI_TEST_FIT points to a saved brms fit; no model is fitted.
test_file <- grep("^--file=", gsub("~+~", " ", commandArgs(), fixed=TRUE), value = TRUE)
test_dir <- if (length(test_file)) dirname(normalizePath(sub("^--file=", "", test_file[1]))) else getwd()
module_file <- file.path(test_dir, "rti_cdf.R")
if (!file.exists(module_file)) stop("FAIL: deterministic marginal CDF module is not implemented")
source(module_file)

passed <- 0L
check <- function(label, expr) {
  force(expr)
  passed <<- passed + 1L
  cat(sprintf("PASS %02d: %s\n", passed, label))
}
near <- function(actual, expected, tolerance = 1e-9) {
  stopifnot(length(actual) == length(expected), all(is.finite(actual)),
            max(abs(actual - expected)) < tolerance)
}
must_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

# Catches wrong normal-Hermite scaling or quadrature weight normalization.
check("normal quadrature integrates Gaussian moments", {
  g <- rti_gauss_normal(40L)
  near(sum(g$weights), 1)
  near(sum(g$weights * g$nodes), 0)
  near(sum(g$weights * g$nodes^2), 1)
  near(sum(g$weights * g$nodes^4), 3)
})

# Catches integrating only a mean draw or replacing the Student residual by Gaussian.
check("posterior mixture equals exact Student mixture with zero RE variance", {
  y <- c(-1.2, 0, 2.4)
  mu <- matrix(c(-0.3, 0.6), 2, 3)
  sigma <- matrix(c(0.7, 1.2), 2, 3)
  nu <- c(3, 11)
  got <- rti_cdf_components(y, mu, sigma, nu, matrix(0, 2, 3), 40L)
  lower <- upper <- density <- numeric(3)
  for (i in 1:3) {
    q <- (y[i] - mu[, i]) / sigma[, i]
    lower[i] <- mean(pt(q, df = nu))
    upper[i] <- mean(pt(q, df = nu, lower.tail = FALSE))
    density[i] <- mean(dt(q, df = nu) / sigma[, i])
  }
  near(exp(got$log_p_lower), lower)
  near(exp(got$log_p_upper), upper)
  near(exp(got$log_density), density)
})

# Adaptive integration is independent of Gaussian quadrature implementation.
check("Gaussian RE plus Student residual agrees with adaptive integration", {
  y <- c(-2.1, 0.25, 3.6)
  mus <- c(-0.3, 0.6)
  sigmas <- c(0.7, 1.2)
  nus <- c(4, 9)
  variances <- c(0.8, 0.35)
  got <- rti_cdf_components(y, matrix(mus, 2, 3), matrix(sigmas, 2, 3),
                            nus, matrix(variances, 2, 3), 120L)
  for (i in seq_along(y)) {
    exact <- sapply(c("lower", "upper", "density"), function(kind) {
      mean(vapply(1:2, function(j) integrate(function(b) {
        q <- (y[i] - mus[j] - b) / sigmas[j]
        v <- switch(kind, lower = pt(q, nus[j]),
                    upper = pt(q, nus[j], lower.tail = FALSE),
                    density = dt(q, nus[j]) / sigmas[j])
        v * dnorm(b, sd = sqrt(variances[j]))
      }, -Inf, Inf, rel.tol = 1e-11)$value, numeric(1)))
    })
    near(exp(unlist(got[i, c("log_p_lower", "log_p_upper", "log_density")])),
         exact, 2e-8)
  }
})

# Catches accidentally discarding new-subject RE variance.
check("large-df limit includes Gaussian RE variance", {
  y <- c(-3.2, -0.5, 1.3, 5.6)
  got <- rti_cdf_components(y, matrix(0.4, 1, 4), matrix(0.8, 1, 4),
                            1e9, matrix(1.7, 1, 4), 100L)
  sd_total <- sqrt(0.8^2 + 1.7)
  near(exp(got$log_p_lower), pnorm(y, 0.4, sd_total), 1e-8)
  near(exp(got$log_density), dnorm(y, 0.4, sd_total), 1e-8)
})

# Catches upper-tail cancellation and empirical-CDF clipping at about 4.75.
check("extreme tails remain finite and exceed absolute z=4.75", {
  got <- rti_cdf_components(c(-1e7, 1e7), matrix(0, 1, 2), matrix(1, 1, 2),
                            7, matrix(0.3, 1, 2), 80L)
  stopifnot(all(is.finite(got$z)), got$z[1] < -10, got$z[2] > 10,
            got$log_p_upper[2] < -50, got$log_p_lower[1] < -50)
  near(got$z[1], -got$z[2], 1e-10)
})

# A symmetric CDF can be exactly 0.5 while density quadrature is badly wrong.
# This catches treating z convergence as sufficient for likelihood comparison.
check("symmetric z convergence cannot conceal an unconverged log density", {
  stopifnot(exists("rti_cdf_checked", mode = "function"))
  for (s in c(0.1, 0.03)) {
    got <- rti_cdf_checked(0, matrix(0, 1, 1), matrix(s, 1, 1), 10,
                           matrix(1, 1, 1), max_nodes = 80L)
    exact <- log(integrate(function(b) dt(-b / s, 10) / s * dnorm(b),
                           -Inf, Inf, rel.tol = 1e-11)$value)
    stopifnot(got$quad_z_diff < 1e-10, got$quad_log_score_diff > 0.01,
              abs(got$log_density - exact) > 0.4, !got$numerical_ok,
              got$quad_nodes_used == 80L)
  }
})

# Independently integrate the successful retry, and keep an extreme case flagged.
check("adaptive quadrature retries failed rows and respects the node limit", {
  got <- rti_cdf_checked(c(0, 0, 0), matrix(0, 1, 3),
                         matrix(c(0.2, 0.03, 0.8), 1, 3), 10,
                         matrix(1, 1, 3), max_nodes = 320L)
  exact <- log(integrate(function(b) dt(-b / 0.2, 10) / 0.2 * dnorm(b),
                         -Inf, Inf, rel.tol = 1e-11)$value)
  stopifnot(identical(got$numerical_ok, c(TRUE, FALSE, TRUE)),
            identical(got$quad_nodes_used, c(320L, 320L, 80L)),
            got$quad_log_score_diff[1] < 0.01,
            got$quad_log_score_diff[2] > 0.01)
  near(got$log_density[1], exact, 1e-5)
  must_error(rti_cdf_checked(0, matrix(0, 1, 1), matrix(1, 1, 1),
                             10, matrix(1, 1, 1), log_score_tolerance = 0))
})

# Literal expectations use distinct platform variances and opposite slope covariances.
check("platform RE design preserves unequal variances and correlations", {
  cn <- c("PlatformRS6", "PlatformConnect", "Time_since_baseline_1y")
  cv <- matrix(c(.16, .02, .03, .02, .64, -.06, .03, -.06, .04), 3, 3)
  cov_draws <- array(0, c(2, 3, 3), dimnames = list(NULL, cn, cn))
  cov_draws[1, , ] <- cv
  cov_draws[2, , ] <- 2 * cv
  nd <- data.frame(Platform = c("RS6", "Connect", "RS6"),
                   Time_since_baseline_1y = c(2, 2, 0))
  got <- rti_re_variance(cov_draws, nd)
  stopifnot(identical(dim(got), c(2L, 3L)))
  near(as.vector(got), c(.44, .88, .56, 1.12, .16, .32))
})

check("intercept/slope RE design uses the covariance cross term", {
  cn <- c("Time_since_baseline_1y", "Intercept")
  cv <- array(c(.04, .03, .03, .16), c(1, 2, 2), dimnames = list(NULL, cn, cn))
  nd <- data.frame(Time_since_baseline_1y = c(0, 2, -2))
  near(as.vector(rti_re_variance(cv, nd)), c(.16, .44, .20))
})

check("invalid parameters and unsupported RE structures are rejected", {
  m <- matrix(1, 1, 1)
  must_error(rti_cdf_components(1, m, -m, 3, m, 40L))
  must_error(rti_cdf_components(1, m, m, NA_real_, m, 40L))
  must_error(rti_cdf_components(1, m, m, 3, -m, 40L))
  must_error(rti_cdf_components(Inf, m, m, 3, m, 40L))
  must_error(rti_cdf_components(1, m, m, 3, m, 1L))
  cv <- array(1, c(1, 1, 1), dimnames = list(NULL, "Age", "Age"))
  must_error(rti_re_variance(cv, data.frame(Time_since_baseline_1y = 0)))
  cv <- array(c(1, 2, 2, 1), c(1, 2, 2),
              dimnames = list(NULL, c("Intercept", "Time_since_baseline_1y"),
                              c("Intercept", "Time_since_baseline_1y")))
  must_error(rti_re_variance(cv, data.frame(Time_since_baseline_1y = 0)))
})

check("exact quadrature cannot conceal noisy posterior probability chains", {
  set.seed(651)
  mu<-matrix(rep(c(-3,-1,1,3),each=100)+rnorm(400,sd=.3),400,1)
  got<-rti_cdf_checked(0,mu,matrix(1,400,1),rep(7,400),matrix(0,400,1))
  stopifnot(got$numerical_ok,got$quad_z_diff<1e-10)
  stopifnot(exists('rti_probability_precision',mode='function'))
  lower<-matrix(pt(-mu[,1],7,log.p=TRUE),400,1)
  upper<-matrix(pt(-mu[,1],7,lower.tail=FALSE,log.p=TRUE),400,1)
  precision<-rti_probability_precision(lower,upper,chains=4L,mcse_limit=.03)
  expected<-as.numeric(posterior::mcse_mean(matrix(exp(lower),100,4)))/
    dnorm(qnorm(mean(exp(lower))))
  near(precision$posterior_mcse_z,expected,1e-9)
  stopifnot(!precision$posterior_mc_ok,precision$posterior_mcse_z>.03,
            precision$posterior_mc_status=='posterior_MC_unresolved')
})

check("posterior precision preserves log tails and rejects undefined MCSE", {
  set.seed(654)
  logsmall<-matrix(-1000+rnorm(4000,sd=.1),4000,1)
  got<-rti_probability_precision(matrix(0,4000,1),logsmall,chains=4L,mcse_limit=.03)
  stopifnot(got$posterior_mc_ok,is.finite(got$posterior_mcse_z),got$posterior_mcse_z>0)
  undefined<-rti_probability_precision(matrix(log(.5),400,1),matrix(log(.5),400,1),4L)
  stopifnot(!undefined$posterior_mc_ok,undefined$posterior_mc_status=='posterior_MC_undefined')
  must_error(rti_probability_precision(matrix(log(.5),403,1),matrix(log(.5),403,1),4L))
})

check("normative release preserves estimates but masks unresolved MC precision", {
  stopifnot(exists('rti_apply_marginal_precision',mode='function'))
  p<-data.frame(z=c(1,2,3,4),numerical_ok=c(TRUE,TRUE,FALSE,TRUE),
    posterior_mc_ok=c(TRUE,FALSE,TRUE,FALSE),posterior_mcse_z=c(.01,.2,.01,NA),
    posterior_mc_status=c('ok','posterior_MC_unresolved','ok','posterior_MC_undefined'))
  got<-rti_apply_marginal_precision(p)
  stopifnot(identical(got$z_estimate,c(1,2,3,4)),identical(is.finite(got$z),c(TRUE,FALSE,FALSE,FALSE)),
    identical(got$numerical_ok,p$numerical_ok),got$score_status[2]=='posterior_MC_unresolved',
    got$score_status[3]=='numerical_unresolved',got$score_status[4]=='posterior_MC_undefined')
})

adapter_fit<-Sys.getenv('RTI_TEST_FIT','')
if(nzchar(adapter_fit)) {
  check('actual brms adapter preserves complete chains and singleton scoring rows', {
    fit<-readRDS(adapter_fit)
    nd<-fit$data[c(1L,6L,14L),,drop=FALSE]
    nd$NM_raw<-expm1(nd$NM_response);nd$row_id_global<-c(1L,6L,14L)
    a<-rti_predict_marginal(fit,nd,chunk_size=2L)
    b<-rti_predict_marginal(fit,nd,chunk_size=1L)
    near(a$z,b$z,1e-12);near(a$posterior_mcse_z,b$posterior_mcse_z,1e-12)
    stopifnot(all(a$posterior_n_chains==posterior::nchains(fit)),
      all(a$posterior_iterations_per_chain==posterior::niterations(fit)))
    n<-posterior::ndraws(fit)
    subset<-rti_predict_marginal(fit,nd,draw_ids=rev(seq_len(n%/%2)),chunk_size=2L)
    stopifnot(all(subset$posterior_mc_status=='not_assessed_draw_subset'),
      !any(subset$posterior_mc_ok),all(is.na(subset$posterior_mcse_z)))
  })
}

if (nzchar(adapter_fit)) {
  fit <- readRDS(adapter_fit)
  nd <- fit$data[1:4, , drop = FALSE]
  nd$NM_raw <- expm1(nd$NM_response)
  nd$row_id_global <- c(17L, 4L, 102L, 9L)
  ids <- 1:200
  check("non-Gaussian or specially weighted grouping structures are refused", {
    bad <- fit
    bad$formula <- brms::bf(NM_response ~ 1 +
      (1 + Time_since_baseline_1y | gr(Subject_ID, dist = "student")),
      sigma ~ Platform, family = brms::student())
    must_error(.rti_validate_fit_structure(bad))
    bad$formula <- brms::bf(NM_response ~ 1 +
      (1 + Time_since_baseline_1y | gr(Subject_ID, cov = A)),
      sigma ~ Platform, family = brms::student())
    must_error(.rti_validate_fit_structure(bad))
    bad$formula <- brms::bf(NM_response ~ 1 +
      (1 + Time_since_baseline_1y | Subject_ID),
      sigma ~ Platform + (1 | Subject_ID), family = brms::student())
    must_error(.rti_validate_fit_structure(bad))
    bad$formula <- brms::bf(NM_response ~ 1 +
      (1 + Time_since_baseline_1y | Subject_ID),
      sigma ~ Platform, nu ~ Platform, family = brms::student())
    must_error(.rti_validate_fit_structure(bad))
  })
  result <- rti_predict_marginal(fit, nd, draw_ids = ids, chunk_size = 3L,
    quad_nodes = 80L, check_nodes = 160L, max_nodes = 1280L,
    z_tolerance = 1e-8, log_score_tolerance = 1e-8)
  check("saved-fit API preserves row order and converges at requested precision", {
    stopifnot(identical(result$row_id_global, nd$row_id_global),
              all(result$numerical_ok), all(result$quad_z_diff < 0.02),
              all(result$quad_log_score_diff <= 0.01))
    near(result$p_lower + result$p_upper, rep(1, 4), 1e-12)
  })
  check("saved-fit score includes log1p Jacobian and integrates random effects", {
    mu <- brms::posterior_linpred(fit, newdata = nd, re_formula = NA, draw_ids = ids)
    sigma <- brms::posterior_linpred(fit, newdata = nd, re_formula = NA,
                                    dpar = "sigma", transform = TRUE, draw_ids = ids)
    nu <- as.numeric(posterior::as_draws_matrix(fit, variable = "nu")) [ids]
    vc <- brms::VarCorr(fit, summary = FALSE)$Subject_ID$cov[ids, , , drop = FALSE]
    variance <- rti_re_variance(vc, nd)
    exact <- rti_cdf_components(nd$NM_response, mu, sigma, nu, variance, 640L)
    near(result$log_score_raw, exact$log_density - log1p(nd$NM_raw), 1e-8)
    near(result$z, exact$z, 1e-8)
    without_re <- rti_cdf_components(nd$NM_response, mu, sigma, nu,
                                     matrix(0, length(ids), nrow(nd)), 80L)
    stopifnot(max(abs(result$z - without_re$z)) > 0.05)
    # Independent analytic-CDF Monte Carlo integrates only the Gaussian RE;
    # this has lower simulation variance than drawing Student residuals as well.
    set.seed(240921)
    mc_n <- 1000000L
    mc <- vapply(seq_len(nrow(nd)), function(i) {
      d <- sample.int(length(ids), mc_n, replace = TRUE)
      b <- rnorm(mc_n, sd = sqrt(variance[d, i]))
      mean(pt((nd$NM_response[i] - mu[d, i] - b) / sigma[d, i], nu[d]))
    }, numeric(1))
    near(result$p_lower, mc, 0.002)
    cat("Saved-fit results (first 4 rows, 200 draws):\n")
    print(result, digits = 12, row.names = FALSE)
    cat("Independent 1,000,000-draw RE integration absolute CDF differences:\n")
    print(abs(result$p_lower - mc), digits = 8)
    cat("Without-RE z (must differ):\n")
    print(without_re$z, digits = 8)
  })
  check("chunk and posterior draw ordering do not change results", {
    reordered <- rti_predict_marginal(fit, nd, draw_ids = rev(ids), chunk_size = 1L,
      quad_nodes = 80L, check_nodes = 160L, max_nodes = 1280L,
      z_tolerance = 1e-8, log_score_tolerance = 1e-8)
    near(result$z, reordered$z, 1e-11)
    near(result$log_score_raw, reordered$log_score_raw, 1e-11)
  })
  check("insufficient quadrature is flagged and invalid raw responses refused", {
    bad_quad <- rti_predict_marginal(fit, nd, ids, quad_nodes = 2L, check_nodes = 4L,
                                    z_tolerance = 1e-12, max_nodes = 4L)
    stopifnot(any(!bad_quad$numerical_ok))
    invalid <- nd
    invalid$NM_raw[1] <- -1
    must_error(rti_predict_marginal(fit, invalid, ids))
  })
}
cat(sprintf("\nAll %d numerical checks passed.\n", passed))
