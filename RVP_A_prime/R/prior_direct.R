#!/usr/bin/env Rscript
# Independent draws from the exact primitive priors of the normative model.
# No HMC, posterior refit, rejection of extremes, or prior change is performed.

validate_prior_model <- function(fit) {
  code <- fit$model
  d <- brms::make_standata(fit)
  splines <- as.integer(sub("^nb_", "", grep("^nb_[0-9]+$",names(d),value=TRUE)))
  stopifnot(length(splines)<=2L, identical(splines,seq_along(splines)))
  d$n_splines <- length(splines)
  if(!length(splines)) d$Ks <- 0L
  compact <- function(x) gsub("[[:space:]]+", "", x)
  lines <- trimws(strsplit(code, "\n", fixed = TRUE)[[1]])
  expected <- c(
    "lprior += normal_lpdf(b | 0, 0.5);",
    "lprior += student_t_lpdf(Intercept | 3, 0, 2.5);",
    if(length(splines)) "lprior += normal_lpdf(bs | 0, 0.5);",
    sprintf("lprior += exponential_lpdf(sds_%d | 1);",splines),
    if(d$Kc_phi>0L) "lprior += normal_lpdf(b_phi | 0, 0.4);",
    "lprior += student_t_lpdf(Intercept_phi | 3, 0, 2.5);",
    "lprior += gamma_lpdf(kappa | 0.01, 0.01);",
    "lprior += exponential_lpdf(sd_1 | 1);",
    "lprior += lkj_corr_cholesky_lpdf(L_1 | 2);")
  actual <- lines[grepl("^lprior[[:space:]]*\\+=", lines)]
  if (!identical(sort(compact(actual)), sort(compact(expected))))
    stop("Stan prior assertion failed: unsupported primitive prior density")
  targets <- c("target += xbeta_lpdf(Y | mu, phi, kappa);", "target += lprior;",
    sprintf("target += std_normal_lpdf(zs_%d_1);",splines),
    "target += std_normal_lpdf(to_vector(z_1));")
  if (!identical(sort(compact(lines[grepl("^target[[:space:]]*\\+=", lines)])),
                 sort(compact(targets))))
    stop("Stan target assertion failed: unsupported additional density")
  mu_terms <- c("Intercept", "Xc * b", if(length(splines)) "Xs * bs",
                sprintf("Zs_%d_1 * s_%d_1",splines,splines))
  required <- c(
    "Xc[, i - 1] = X[, i] - means_X[i - 1];",
    if(d$Kc_phi>0L) "Xc_phi[, i - 1] = X_phi[, i] - means_X_phi[i - 1];",
    sprintf("s_%d_1 = sds_%d[1] * zs_%d_1;",splines,splines,splines),
    "return transpose(diag_pre_multiply(SD, L) * z);",
    "r_1 = scale_r_cor(z_1, sd_1, L_1);",
    paste0("mu += ",paste(mu_terms,collapse=" + "),";"),
    if(d$Kc_phi>0L) "phi += Intercept_phi + Xc_phi * b_phi;" else "phi += Intercept_phi;",
    "mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];",
    "mu = inv_logit(mu);", "phi = exp(phi);")
  if (!all(vapply(compact(required), grepl, logical(1), x=compact(code), fixed=TRUE)))
    stop("Stan transformation assertion failed: unsupported linear predictor")
  stopifnot(d$M_1 == 2, d$NC_1 == 1,
    all(vapply(splines,function(i) identical(d[[paste0("nb_",i)]],1L),logical(1))),
    d$Kc == d$K-1L, d$Kc_phi == d$K_phi-1L, d$N == nrow(fit$data),
    all(d$X[,1] == 1), all(d$X_phi[,1] == 1),
    identical(as.numeric(d$Y), as.numeric(fit$data$NM_response)),
    identical(as.integer(d$J_1), as.integer(factor(fit$data$Subject_ID))))
  d$Xc <- sweep(d$X[,-1,drop=FALSE], 2, colMeans(d$X[,-1,drop=FALSE]), "-")
  d$Xc_phi <- if(d$Kc_phi>0L) sweep(d$X_phi[,-1,drop=FALSE], 2,
                    colMeans(d$X_phi[,-1,drop=FALSE]), "-") else matrix(0,d$N,0)
  audit <- data.frame(check=c("exact_prior_density_assertions", "exact_target_density_assertions",
    "exact_predictor_assertions",
    "subject_design_index_matches", "centered_X_max_abs_column_mean",
    "centered_X_phi_max_abs_column_mean", "n_observations", "n_subjects"),
    value=c("PASS", "PASS", "PASS",
      "PASS", max(abs(colMeans(d$Xc))), if(d$Kc_phi>0L) max(abs(colMeans(d$Xc_phi))) else 0, d$N, d$N_1))
  list(data=d, audit=audit)
}

draw_prior_parameters <- function(n, d) {
  normal_matrix <- function(k, sd=1) matrix(rnorm(n*k,sd=sd), n, k)
  p <- list(b=normal_matrix(d$Kc,.5), bs=normal_matrix(d$Ks,.5),
    b_phi=normal_matrix(d$Kc_phi,.4), Intercept=2.5*rt(n,3),
    Intercept_phi=2.5*rt(n,3),
    # Gamma boost identity is exact and keeps tiny positive kappa in log form.
    log_kappa=log(rgamma(n,shape=1.01,rate=.01))+log(runif(n))/.01,
    sds1=if(d$n_splines>=1L) rexp(n), sds2=if(d$n_splines>=2L) rexp(n), sd1=rexp(n), sd2=rexp(n),
    rho=2*rbeta(n,2,2)-1,
    zs1=if(d$n_splines>=1L) normal_matrix(d$knots_1[1]),
    zs2=if(d$n_splines>=2L) normal_matrix(d$knots_2[1]))
  z1 <- normal_matrix(d$N_1)
  z2 <- normal_matrix(d$N_1)
  p$r1 <- z1*p$sd1
  p$r2 <- (z1*p$rho + z2*sqrt(1-p$rho^2))*p$sd2
  p
}

direct_linear_predictors <- function(p, d) {
  eta <- p$b %*% t(d$Xc) + p$Intercept
  if(d$n_splines>=1L) eta <- eta + p$bs %*% t(d$Xs) + (p$zs1*p$sds1) %*% t(d$Zs_1_1)
  if(d$n_splines>=2L) eta <- eta + (p$zs2*p$sds2) %*% t(d$Zs_2_1)
  eta <- eta + sweep(p$r1[,d$J_1,drop=FALSE],2,d$Z_1_1,"*") +
    sweep(p$r2[,d$J_1,drop=FALSE],2,d$Z_1_2,"*")
  list(eta=eta, log_phi=p$b_phi %*% t(d$Xc_phi)+p$Intercept_phi)
}

softplus <- function(x) pmax(x,0)+log1p(exp(-abs(x)))

beta_log_draw <- function(eta, log_phi) {
  # Avoid 1-plogis(eta): it loses the second shape at eta around 37.
  eta <- as.vector(eta); log_phi <- as.vector(log_phi)
  stopifnot(length(eta)==length(log_phi), all(is.finite(eta)), all(is.finite(log_phi)))
  n <- length(eta)
  lm <- -softplus(-eta); lmc <- -softplus(eta)
  la <- log_phi+lm; lb <- log_phi+lmc
  lt <- l1t <- rep(NA_real_,n)
  approximation <- rep(FALSE,n)
  # With both shapes >=1 and moderate mean/precision, base R's beta RNG is
  # accurate without latent endpoint saturation. Other cases use the exact
  # beta = G_a/(G_a+G_b) identity, preserving both logarithmic tails.
  regular <- la>=0 & lb>=0 & pmax(la,lb)<=log(1e6) & abs(eta)<=4
  if (any(regular)) {
    t <- rbeta(sum(regular),exp(la[regular]),exp(lb[regular]))
    lt[regular] <- log(t); l1t[regular] <- log1p(-t)
    if (any(t==0 | t==1)) stop("Unexpected endpoint from moderate-shape beta RNG")
  }
  use_gamma <- !regular
  if (any(use_gamma)) {
    log_gamma <- function(ls) {
      out <- numeric(length(ls))
      tiny <- ls < -700
      # log(G) magnitude would exceed roughly 1e304; the zero-shape limit
      # is explicitly flagged by the caller. No such draws are silently lost.
      out[tiny] <- -Inf
      huge <- ls > 700
      # Relative SD is <exp(-350); retain the logarithmic mean. Applying
      # this separately to each gamma also handles a tiny opposite shape.
      out[huge] <- ls[huge]
      small <- !tiny & ls < 0
      if(any(small)) out[small] <- log(rgamma(sum(small),exp(ls[small])+1)) +
        log(runif(sum(small)))/exp(ls[small])
      large <- !tiny & !small & !huge
      if(any(large)) out[large] <- log(rgamma(sum(large),exp(ls[large])))
      out
    }
    ga <- log_gamma(la[use_gamma]); gb <- log_gamma(lb[use_gamma])
    both_tiny <- is.infinite(ga) & is.infinite(gb)
    delta <- gb-ga
    tl <- -softplus(delta); t1l <- -softplus(-delta)
    if(any(both_tiny)) {
      upper <- runif(sum(both_tiny)) < plogis(eta[use_gamma][both_tiny])
      tl[both_tiny] <- ifelse(upper,0,-Inf)
      t1l[both_tiny] <- ifelse(upper,-Inf,0)
    }
    lt[use_gamma] <- tl; l1t[use_gamma] <- t1l
    approximation[use_gamma] <- abs(la[use_gamma]) > 700 | abs(lb[use_gamma]) > 700
  }
  stopifnot(all(!is.na(lt)), all(!is.na(l1t)), all(lt<=0), all(l1t<=0))
  list(log_t=lt,log1m_t=l1t,approximation=approximation,
    used_log_gamma=use_gamma, latent_numeric_zero=exp(lt)==0,
    latent_numeric_one=exp(lt)==1,
    shape_exp_underflow=(exp(la)==0 | exp(lb)==0),
    shape_exp_overflow=(la>log(.Machine$double.xmax) | lb>log(.Machine$double.xmax)),
    mu_numeric_endpoint=(plogis(eta)==0 | plogis(eta)==1))
}

censor_xbeta <- function(log_t,log1m_t,log_kappa) {
  stopifnot(length(log_t)==length(log1m_t), length(log_t)==length(log_kappa),
    all(is.finite(log_kappa)))
  ld <- softplus(log(2)+log_kappa)
  threshold <- log_kappa-ld
  floor <- log_t <= threshold
  ceiling <- log1m_t <= threshold
  stopifnot(!any(floor & ceiling))
  y <- numeric(length(log_t)); y[ceiling] <- 1
  interior <- !floor & !ceiling
  lower_half <- interior & log_t <= -log(2)
  if(any(lower_half)) {
    a <- log_t[lower_half]+ld[lower_half]
    y[lower_half] <- exp(a)*(-expm1(log_kappa[lower_half]-a))
  }
  upper_half <- interior & !lower_half
  if(any(upper_half)) {
    a <- log1m_t[upper_half]+ld[upper_half]
    y[upper_half] <- 1-exp(a)*(-expm1(log_kappa[upper_half]-a))
  }
  stopifnot(all(is.finite(y)), all(y>=0 & y<=1))
  list(y=y,floor=floor,ceiling=ceiling,numeric_zero=interior & y==0,
       numeric_one=interior & y==1)
}

quantile_mc_interval <- function(independent_y, probs) {
  # Each y must come from a DIFFERENT independent parameter/dataset draw.
  # Inverting Binomial(n,p) gives conservative 95% quantile confidence bounds
  # without a density estimate or a normal approximation near boundary atoms.
  y <- sort(independent_y); n <- length(y)
  lower_order <- qbinom(.025,n,probs)
  upper_order <- qbinom(.975,n,probs)+1L
  data.frame(probability=probs,
    lower=ifelse(lower_order<1L,0,y[pmax(1L,lower_order)]),
    upper=ifelse(upper_order>n,1,y[pmin(n,upper_order)]),
    lower_order=lower_order,upper_order=upper_order)
}

run_prior_direct <- function(fit_or_path, outdir, draws=40000L, seed=20260921L, batch_size=500L,
                             training_data=NULL) {
  suppressPackageStartupMessages(library(brms))
  out_dir <- outdir
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  fit <- if(is.character(fit_or_path)) read_accepted_fit(fit_or_path) else fit_or_path
  if(!inherits(fit,"brmsfit")) stop("Prior audit requires the newly fitted primary brmsfit")
  if("NM_baseline_response" %in% all.vars(fit$formula$formula))
    stop("Direct prior audit is for the primary population specification")
  # brms drops columns unused by the formula. Retain platform groups when a
  # single represented platform has correctly been omitted from the model.
  if(is.null(training_data)) training_data <- fit$data
  if(!is.data.frame(training_data) || nrow(training_data)!=nrow(fit$data) ||
     !all(names(fit$data)%in%names(training_data)) ||
     !all(vapply(names(fit$data),function(nm)
       identical(as.character(fit$data[[nm]]),as.character(training_data[[nm]])),logical(1))))
    stop("Prior summary training rows do not match the fitted model")
  if(!"Platform"%in%names(training_data) || anyNA(training_data$Platform))
    stop("Training data with Platform are required for prior summaries")
  validation <- validate_prior_model(fit); d <- validation$data
  stopifnot(draws>=2000, draws%%batch_size==0)
  write.csv(validation$audit,file.path(out_dir,"design_and_prior_validation.csv"),row.names=FALSE)
  writeLines(fit$model,file.path(out_dir,"validated_model.stan"))
  platform <- as.character(training_data$Platform)
  groups <- c(list(All=seq_len(d$N)),split(seq_len(d$N),platform))
  probs <- c(.01,.025,.05,.25,.5,.75,.95,.975,.99)
  all_y <- matrix(NA_real_,draws,d$N)
  rates <- list(); batch_q <- list(); diagnostics <- list(); parameter_audit <- list()
  set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
  for (b in seq_len(draws/batch_size)) {
    ii <- ((b-1L)*batch_size+1L):(b*batch_size)
    p <- draw_prior_parameters(batch_size,d)
    lp <- direct_linear_predictors(p,d)
    z <- beta_log_draw(lp$eta,lp$log_phi)
    x <- censor_xbeta(z$log_t,z$log1m_t,rep(p$log_kappa,d$N))
    yy <- matrix(x$y,batch_size,d$N)
    ff <- matrix(x$floor,batch_size,d$N); cc <- matrix(x$ceiling,batch_size,d$N)
    nz <- matrix(x$numeric_zero,batch_size,d$N); no <- matrix(x$numeric_one,batch_size,d$N)
    all_y[ii,] <- yy
    for (g in names(groups)) {
      j <- groups[[g]]; obs <- d$Y[j]
      rates[[length(rates)+1L]] <- data.frame(draw=ii,platform=g,
        mean=rowMeans(yy[,j,drop=FALSE]), floor=rowMeans(ff[,j,drop=FALSE]),
        ceiling=rowMeans(cc[,j,drop=FALSE]),
        below_0_5=rowMeans(yy[,j,drop=FALSE]<.5),
        below_observed_min=rowMeans(yy[,j,drop=FALSE]<min(obs)),
        below_observed_p025=rowMeans(yy[,j,drop=FALSE]<quantile(obs,.025,names=FALSE)),
        numeric_zero=rowMeans(nz[,j,drop=FALSE]), numeric_one=rowMeans(no[,j,drop=FALSE]))
      batch_q[[length(batch_q)+1L]] <- data.frame(batch=b,platform=g,probability=probs,
        quantile=as.numeric(quantile(yy[,j,drop=FALSE],probs)))
    }
    diagnostic <- c(observations=length(x$y), independent_parameter_draws=batch_size,
      kappa_exp_underflow=sum(exp(p$log_kappa)==0),
      kappa_subnormal_or_underflow=sum(p$log_kappa<log(.Machine$double.xmin)),
      log_gamma_beta_draws=sum(z$used_log_gamma), approximated_beta_draws=sum(z$approximation),
      beta_shape_exp_underflow=sum(z$shape_exp_underflow),
      beta_shape_exp_overflow=sum(z$shape_exp_overflow),
      naive_mu_endpoint=sum(z$mu_numeric_endpoint),
      naive_latent_beta_zero=sum(z$latent_numeric_zero),
      naive_latent_beta_one=sum(z$latent_numeric_one),
      censored_floor=sum(x$floor),censored_ceiling=sum(x$ceiling),
      interior_numeric_zero=sum(x$numeric_zero), interior_numeric_one=sum(x$numeric_one))
    diagnostics[[b]] <- data.frame(batch=b,metric=names(diagnostic),count=as.numeric(diagnostic))
    audit_columns <- list(draw=ii,Intercept=p$Intercept,Intercept_phi=p$Intercept_phi,
      b_first=p$b[,1],bs_first=if(ncol(p$bs)) p$bs[,1],
      b_phi_first=if(ncol(p$b_phi)) p$b_phi[,1],
      log_kappa=p$log_kappa,kappa=exp(p$log_kappa),sds1=p$sds1,sds2=p$sds2,
      sd1=p$sd1,sd2=p$sd2,rho=p$rho,
      zs1_first=if(d$n_splines>=1L) p$zs1[,1],zs2_first=if(d$n_splines>=2L) p$zs2[,1])
    parameter_audit[[b]] <- as.data.frame(audit_columns[lengths(audit_columns)>0L])
    if(b%%10==0) cat("Independent prior draws:",max(ii),"/",draws,"\n")
  }
  rates <- do.call(rbind,rates); batch_q <- do.call(rbind,batch_q)
  diagnostics <- do.call(rbind,diagnostics); parameter_audit <- do.call(rbind,parameter_audit)
  diag_summary <- aggregate(count~metric,diagnostics,sum)
  write.csv(diag_summary,file.path(out_dir,"numerical_diagnostics.csv"),row.names=FALSE)
  write.csv(diagnostics,file.path(out_dir,"numerical_diagnostics_by_batch.csv"),row.names=FALSE)
  saveRDS(rates,file.path(out_dir,"independent_draw_metrics.rds"),compress="gzip")
  saveRDS(parameter_audit,file.path(out_dir,"primitive_parameter_audit_draws.rds"),compress="gzip")
  write.csv(batch_q,file.path(out_dir,"independent_batch_quantiles.csv"),row.names=FALSE)
  metrics <- setdiff(names(rates),c("draw","platform"))
  rate_summary <- do.call(rbind,lapply(names(groups),function(g) {
    rows <- rates[rates$platform==g,,drop=FALSE]
    do.call(rbind,lapply(metrics,function(m) {
      v <- rows[[m]]; se <- sd(v)/sqrt(length(v)); q <- quantile(v,c(.025,.5,.975))
      data.frame(platform=g,metric=m,mean=mean(v),MCSE=se,
        MC95_lower=max(0,mean(v)-1.96*se),MC95_upper=min(1,mean(v)+1.96*se),
        prior_dataset_p025=q[1],prior_dataset_p50=q[2],prior_dataset_p975=q[3],
        independent_draws=length(v),observations_per_draw=length(groups[[g]]),row.names=NULL)
    }))
  }))
  marginal <- list()
  quantiles <- do.call(rbind,lapply(names(groups),function(g) {
    j <- groups[[g]]
    pooled <- as.numeric(quantile(all_y[,j,drop=FALSE],probs))
    se <- vapply(probs,function(pr) sd(batch_q$quantile[batch_q$platform==g &
                         batch_q$probability==pr])/sqrt(draws/batch_size),0.0)
    # Uniformly choose one eligible observation from each independent dataset;
    # these are iid draws from the relevant platform mixture prior predictive.
    selected <- sample(j,draws,replace=TRUE)
    yy <- all_y[cbind(seq_len(draws),selected)]
    ci <- quantile_mc_interval(yy,probs)
    marginal[[g]] <<- data.frame(draw=seq_len(draws),platform=g,observation=selected,y=yy)
    data.frame(platform=g,probability=probs,prior_predictive_quantile=pooled,
      independent_marginal_quantile=as.numeric(quantile(yy,probs)),
      MC95_lower=ci$lower,MC95_upper=ci$upper,
      independent_marginal_draws=draws,batch_MCSE=se,
      observed_quantile=as.numeric(quantile(d$Y[j],probs)),independent_batches=draws/batch_size)
  }))
  saveRDS(do.call(rbind,marginal),file.path(out_dir,"independent_marginal_outcomes.rds"),compress="gzip")
  observed <- do.call(rbind,lapply(names(groups),function(g) {
    y <- d$Y[groups[[g]]]
    data.frame(platform=g,n=length(y),min=min(y),mean=mean(y),sd=sd(y),median=median(y),
      max=max(y),floor=mean(y==0),ceiling=mean(y==1),below_0_5=mean(y<.5))
  }))
  parameter_summary <- do.call(rbind,lapply(setdiff(names(parameter_audit),"draw"),function(nm) {
    v <- parameter_audit[[nm]]
    data.frame(parameter=nm,mean=mean(v),sd=sd(v),p001=quantile(v,.001),p01=quantile(v,.01),
      p025=quantile(v,.025),p50=median(v),p975=quantile(v,.975),p99=quantile(v,.99),
      p999=quantile(v,.999),min=min(v),max=max(v),row.names=NULL)
  }))
  write.csv(rate_summary,file.path(out_dir,"prior_predictive_rates.csv"),row.names=FALSE)
  write.csv(quantiles,file.path(out_dir,"prior_predictive_quantiles.csv"),row.names=FALSE)
  write.csv(observed,file.path(out_dir,"observed_comparison.csv"),row.names=FALSE)
  write.csv(parameter_summary,file.path(out_dir,"primitive_parameter_summary.csv"),row.names=FALSE)
  snapshot <- list(seed=seed,draws=draws,batch_size=batch_size,N=d$N,subjects=d$N_1,
    method="independent_exact_primitive_prior_draws",
    source_fit_sha256=digest::digest(fit,algo="sha256"),
    source_model_sha256=digest::digest(fit$model,algo="sha256"),
    sessionInfo=capture.output(sessionInfo()))
  saveRDS(snapshot,file.path(out_dir,"reproducibility_manifest.rds"))
  capture.output(str(snapshot),file=file.path(out_dir,"reproducibility_manifest.txt"))
  count <- function(m) diag_summary$count[diag_summary$metric==m]
  allmean <- function(m) rate_summary$mean[rate_summary$platform=="All" & rate_summary$metric==m]
  notes <- c("Direct prior predictive validation of the primary RVP A-prime model",
    sprintf("%s independent parameter draws, seed %s; %s rows and %s shared subjects per dataset.",
      draws,seed,d$N,d$N_1),
    "Stan primitive densities, predictor transformations, and the fitted make_standata design asserted.",
    "X and X_phi use Stan column centering; Xs and spline bases are used unchanged. Subject intercept/slope draws are shared across all visits.",
    "b and bs: N(0,.5); b_phi: N(0,.4); centered intercepts: t_3(0,2.5); smooth and subject SDs: Exp(1).",
    "kappa: Gamma(.01,rate=.01), sampled using the exact log-gamma boost identity. LKJ(2) in dimension two: rho=2*Beta(2,2)-1.",
    "Moderate beta shapes use base R rbeta; extremes use the exact independent log-gamma ratio identity, retaining log(T) and log(1-T).",
    "Censoring is performed on logarithmic thresholds before converting outcomes to doubles. Numeric interior values rounding to zero/one are not counted as modeled boundary atoms.",
    "Independent sampling has no HMC divergence/Rhat diagnostic.",
    sprintf("kappa double underflows: %s/%s; retained as finite log-kappa, so no draw is omitted.",
      count("kappa_exp_underflow"),draws),
    sprintf("Explicit beta limiting approximations: %s/%s observations (flagged if either log-shape is outside [-700,700]).",
      count("approximated_beta_draws"),draws*d$N),
    sprintf("Interior outcome double rounding: %s to zero; %s to one. Log-domain censoring distinguishes these from atoms.",
      count("interior_numeric_zero"),count("interior_numeric_one")),
    "Independent datasets are the Monte Carlo units for rate MCSE and 95% MC intervals; observations within a dataset are dependent.",
    "Dataset 2.5/50/97.5 percentiles describe prior predictive variation, not uncertainty in the prior mean.",
    "Pooled outcome quantile points use all draws/rows. Quantile MC95 bounds invert the exact binomial distribution for a separate marginal sample with one uniformly selected eligible row per independent dataset; these bounds allow atoms. The saved independent marginal draws make this check reproducible. Batch quantile MCSE is an additional approximate diagnostic.",
    "Outcome quantiles are on the double-precision A-prime scale: a reported 0/1 quantile can reflect numeric rounding. Only the log-domain floor/ceiling indicators estimate the modeled boundary atoms.",
    sprintf("Observed A-prime spans %.6f to %.6f, with mean %.6f. Prior predictive mean %.6f; %.2f%% below .5; %.2f%% below the observed overall minimum.",
      min(d$Y),max(d$Y),mean(d$Y),allmean("mean"),100*allmean("below_0_5"),100*allmean("below_observed_min")),
    sprintf("Modeled prior floor and ceiling rates are %.2f%% and %.2f%%; observed rates are %.2f%% and %.2f%%.",
      100*allmean("floor"),100*allmean("ceiling"),100*mean(d$Y==0),100*mean(d$Y==1)),
    "These proper priors are extremely broad relative to the observed high A-prime scores and permit substantial low/boundary outcomes. Prior simulation does not establish physiological plausibility or justify the priors scientifically.",
    "The direct audit evaluates the specified primary prior without posterior HMC for prior draws.")
  writeLines(notes,file.path(out_dir,"README.txt"))
  if(count("kappa_exp_underflow")+count("interior_numeric_zero")+count("interior_numeric_one")>0)
    warning("Representational saturation detected and audited; log-domain censoring retained all draws. See numerical_diagnostics.csv.")
  if(count("approximated_beta_draws")>0)
    warning("Explicit extreme-shape beta limiting approximations were required. See numerical_diagnostics.csv.")
  cat(paste(notes,collapse="\n"),"\n")
  invisible(list(rates=rate_summary,quantiles=quantiles,diagnostics=diag_summary))
}
