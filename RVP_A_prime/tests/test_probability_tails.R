# Real brms prediction adapter with constructed draws; no MCMC is run here.
args <- gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=','',args[grepl('^--file=',args)][1]))))
source(file.path(dirname(root),'RTI_movement','synthetic_input.R'))
source(file.path(root,'R','folds.R'))
source(file.path(root,'R','prepare_data.R'))
source(file.path(root,'R','model_spec.R'))
source(file.path(root,'R','precision_helpers.R'))
source(file.path(root,'R','calibration.R'))
local({
  old_warn <- getOption('warn');options(warn=2)
  on.exit(options(warn=old_warn),add=TRUE)
  # For Beta(1, .01), the upper boundary mass is exactly delta^.01.
  # A tiny boundary width can still carry substantial probability.
  for(kappa in c(1e-10,1e-30,1e-100,1e-300)) {
    log_mass <- .01*(log(kappa)-log1p(2*kappa))
    log_lower <- log(-expm1(log_mass))
    upper <- xbeta_bounds_logit(1,log(100),1.01,kappa)
    lower <- xbeta_bounds_logit(0,-log(100),1.01,kappa)
    stopifnot(abs(upper$log_mass-log_mass)<1e-12,
      abs(upper$log_lower-log_lower)<1e-12,
      abs(lower$log_upper-log_lower)<1e-12,
      abs(lower$log_mass-log_mass)<1e-12)
  }
  work <- tempfile('aprime_tail_test_');dir.create(work)
  on.exit(unlink(work,recursive=TRUE),add=TRUE)
  input <- file.path(work,'input.csv')
  write.csv(synthetic_normative_input(),input,row.names=FALSE)
  prepared <- prepare_data(input,NULL,file.path(work,'prepared'))
  model <- build_model(prepared$hc)
  fit <- brms::brm(model$formula,data=prepared$hc,family=model$family,
    prior=model$prior,empty=TRUE)
  sd <- brms::standata(fit)
  cn <- c(paste0('b_',colnames(sd$X)),paste0('bs_',colnames(sd$Xs)),
    paste0('b_phi_',colnames(sd$X_phi)),paste0('s_sAge_baseline_10y_1[',1:3,']'),
    paste0('s_sTime_since_baseline_1y_1[',1:2,']'),'kappa',
    'sd_Subject_ID__Intercept','sd_Subject_ID__Time_since_baseline_1y',
    'cor_Subject_ID__Intercept__Time_since_baseline_1y')
  loadNamespace('rstan')
  make_fit <- function(phi=45,reflect=FALSE) {
    set.seed(2727)
    samples <- lapply(1:4,function(ch) {
      s <- setNames(lapply(cn,function(nm)rnorm(250,0,.001)),cn)
      s$b_Intercept <- rnorm(250,qlogis(.8),.005)
      s$b_phi_Intercept <- rnorm(250,log(phi),.005)
      s$kappa <- exp(rnorm(250,log(.001),.005))
      s$sd_Subject_ID__Intercept <- exp(rnorm(250,log(.05),.005))
      s$sd_Subject_ID__Time_since_baseline_1y <- exp(rnorm(250,log(.01),.005))
      s$cor_Subject_ID__Intercept__Time_since_baseline_1y <- rnorm(250,.2,.005)
      if(reflect)for(nm in setdiff(grep('^(b_|bs_|s_)',cn,value=TRUE),
        grep('^b_phi_',cn,value=TRUE)))s[[nm]] <- -s[[nm]]
      s
    })
    f <- fit
    f$fit <- new('stanfit',mode=0L,sim=list(samples=samples,chains=4L,
      n_save=rep(250L,4),warmup2=rep(0L,4)))
    f
  }
  nd <- fit$data[1:2,,drop=FALSE];nd$row_id_global <- 1:2;nd$NM_response <- 1
  # These independent log-tail reference values were checked at 32--256 nodes.
  for(phi in c(45,100)) {
    f <- make_fit(phi)
    ans <- tryCatch(score_cdf_quadrature(f,nd,'normative_full'),error=identity)
    if(inherits(ans,'error'))stop('Finite upper-tail score was rejected: ',conditionMessage(ans))
    want <- if(phi==45)8.667864089469923 else 12.60848750772812
    stopifnot(abs(ans$z[1]-want)<1e-8,all(ans$precision_status=='ok'),
      all(is.finite(ans$posterior_mcse_z)),all(ans$posterior_mcse_z<.03))
    if(phi==45)stopifnot(abs(ans$posterior_mcse_z[1]-.0010068727)<1e-9)
    # Reflecting the beta distribution preserves precision and reverses z.
    reflected <- nd;reflected$NM_response <- 0
    left <- score_cdf_quadrature(make_fit(phi,TRUE),reflected,'normative_full',chunk=1L)
    stopifnot(max(abs(left$z+ans$z))<1e-8,all(left$precision_status=='ok'))
    # Calibration must retain the boundary mass even if printed CDF bounds equal 1.
    moment <- pit_moments(ans$p_lower,ans$p_upper,ans$log_lower,ans$log_mass,ans$log_upper)
    a <- qnorm(ans$log_mass,lower.tail=FALSE,log.p=TRUE)
    mills <- exp(dnorm(a,log=TRUE)-ans$log_mass)
    stopifnot(max(abs(moment$m-mills))<1e-10,
      max(abs(moment$m2-(1+a*mills)))<1e-9,all(moment$low==0),all(moment$high==1))
    mirror <- pit_moments(left$p_lower,left$p_upper,left$log_lower,left$log_mass,left$log_upper)
    stopifnot(max(abs(moment$m+mirror$m))<1e-8,max(abs(moment$m2-mirror$m2))<1e-7)
  }
  # Carry extreme endpoint probabilities through the exported CSV and bootstrap.
  scores <- rbind(ans,left);scores$row_id_global <- 1:4
  scores$role <- 'normative_OOF_fold1'
  results <- file.path(work,'results');dir.create(results)
  write.csv(scores,file.path(results,'normative_score_precision.csv'),row.names=FALSE)
  write.csv(data.frame(row_id_global=1:4,Subject_ID=c('001','002','003','004'),
    Diagnosis='HC',Platform='RS6',Practice_bin='0_baseline',Time_since_baseline_years=0),
    file.path(results,'RVP_A_prime_scores.csv'),row.names=FALSE)
  calibrate_scores(work,list(seed=2727L,calibration_bootstraps=32L))
  saved <- read.csv(file.path(results,'OOF_PIT_moments.csv'))
  want_moments <- rbind(moment,mirror)
  stopifnot(max(abs(as.matrix(saved[names(want_moments)])-as.matrix(want_moments)))<1e-9)
  summary <- read.csv(file.path(results,'OOF_calibration.csv'))
  stopifnot(all(is.finite(as.matrix(summary[vapply(summary,is.numeric,logical(1))]))),
    all(summary$n_subjects==4),all(summary$low_tail==.5),all(summary$high_tail==.5))
  # Ordinary and very narrow PIT intervals agree with adaptive integration.
  a <- c(0,.01,.2,.9,.999);b <- c(.01,.04,.2+1e-12,.93,1)
  moments <- pit_moments(a,b)
  for(i in seq_along(a))for(power in 1:2) {
    reference <- integrate(function(u)qnorm(a[i]+(b[i]-a[i])*u)^power,0,1,
      rel.tol=1e-9)$value
    stopifnot(abs(moments[i,power]-reference)<1e-7)
  }
  # Interior probabilities and tiny boundary widths agree with adaptive integration.
  for(tiny_boundary in c(FALSE,TRUE)) {
    f <- make_fit(if(tiny_boundary).01 else 45)
    if(tiny_boundary)for(ch in 1:4)f$fit@sim$samples[[ch]]$kappa <- rep(1e-30,250)
    nd$NM_response <- if(tiny_boundary)c(1,1) else c(.8,.99)
    ans <- score_cdf_quadrature(f,nd,'normative_full')
    stopifnot(all(ans$precision_status=='ok'))
    dr <- posterior::as_draws_matrix(f)
    pr <- brms::prepare_predictions(f,newdata=nd,re_formula=NA)
    for(j in 1:2) {
      eta <- brms::get_dpar(pr,'mu',i=j,inv_link=FALSE)
      phi <- brms::get_dpar(pr,'phi',i=j);kappa <- brms::get_dpar(pr,'kappa',i=j)
      t <- nd$Time_since_baseline_1y[j]
      re_sd <- sqrt(dr[,'sd_Subject_ID__Intercept']^2+t^2*dr[,'sd_Subject_ID__Time_since_baseline_1y']^2+
        2*t*dr[,'cor_Subject_ID__Intercept__Time_since_baseline_1y']*
        dr[,'sd_Subject_ID__Intercept']*dr[,'sd_Subject_ID__Time_since_baseline_1y'])
      integrand <- function(u)vapply(u,function(v) {
        if(tiny_boundary)tail <- pbeta(kappa/(1+2*kappa),
          plogis(-eta-re_sd*v)*phi,plogis(eta+re_sd*v)*phi)/2 else
          tail <- pbeta((nd$NM_response[j]+kappa)/(1+2*kappa),
            plogis(eta+re_sd*v)*phi,plogis(-eta-re_sd*v)*phi,lower.tail=FALSE)
        mean(tail)*dnorm(v)
      },numeric(1))
      upper <- integrate(integrand,-12,12,rel.tol=1e-10,abs.tol=1e-100)$value
      stopifnot(abs(ans$z[j]-qnorm(upper,lower.tail=FALSE))<1e-7)
    }
  }
  # Saturating an inverse-logit must not destroy a finite predictor before scoring.
  f <- make_fit()
  for(ch in 1:4)f$fit@sim$samples[[ch]]$b_Intercept <- rep(40,250)
  nd$NM_response <- 1
  ans <- score_cdf_quadrature(f,nd,'normative_full')
  stopifnot(all(is.finite(ans$z)),all(ans$precision_status=='ok'))
  cat('PASS: A-prime extreme tails, reflected endpoints, MCSE, calibration and actual brms adapters.\n')
})
