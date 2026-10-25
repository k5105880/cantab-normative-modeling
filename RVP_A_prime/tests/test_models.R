#!/usr/bin/env Rscript
args <- gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]))))
check <- function(x, label) {
  if (!isTRUE(x)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS:", label, "\n")
}
check(file.exists(file.path(root, "R/model_spec.R")), "model library exists")
source(file.path(root, "R/model_spec.R"))
source(file.path(root, "R/fit_models.R"))
suppressPackageStartupMessages(library(brms))
set.seed(2801)
d <- data.frame(Subject_ID=factor(rep(seq_len(40),each=4)),
  Age_baseline_10y=rep(seq(-2,2,length.out=40),each=4),
  Time_since_baseline_1y=rep(c(0,1,2,3),40)+runif(160,0,.9),
  FSIQ_15=rnorm(160), Gender=factor(rep(c("F","M"),each=4,length.out=160)),
  Practice_bin=factor(rep(c("0_baseline","1_second","2_third","3_fourth_or_later"),40)),
  Platform=factor(sample(c("Connect","RS6"),160,TRUE)),
  NM_response=runif(160,.65,.99),row_id_global=seq_len(160),Diagnosis="HC")
d$PlatformConnect_indicator <- as.numeric(d$Platform=="Connect")
m <- build_model(d)
v <- all.vars(m$formula$formula)
check(identical(m$family$family,"xbeta"), "censored extended beta family")
check("FSIQ_15" %in% v,"FSIQ predictor included")
s <- brms::make_standata(m$formula,data=d,family=m$family,prior=m$prior)
check(qr(cbind(s$X,s$Xs))$rank==ncol(s$X)+ncol(s$Xs),"full rank design")
check(s$M_1==2 && s$NC_1==1,"correlated intercept and slope")
good <- data.frame(n_divergent=0,n_max_treedepth=0,max_rhat=1.005,
  min_ess_bulk=900,min_ess_tail=1000,min_EBFMI=.5)
check(isTRUE(accept_sampler_diagnostics(good)),"adequate sampler accepted")
for (nm in names(good)) {
  bad <- good; bad[[nm]] <- switch(nm,n_divergent=1,n_max_treedepth=1,
    max_rhat=1.01,min_ess_bulk=399,min_ess_tail=399,min_EBFMI=.3)
  check(inherits(try(accept_sampler_diagnostics(bad),silent=TRUE),"try-error"),paste(nm,"failure rejected"))
}
check(inherits(try(accept_draw_shape(999L,4L,list(iter=2500L,warmup=1500L,chains=4L)),silent=TRUE),"try-error"),
      "incomplete chains rejected")
f <- tempfile(fileext=".rds")
saveRDS(list(identity="wrong",fit=list(),fit_sha256=model_hash(list()),
  metadata=list(),metadata_sha256=model_hash(list()),diagnostics=good,accepted=TRUE),f)
check(inherits(try(read_accepted_fit(f,"expected"),silent=TRUE),"try-error"),"cache identity mismatch rejected")
bad <- list(identity="expected",fit=list(),fit_sha256="invalid",metadata=list(),
  metadata_sha256=model_hash(list()),diagnostics=good,accepted=TRUE)
saveRDS(bad,f)
before <- digest::digest(file=f,algo="sha256")
check(inherits(try(read_accepted_fit(f,"expected"),silent=TRUE),"try-error"),"corrupt fit hash rejected")
check(identical(before,digest::digest(file=f,algo="sha256")),"cache rejection is read-only")
unlink(f)
source(file.path(root,"R/prior_direct.R"))
# The prior simulation must accept every identifiable design allowed by build_model,
# including data-dependent linear terms and a constant IQ predictor.
for (scenario in c("full", "linear_time", "linear_age", "linear_both", "constant_iq", "intercept_phi")) {
  x <- d
  if (scenario %in% c("linear_time", "linear_both")) {
    x <- droplevels(x[rep(c(TRUE,TRUE,TRUE,FALSE),40),])
    x$Time_since_baseline_1y <- rep(c(0,1,2),40)
    x$Time_since_baseline_1y[seq(6,120,6)] <- 3
  }
  if (scenario %in% c("linear_age", "linear_both"))
    x$Age_baseline_10y <- (as.integer(x$Subject_ID) %% 4)-2
  if (scenario=="constant_iq") x$FSIQ_15 <- 0
  if (scenario=="intercept_phi") {
    x$Platform <- factor("Connect"); x$Practice_bin <- factor("0_baseline")
    x$PlatformConnect_indicator <- 1
  }
  model <- build_model(x)
  model_design_audit(model,x)
  # Generate the real brms model and Stan data without compilation or sampling.
  fit <- brms::brm(model$formula,data=x,prior=model$prior,family=model$family,empty=TRUE)
  directory <- tempfile(paste0("prior_",scenario,"_"))
  result <- suppressWarnings(run_prior_direct(fit,directory,draws=2000L,batch_size=500L,training_data=x))
  check(all(is.finite(result$rates$mean)) && all(is.finite(result$quantiles$prior_predictive_quantile)),
        paste("complete direct prior simulation:",scenario))
  check(setequal(unique(result$rates$platform),c("All",as.character(unique(x$Platform)))),
        paste("platform summaries preserved:",scenario))
  wrong_rows <- x[rev(seq_len(nrow(x))),]
  check(inherits(try(run_prior_direct(fit,directory,draws=2000L,batch_size=500L,
                                     training_data=wrong_rows),silent=TRUE),"try-error"),
        paste("mismatched prior-summary rows rejected:",scenario))
  # A changed prior must still be rejected, rather than silently simulated with
  # the original primitive distribution.
  changed <- fit
  changed$model <- sub("normal_lpdf(b | 0, 0.5)","normal_lpdf(b | 0, 2)",fit$model,fixed=TRUE)
  check(inherits(try(validate_prior_model(changed),silent=TRUE),"try-error"),
        paste("changed prior rejected:",scenario))
  unlink(directory,recursive=TRUE)
}
cat("ALL MODEL TESTS PASSED\n")
