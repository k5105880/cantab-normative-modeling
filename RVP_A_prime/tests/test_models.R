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
cat("ALL MODEL TESTS PASSED\n")
