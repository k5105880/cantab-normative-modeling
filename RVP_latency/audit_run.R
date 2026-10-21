# Recompute numerical acceptance and verify row identities in completed model artifacts.
audit_normative_run <- function(run,out,input) {
  out <- cantab_output_directory(out,dirname(script_dir))
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  raw <- read.csv(input,check.names=FALSE,colClasses=c(Subject_ID='character'),fileEncoding='UTF-8-BOM',na.strings=c('','NA','NaN'),stringsAsFactors=FALSE)
  d <- read.csv(file.path(run,'RVP_latency_scores.csv'),check.names=FALSE,colClasses=c(Subject_ID='character'),na.strings=c('','NA','NaN'),stringsAsFactors=FALSE)
  object <- readRDS(file.path(run,'normative_objects.rds'))
  checks <- data.frame(check=character(),pass=logical())
  check <- function(name,value) {checks<<-rbind(checks,data.frame(check=name,pass=isTRUE(value)));if(!isTRUE(value))stop('Audit failed: ',name)}
  keys<-intersect(c('Subject_ID','Timepoint','Platform','Diagnosis','StudyDate'),names(raw))
  check('input rows and identities',nrow(raw)==nrow(d)&&identical(raw[keys],d[keys]))
  check('raw outcome preserved',isTRUE(all.equal(raw[['RVP_Median_latency']],d[['RVP_Median_latency']])))
  hc<-as.character(d$Diagnosis)%in%'HC'
  expected<-d$normative_z_full_reference;expected[hc]<-d$normative_z_cross_validation[hc]
  check('held-out score routing',identical(expected,d$normative_z))
  check('explicit model source',identical(d$score_source,ifelse(hc,'HC_OOF','full_HC_new_subject')))
  md<-prepare_model_data(data,outcome_specs[[1L]],'gamma_raw')
  md$CV_fold<-make_grouped_folds(md,cv_K,seed_main+100L)$folds
  check('raw training data reconstructed',identical(md,object$selected_model_dat))
  check('full fit training identities',identical(as.integer(object$selected_fit$data$row_id_global),md$row_id_global))
  dg<-rvp_add_bfmi(diagnose_brms_fit(object$selected_fit,'RVP_latency','gamma_raw','full'),object$selected_fit)
  check('actual full fit diagnostics',all(dg$validation_pass))
  kf<-object$selected_kfold;check('five saved fold fits',nrow(kf$fits)==5L)
  omitted<-integer()
  for(k in seq_len(nrow(kf$fits))) {
    f<-kf$fits[[k,'fit']];ix<-as.integer(kf$fits[[k,'omitted']]);omitted<-c(omitted,ix)
    check(paste('fold',k,'training rows'),identical(as.integer(f$data$row_id_global),md$row_id_global[-ix]))
    check(paste('fold',k,'held-out subjects'),!any(as.character(md$Subject_ID[ix])%in%as.character(f$data$Subject_ID)))
    dg<-rvp_add_bfmi(diagnose_brms_fit(f,'RVP_latency','gamma_raw',paste0('fold',k)),f)
    check(paste('fold',k,'actual diagnostics'),all(dg$validation_pass))
  }
  check('complete held-out partition',identical(sort(omitted),seq_len(nrow(md))))
  check('saved OOF score matches export',isTRUE(all.equal(as.numeric(object$oof$z),d$normative_z_cross_validation[match(object$oof$row_id_global,d$row_id_global)])))
  eligible<-d$model_eligible
  check('eligible scores accepted',is.logical(eligible)&&!anyNA(eligible)&&all(is.finite(d$normative_z[eligible]))&&all(d$precision_status[eligible]=='ok'))
  check('excluded rows withheld',all(is.na(d$normative_z[!eligible]))&&all(d$precision_status[!eligible]=='ineligible_missing_or_invalid_input'))
  for(label in c('full_population','oof_population')) {
    precision<-read.csv(file.path(run,paste0('RVP_latency_gamma_raw_',label,'_state_cdf_diagnostics.csv')))
    expected_ids<-if(label=='full_population')d$row_id_global[eligible]else md$row_id_global
    check(paste(label,'precision row coverage'),!anyDuplicated(precision$row_id_global)&&identical(sort(precision$row_id_global),sort(expected_ids)))
    target<-if(label=='full_population')d$normative_z_full_reference else d$normative_z_cross_validation
    check(paste(label,'precision scores'),isTRUE(all.equal(precision$z,target[match(precision$row_id_global,d$row_id_global)])))
    check(paste(label,'integration and posterior precision'),all(precision$numerical_ok)&&all(precision$posterior_precision_ok)&&all(is.finite(precision$posterior_mcse_z))&&all(precision$posterior_mcse_z<=.02))
    expected_draws<-if(label=='full_population')rep(nrow(posterior::as_draws_df(object$selected_fit)),nrow(precision))else
      vapply(seq_len(nrow(kf$fits)),function(k)nrow(posterior::as_draws_df(kf$fits[[k,'fit']])),integer(1))[precision$fold]
    check(paste(label,'all posterior draws'),!anyNA(expected_draws)&&all(precision$posterior_draws==expected_draws))
  }
  calibration<-read.csv(file.path(run,'OOF_calibration.csv'))
  check('held-out calibration coverage',nrow(calibration)==1L&&calibration$n==sum(eligible&hc))
  write.csv(checks,file.path(out,'audit_checks.csv'),row.names=FALSE)
  writeLines(paste(nrow(checks),'checks passed'),file.path(out,'audit_status.txt'))
  invisible(checks)
}
