# Validate the actual full and held-out model matrices without sampling.
preflight_models <- function(input_file,output_dir) {
  e <- normative_context(input_file,output_dir,prepare_run=TRUE)
  previous <- setwd(e$run_dir);on.exit(setwd(previous),add=TRUE)
  sp<-e$outcome_specs[[1L]];family<-sp$primary_family
  md<-e$prepare_model_data(e$data,sp,family)
  folds<-e$make_grouped_folds(md,e$cv_K,e$seed_main+100L)$folds
  write.csv(unique(data.frame(Subject_ID=as.character(md$Subject_ID),fold=folds)),
    'primary_subject_folds.csv',row.names=FALSE)
  auditfile<-file.path(e$run_dir,'standata_rank.csv')
  if(file.exists(auditfile))unlink(auditfile)
  check<-function(def,dat,label,training_folds=NULL,sample_prior='no') {
    make<-function(d)brms::make_standata(def$formula,data=d,family=def$family,prior=def$priors,sample_prior=sample_prior)
    e$cantab_audit_standata(make(dat),paste(label,'full'),auditfile)
    if(!is.null(training_folds))for(k in sort(unique(training_folds))) {
      train<-dat[training_folds!=k,,drop=FALSE];test<-dat[training_folds==k,,drop=FALSE]
      stopifnot(!any(as.character(train$Subject_ID)%in%as.character(test$Subject_ID)))
      e$cantab_audit_standata(make(train),paste(label,'fold',k),auditfile)
    }
    code<-brms::make_stancode(def$formula,data=dat,family=def$family,prior=def$priors,sample_prior=sample_prior)
    file<-file.path(e$run_dir,paste0(label,'.stan'));writeLines(code,file)
    cmdstanr::cmdstan_model(file,compile=FALSE)$check_syntax(quiet=TRUE)
  }
  definition<-e$build_model_definition(md,family,outcome_prefix=sp$outcome_prefix)
  check(definition,md,'normative',folds)
  fallback<-e$build_model_definition(md,family,outcome_prefix=sp$outcome_prefix,force_random_intercept=TRUE)
  check(fallback,md,'random_intercept')
  check(definition,md,'prior_predictive',sample_prior='only')
  cat('PREFLIGHT_PASS: six main designs, random-intercept fallback and prior predictive design; three Stan programs; no MCMC.\n')
  invisible(e$run_dir)
}
