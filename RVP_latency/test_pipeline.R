# Orchestration regression with explicit model/scoring doubles; no fitting.
here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(here,'core.R'));source(file.path(here,'synthetic_input.R'))
local({
 source(file.path(dirname(here),'R/input_data.R'),local=environment())
 source(file.path(dirname(here),'R/results.R'),local=environment())
 work<-tempfile('normative_pipeline_');dir.create(work);old<-getwd()
 on.exit({setwd(old);unlink(work,recursive=TRUE)},add=TRUE)
 for(controls_only in c(FALSE,TRUE)) {
  out<-file.path(work,if(controls_only)'controls_only'else'mixed');dir.create(out)
  dat<-synthetic_normative_input();if(controls_only)dat<-dat[dat$Diagnosis=='HC',]
  dat[['RVP_Median_latency']][1]<-NA_real_;dat[['RTI_Median_five-choice_movement_time']][1]<-NA_real_
  # A valid source value exceeds write.csv's default numeric output precision.
  dat[['RVP_Median_latency']]<-as.character(dat[['RVP_Median_latency']])
  dat[['RTI_Median_five-choice_movement_time']]<-as.character(dat[['RTI_Median_five-choice_movement_time']])
  dat[['RVP_Median_latency']][2]<-'350.12345678901234'
  dat[['RTI_Median_five-choice_movement_time']][2]<-'350.12345678901234'
  dat$normative_z<-'old_input_metadata'
  input<-file.path(out,'input.csv');write.csv(dat,input,row.names=FALSE,na='')
  e<-normative_context(input,out,prepare_run=TRUE);setwd(out);dir.create(e$checkpoint_dir)
  e$fit_normative_model<-function(model_dat,outcome_prefix,family_id,seed) {
    stopifnot(seed==2235L)
    list(fit=list(data=model_dat),definition=list(model_label='random_slope'))
  }
  e$run_grouped_kfold<-function(fit,folds,save_fits,seed) {
    stopifnot(seed==11235L,save_fits)
    list(fits=lapply(1:5,function(k)list(fit=list(data=fit$data[folds!=k,]),omitted=which(folds==k))))
  }
  e$diagnose_brms_fit<-function(...)data.frame(validation_pass=TRUE)
  prefix<-if(basename(here)=='RTI_movement')'rti'else'rvp'
  assign(paste0(prefix,'_add_bfmi'),function(d,fit)d,e)
  assign(paste0(prefix,'_validate_kfold'),function(...)TRUE,e)
  e$extract_kfold_summary<-function(...)data.frame(elpd_kfold=-100,status='ok')
  e$score_reference_population<-function(fit,pred_dat,family_id)pred_dat$row_id_global/1000
  e$score_cross_validation<-function(kf_obj,family_id,model_dat,rpit_seed) {
    stopifnot(rpit_seed==32234L)
    z<--model_dat$row_id_global/1000;y<-model_dat$NM_response
    list(row_id_global=model_dat$row_id_global,z=z,pit_randomized=pnorm(z),z_rpit=z,y=y,
      yrep=matrix(rep(y,each=16L),16L))
  }
  e$run_prior_predictive_check<-function(model_dat,definition,spec,selected_family,outcome_index) {
    stopifnot(outcome_index==1L);data.frame(status='ok')
  }
  result<-e$run_normative_pipeline()
  h<-result$Diagnosis=='HC';eligible<-result$model_eligible
  stopifnot(identical(result$normative_z[h],result$normative_z_cross_validation[h]),
    identical(result$normative_z[!h],result$normative_z_full_reference[!h]),
    identical(result$score_source,ifelse(h,'HC_OOF','full_HC_new_subject')),
    all(result$precision_status[eligible]=='ok'),all(is.na(result$normative_z[!eligible])),
    file.exists(file.path(out,paste0(basename(here),'_scores.csv'))),
    file.exists(file.path(out,'OOF_calibration.csv')),
    !any(grepl('conditional|sensitivity|Zdelta|noFSIQ',list.files(out))))
  saved<-read.csv(file.path(out,paste0(basename(here),'_scores.csv')),colClasses='character',check.names=FALSE,na.strings=NULL)
  stopifnot(identical(saved[['RVP_Median_latency']][2],'350.12345678901234'),
    identical(saved[['RTI_Median_five-choice_movement_time']][2],'350.12345678901234'),
    !any(saved$normative_z=='old_input_metadata'))
  exported<-cantab_export_results(basename(here),input,out)
  common<-read.csv(exported,colClasses='character',check.names=FALSE,na.strings=NULL)
  stopifnot(identical(common$raw_score[2],'350.12345678901234'))
  e$fit_normative_model<-function(...)stop('Checkpoint must be reused')
  stopifnot(identical(result,e$run_normative_pipeline()))
  e$data$FSIQ_15[2]<-e$data$FSIQ_15[2]+.1
  error<-tryCatch(e$run_normative_pipeline(),error=identity)
  stopifnot(inherits(error,'error'),grepl('training rows or folds changed',conditionMessage(error)))
 }
})
cat('PASS: mixed/all-control score routing, required outputs, exact raw-value export, orchestration seeds, checkpoint reuse and changed-training rejection.\n')

# A failed coefficient diagnostic must not disappear in an na.rm summary.
# Synthetic draws keep posterior diagnostics real without running a sampler.
set.seed(445)
diagnostic_draws <- array(rnorm(1000L*4L*5L),c(1000L,4L,5L),
  dimnames=list(NULL,NULL,c('b_Intercept','b_FSIQ_15','L_1[1,1]','L_1[1,2]','L_1[2,1]')))
diagnostic_draws[,,'L_1[1,1]'] <- 1
diagnostic_draws[,,'L_1[1,2]'] <- 0
sampler <- expand.grid(Iteration=1:1000,Chain=1:4,
  Parameter=c('energy__','divergent__','treedepth__'),stringsAsFactors=FALSE)
sampler$Value <- ifelse(sampler$Parameter=='energy__',rnorm(nrow(sampler)),
  ifelse(sampler$Parameter=='treedepth__',3,0))
as_draws_df.normative_diagnostic_fixture <- function(x,...)posterior::as_draws_df(x$draws)
nuts_params.normative_diagnostic_fixture <- function(object,...)object$sampler
local({
  diagnostic <- new.env(parent=globalenv())
  prefix <- if(basename(here)=='RTI_movement')'rti' else 'rvp'
  source(file.path(here,'model_functions.R'),local=diagnostic)
  source(file.path(here,paste0(prefix,'_validation_helpers.R')),local=diagnostic)
  diagnostic$max_treedepth_fit <- 20
  check_draws <- function(draws) {
    fit <- structure(list(draws=posterior::as_draws_array(draws),sampler=sampler),
      class='normative_diagnostic_fixture')
    result <- diagnostic$diagnose_brms_fit(fit,basename(here),'fixture','random_slope')
    diagnostic[[paste0(prefix,'_add_bfmi')]](result,fit)
  }
  # Cholesky's first diagonal and upper triangle are deterministic by design.
  stopifnot(check_draws(diagnostic_draws)$validation_pass)
  stalled <- diagnostic_draws;stalled[,,'b_FSIQ_15'] <- 0
  stopifnot('A constant sampled coefficient must fail diagnostics' =
    !check_draws(stalled)$validation_pass)
  # A lower-triangular Cholesky element is sampled and must not be exempted.
  stalled <- diagnostic_draws;stalled[,,'L_1[2,1]'] <- 0
  stopifnot('A constant sampled Cholesky element must fail diagnostics' =
    !check_draws(stalled)$validation_pass)
})
cat('PASS: nonfinite coefficient diagnostics fail; only deterministic Cholesky entries are exempt.\n')
