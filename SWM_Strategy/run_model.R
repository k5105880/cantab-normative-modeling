#!/usr/bin/env Rscript
a<-commandArgs(FALSE);root<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE)))
# Package libraries are inherited from the calling R session.
args<-commandArgs(TRUE);mode<-if(length(args))sub('^--mode=','',args[1]) else 'check'
if(!mode%in%c('check','fit','score','all'))stop('Use: run_model.R check|fit|score|all')
source(file.path(root,'R/load.R'))
out<-Sys.getenv('NM_OUTPUT_DIR','')
if(!nzchar(out))stop('NM_OUTPUT_DIR is required; use run_outcome() to supply the output directory.')
input_helper<-new.env(parent=globalenv())
source(file.path(dirname(root),'R','input_data.R'),local=input_helper)
out<-input_helper$cantab_output_directory(out,dirname(root))
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
status<-file.path(out,'RUN_STATUS.txt')
main<-function(){
thread_names<-c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
  'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS')
previous_threads<-Sys.getenv(thread_names,unset=NA_character_)
on.exit({
  Sys.unsetenv(names(previous_threads)[is.na(previous_threads)])
  if(any(!is.na(previous_threads)))do.call(Sys.setenv,as.list(previous_threads[!is.na(previous_threads)]))
},add=TRUE)
do.call(Sys.setenv,as.list(setNames(rep('1',length(thread_names)),thread_names)))
ctx<-swm_context(root,out)
tryCatch({
writeLines(c(paste('RUNNING',mode),as.character(Sys.time())),status)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
swm_write_csv(ctx$hc[,c('row_id_global','Subject_ID','Platform','CV_fold','score','grid','Age_band')],out,'healthy_control_folds.csv')
swm_write_csv(ctx$fold_audit,out,'fold_balance_audit.csv')
swm_write_csv(ctx$all[,c('row_id_global','Subject_ID','Diagnosis','Platform','Timepoint',
  'SWM_Strategy','SWM_Strategy_QC_reason','score','index','grid','eligible','CV_fold')],out,'input_quality_control.csv')

# Complete design and syntax preflight also runs before fitting or cache reuse.
{
  def<-swm_model_definition()
  ranks<-list()
  for(k in 0:5) {
    d<-if(k==0)ctx$hc else ctx$hc[ctx$hc$CV_fold!=k,]
    sd<-brms::make_standata(def$formula,data=d,family=def$family,prior=def$prior,stanvars=ctx$stanvars)
    ranks[[k+1L]]<-cbind(cantab_audit_standata(sd,paste0('main_fold',k)),model='normative',fold=k)
    if(k>0)swm_assert_training(list(data=d),d,ctx$hc[ctx$hc$CV_fold==k,])
  }
  swm_write_csv(do.call(rbind,ranks),out,'design_ranks.csv')
  code<-brms::make_stancode(def$formula,data=ctx$hc,family=def$family,prior=def$prior,stanvars=ctx$stanvars)
  path<-file.path(out,'main_generated.stan');writeLines(code,path)
  cmdstanr::cmdstan_model(path,compile=FALSE)$check_syntax()
  swm_write_csv(as.data.frame(def$prior),out,'main_priors.csv')
  writeLines(capture.output(print(def$formula)),file.path(out,'main_formula.txt'))
  if(mode=='check') {
    writeLines(c('CHECK PASSED: input preparation, full control and all five held-out model data, Stan syntax.',
      'MCMC and final score creation have NOT been run.',
      sprintf('Raw rows=%d; eligible control rows=%d; control participants=%d',nrow(ctx$raw),nrow(ctx$hc),length(unique(ctx$hc$Subject_ID)))),status)
    message('CHECK PASSED; no MCMC fitted.');return(invisible(TRUE))
  }
}

load_or_fit<-function(fold=NULL,fit_allowed=FALSE) {
  swm_checked_fit(ctx,fold,fit_allowed)
}

if(mode%in%c('fit','all')) {
  diagnostics<-list()
  for(k in 1:6) {
    a<-load_or_fit(if(k==6)NULL else k,TRUE);diagnostics[[k]]<-a$diagnostics
    rm(a);invisible(gc())
  }
  swm_write_csv(do.call(rbind,diagnostics),out,'all_fit_diagnostics.csv')
  writeLines('FIT COMPLETE: five held-out control fits and one full control fit. Scores have not yet been exported.',status)
}

if(mode%in%c('score','all')) {
  blocks<-list()
  for(k in 1:5) {
    a<-load_or_fit(k,FALSE);test<-ctx$hc[ctx$hc$CV_fold==k,]
    swm_assert_training(a$fit,ctx$hc[ctx$hc$CV_fold!=k,],test)
    blocks[[k]]<-swm_prediction_cache(ctx,a$fit,test,'fold',k,fit_path=a$fit_path)
    rm(a);invisible(gc())
  }
  oof<-swm_combine(blocks,ctx$hc)
  swm_atomic_rds(list(fingerprint=ctx$fingerprint,pred=oof),file.path(out,'held_out_predictions.rds'))
  swm_write_csv(oof$rows,out,'held_out_scores.csv')
  swm_write_csv(swm_calibration(oof),out,'held_out_calibration.csv')
  a<-load_or_fit(NULL,FALSE)
  others<-ctx$all[ctx$all$eligible & !ctx$all$Diagnosis%in%'HC',]
  subjects<-unique(as.character(others$Subject_ID))
  if(length(subjects)) {
    chunks<-split(subjects,ceiling(seq_along(subjects)/ctx$cfg$prediction_chunk_subjects))
    for(k in seq_along(chunks)) {
      d<-others[as.character(others$Subject_ID)%in%chunks[[k]],]
      blocks[[length(blocks)+1]]<-swm_prediction_cache(ctx,a$fit,d,'other_participants',k,fit_path=a$fit_path)
    }
  }
  pred<-swm_combine(blocks,ctx$all[ctx$all$eligible,])
  result<-swm_export(ctx,pred)
  writeLines(c('NORMATIVE_COMPLETE',
    sprintf('All raw rows=%d; normative z scores=%d',nrow(result),sum(is.finite(result$normative_z))),
    'Sampling diagnostics and any divergence-only recovery attempt are recorded in all_fit_diagnostics.csv.',
    'Calibration and sampled probability-mass diagnostics are retained.',
    'Controls: participant-held-out five-fold predictions. Other diagnoses: full-control new-person predictions.'),status)
  message('COMPLETE: ',file.path(out,'SWM_Strategy_scores.csv'))
}

},error=function(e){writeLines(c('FAILED: not a completed result.',conditionMessage(e),as.character(Sys.time())),status);stop(e)})
}
swm_with_lock(out,main())
