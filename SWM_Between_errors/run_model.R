#!/usr/bin/env Rscript
# Sourcing defines the runner; model loading and execution happen only on call.
run_swm<-local({
 locate_root<-function(){
  ofiles<-as.character(unlist(lapply(sys.frames(),function(f)f$ofile),use.names=FALSE))
  ofiles<-ofiles[basename(ofiles)=='run_model.R']
  if(length(ofiles))return(dirname(normalizePath(tail(ofiles,1),mustWork=TRUE)))
  f<-grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(f)){
   file<-sub('^--file=','',f[1])
   if(!file.exists(file))file<-gsub('~+~',' ',file,fixed=TRUE)
   return(dirname(normalizePath(file,mustWork=TRUE)))
  }
  stop('Source run_model.R or execute it with Rscript.')
 }
 root<-locate_root()
 function(action='check',workers=1L){
  if(!action%in%c('check','worker','fit','score','collect','all'))stop('Actions: check | all [workers] | fit [workers] | score | collect')
  previous_options<-options(c('contrasts','cmdstanr_write_stan_file_dir'))
  on.exit(options(previous_options),add=TRUE)
  options(contrasts=c('contr.treatment','contr.poly'))
  thread_names<-c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS')
  previous_threads<-Sys.getenv(thread_names,unset=NA_character_)
  on.exit({Sys.unsetenv(names(previous_threads)[is.na(previous_threads)]);if(any(!is.na(previous_threads)))do.call(Sys.setenv,as.list(previous_threads[!is.na(previous_threads)]))},add=TRUE)
  do.call(Sys.setenv,as.list(setNames(rep('1',length(thread_names)),thread_names)))
  runtime<-new.env(parent=environment())
  sys.source(file.path(root,'R/bootstrap.R'),envir=runtime)
  runtime<-runtime$load_runtime(root,runtime)
  ctx<-runtime$resolve_context(root);message('Run folder: ',ctx$run)
  if(action!='worker') { release<-runtime$acquire_run_lock(ctx);on.exit(release(),add=TRUE) }
  if(action=='check')return(invisible(runtime$preflight(ctx)))
  if(action=='worker')return(invisible(runtime$normative_job(as.character(workers),ctx)))
  if(action=='fit')return(invisible(runtime$run_fits(ctx,as.integer(workers))))
  if(action=='score')return(invisible(runtime$run_score(ctx)))
  if(action=='collect')return(invisible(runtime$collect_diagnostics(ctx)))
  runtime$run_fits(ctx,as.integer(workers));runtime$run_score(ctx);runtime$collect_diagnostics(ctx)
  message('Normative execution complete. Review aggregate calibration and model diagnostics.')
  invisible(ctx)
 }
})
if(sys.nframe()==0L){
 a<-commandArgs(TRUE);action<-if(length(a))a[1]else'check';workers<-if(length(a)>1)a[2]else'1'
 tryCatch(run_swm(action,workers),error=function(e){message('STOP: ',conditionMessage(e));quit(status=1)})
}
