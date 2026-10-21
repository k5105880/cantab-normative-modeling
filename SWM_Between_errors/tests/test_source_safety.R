# Source the configuration and runner without changing the caller's R session.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
thread_names<-c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
  'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS')
original<-list(threads=Sys.getenv(thread_names,unset=NA_character_),search=search(),
  libraries=.libPaths(),contrasts=getOption('contrasts'))
env<-new.env(parent=globalenv())
source(file.path(root,'config.R'),local=env)
source(file.path(root,'run_model.R'),local=env)
stopifnot(identical(original$threads,Sys.getenv(thread_names,unset=NA_character_)),
  identical(original$search,search()),identical(original$libraries,.libPaths()),
  identical(original$contrasts,getOption('contrasts')),
  is.function(env$run_swm),!exists('CONFIG',globalenv(),inherits=FALSE))
source(file.path(root,'R/bootstrap.R'),local=env)
runtime<-env$load_runtime(root)
stopifnot(is.environment(runtime),is.function(runtime$fixed_definition),
  !exists('fixed_definition',globalenv(),inherits=FALSE),
  identical(original$search,search()),identical(original$libraries,.libPaths()),
  identical(original$threads,Sys.getenv(thread_names,unset=NA_character_)))
# Paths inside the source tree must fail before an output directory is created.
for(path in c(root,file.path(root,'outputs'),file.path(dirname(root),'generated'))){
 failure<-tryCatch(runtime$cantab_output_directory(path,dirname(root)),error=identity)
 stopifnot(inherits(failure,'error'),!dir.exists(file.path(root,'outputs')))
}
# A failed direct call also restores temporary execution options and threads.
saved_paths<-Sys.getenv(c('SWM_INPUT','SWM_OUTPUT_RUN'),unset=NA_character_)
Sys.unsetenv(names(saved_paths))
failure<-tryCatch(env$run_swm('check'),error=identity)
Sys.unsetenv(names(saved_paths)[is.na(saved_paths)])
if(any(!is.na(saved_paths)))do.call(Sys.setenv,as.list(saved_paths[!is.na(saved_paths)]))
stopifnot(inherits(failure,'error'),grepl('SWM_INPUT is required',conditionMessage(failure),fixed=TRUE),
  identical(original$threads,Sys.getenv(thread_names,unset=NA_character_)),
  identical(original$contrasts,getOption('contrasts')),
  identical(original$search,search()),identical(original$libraries,.libPaths()))
cat('PASS: source-safe configuration and runner, isolated namespace imports, and external-output guard.\n')
