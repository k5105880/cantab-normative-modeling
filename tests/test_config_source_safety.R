# Reading either configuration must preserve the caller's environment, even on error.
args <- gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE)
root <- dirname(dirname(normalizePath(sub('^--file=', '', grep('^--file=',args,value=TRUE)[1]))))
local({
  thread_names <- c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS')
  env_names <- c(thread_names,'CANTAB_JOB_CORES')
  original <- Sys.getenv(env_names,unset=NA_character_)
  on.exit({
    Sys.unsetenv(names(original)[is.na(original)])
    if(any(!is.na(original)))do.call(Sys.setenv,as.list(original[!is.na(original)]))
  })
  for(model in c('RTI_reaction','SWM_Strategy')) {
    do.call(Sys.setenv,as.list(setNames(rep('7',length(thread_names)),thread_names)))
    Sys.unsetenv(thread_names[2])
    Sys.setenv(CANTAB_JOB_CORES='2')
    expected <- Sys.getenv(env_names,unset=NA_character_)
    e <- new.env(parent=baseenv())
    sys.source(file.path(root,model,'config.R'),e)
    stopifnot(identical(expected,Sys.getenv(env_names,unset=NA_character_)))
    cfg <- if(model=='RTI_reaction')e$rti_config else e$NORMATIVE_CONFIG
    stopifnot(cfg$cores==1L,cfg$chains==4L)
    Sys.setenv(CANTAB_JOB_CORES='invalid')
    expected <- Sys.getenv(env_names,unset=NA_character_)
    err <- tryCatch(sys.source(file.path(root,model,'config.R'),new.env(parent=baseenv())),
      error=identity)
    stopifnot(inherits(err,'error'),
      grepl('CANTAB_JOB_CORES must be a positive integer',conditionMessage(err),fixed=TRUE),
      identical(expected,Sys.getenv(env_names,unset=NA_character_)))
  }
  cat('PASS: Reaction and Strategy configurations preserve thread settings on source and validation failure.\n')
})
