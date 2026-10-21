# Exercise launcher dispatch and resume corruption guards without fitting.
args <- gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE)
here <- dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
local({
  root <- tempfile('reaction launcher ');dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  raw <- expand.grid(Timepoint=1:4,subject=1:12)
  raw$Subject_ID <- sprintf('%03d',raw$subject)
  raw$Diagnosis <- ifelse(raw$subject<=10,'HC','MDD')
  raw$Platform <- ifelse(raw$subject<=5 | raw$subject==11,'RS6','Connect')
  raw$Gender <- ifelse(raw$subject%%2,'F','M')
  raw$FSIQ <- 90+raw$subject
  raw$Age_baseline_years <- 30+raw$subject
  raw$Time_since_baseline_years <- raw$Timepoint-1
  raw$SWM_Between_errors <- 4;raw$SWM_Strategy <- 5
  raw$RVP_A_prime <- .8;raw$RVP_Median_latency <- 400
  raw[['RTI_Median_five-choice_reaction_time']] <- 320
  raw[['RTI_Median_five-choice_movement_time']] <- 250
  raw$subject <- NULL
  input <- file.path(root,'input with spaces.csv');write.csv(raw,input,row.names=FALSE)
  refs <- file.path(root,'reference');dir.create(file.path(refs,'input'),recursive=TRUE)
  write.csv(unique(raw[c('Subject_ID','Diagnosis')]),file.path(refs,'diagnosis_reference.csv'),row.names=FALSE)
  write.csv(data.frame(Subject_ID=sprintf('%03d',1:10),CV_fold=rep(1:5,2)),
            file.path(refs,'input','frozen_subject_folds.csv'),row.names=FALSE)
  calls <- list()
  e <- new.env(parent=globalenv())
  e$system2 <- function(command,args,...) {
    stopifnot(Sys.getenv("CANTAB_JOB_CORES")=="1",
      Sys.getenv("OMP_NUM_THREADS")=="1",Sys.getenv("MAKEFLAGS")=="-j1")
    calls[[length(calls)+1L]] <<- list(command=command,args=args)
    stop('TEST_NO_MCMC',call.=FALSE)
  }
  source(file.path(here,'run.R'),local=e)
  stopifnot(length(calls)==0L)
  global_before <- ls(envir=.GlobalEnv,all.names=TRUE)
  environment_before <- Sys.getenv()
  engine <- new.env(parent=globalenv())
  source(file.path(here,'run_analysis.R'),local=engine)
  source(file.path(here,'rti_core.R'),local=engine)
  stopifnot(exists('main',envir=engine,inherits=FALSE),
    exists('cantab_audit_standata',envir=engine,inherits=FALSE),
    identical(global_before,ls(envir=.GlobalEnv,all.names=TRUE)),
    identical(environment_before,Sys.getenv()),length(calls)==0L)
  out <- file.path(root,'fresh run')
  stopped <- function(expr,pattern) {
    err <- tryCatch(force(expr),error=identity)
    stopifnot(inherits(err,'error'),grepl(pattern,conditionMessage(err),fixed=TRUE))
  }
  blocked_output <- file.path(dirname(here),'reaction_test_output_must_not_exist')
  stopped(e$run_outcome(input,blocked_output,'check'),'outside')
  stopifnot(length(calls)==0L,!file.exists(blocked_output))
  env_before <- Sys.getenv(c('CANTAB_JOB_CORES','R_LIBS','OMP_NUM_THREADS','MAKEFLAGS'),unset=NA_character_)
  stopped(e$run_outcome(input,out,'fresh',reference_dir=refs),'TEST_NO_MCMC')
  stopifnot(length(calls)==1L,any(grepl('--run',calls[[1]]$args,fixed=TRUE)),
    basename(calls[[1]]$command) %in% c('Rscript','Rscript.exe'),
    !dir.exists(file.path(out,'.pipeline_lock')),
    identical(env_before,Sys.getenv(names(env_before),unset=NA_character_)))
  prepared <- read.csv(file.path(out,'analysis_input.csv'),colClasses='character',check.names=FALSE)
  stopifnot(identical(prepared$Subject_ID,raw$Subject_ID),!('Age_years' %in% names(raw)))
  dir.create(file.path(out,'model'))
  saveRDS(list(),file.path(out,'model','manifest.rds'))
  stopped(e$run_outcome(input,out,'resume'),'TEST_NO_MCMC')
  stopifnot(length(calls)==2L,any(grepl('--resume',calls[[2]]$args,fixed=TRUE)))
  original <- readLines(input)
  cat('\n',file=input,append=TRUE)
  stopped(e$run_outcome(input,out,'resume'),'input changed')
  stopifnot(length(calls)==2L)
  writeLines(original,input)
  fold_file <- file.path(out,'frozen_subject_folds.csv')
  cat('\n',file=fold_file,append=TRUE)
  stopped(e$run_outcome(input,out,'resume'),'prepared input or folds changed')
  stopifnot(length(calls)==2L)
  cat('PASS: launcher/engine source isolation, repository output rejection, fresh/resume Rscript dispatch, spaces, textual IDs, environment/lock cleanup and changed-input/fold rejection; fitting intercepted.\n')
})
