# Source this file in R/RStudio, then call run_outcome(). Sourcing never runs MCMC.
run_outcome <- local({
  frames <- vapply(sys.frames(),function(x)if(is.null(x$ofile))'' else x$ofile,character(1))
  here <- dirname(normalizePath(tail(frames[nzchar(frames)],1),winslash='/',mustWork=TRUE))
  hash <- function(x) digest::digest(file=x,algo='sha256')
  status <- function(path) if(file.exists(path))readLines(path,n=1L,warn=FALSE) else ''
  function(input_file, output_dir, mode=c('check','fresh','resume'), reference_dir=NULL) {
    mode <- match.arg(mode)
    # One worker; all four sampling chains and iteration counts are preserved.
    cores <- 1L
    packages <- c('brms','cmdstanr','posterior','digest','jsonlite','mgcv','dplyr','loo')
    missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
    if(length(missing))stop('Missing or unloadable R packages: ',paste(missing,collapse=', '))
    input_file <- normalizePath(input_file,winslash='/',mustWork=TRUE)
    if(!is.null(reference_dir))reference_dir <- normalizePath(reference_dir,winslash='/',mustWork=TRUE)
    helper <- file.path(dirname(here),'R','input_data.R')
    prep <- new.env(parent=globalenv());sys.source(helper,envir=prep)
    output_dir <- prep$cantab_output_directory(output_dir,repository=dirname(here))
    if(mode=='resume')output_dir <- normalizePath(output_dir,winslash='/',mustWork=TRUE) else {
      if(file.exists(output_dir))stop('output_dir already exists; choose a new directory or use mode="resume".')
      dir.create(output_dir,recursive=TRUE)
      output_dir <- normalizePath(output_dir,winslash='/',mustWork=TRUE)
    }
    lock <- file.path(output_dir,'.pipeline_lock')
    if(!dir.create(lock,showWarnings=FALSE))stop('Run is locked. Confirm no run is active before removing ',lock)
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    writeLines(as.character(Sys.getpid()),file.path(lock,'pid'))
    completed <- FALSE
    progress_started <- mode!='resume'
    on.exit(if(!completed && progress_started)writeLines('FAILED_OR_INTERRUPTED',file.path(output_dir,'PIPELINE_STATUS.txt')),add=TRUE)
    child_env <- c(CANTAB_JOB_CORES='1',R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),
      OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
      VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',STAN_NUM_THREADS='1',
      MAKEFLAGS='-j1',CMAKE_BUILD_PARALLEL_LEVEL='1')
    env_names <- names(child_env)
    old_env <- Sys.getenv(env_names,unset=NA_character_)
    on.exit({
      Sys.unsetenv(env_names[is.na(old_env)])
      if(any(!is.na(old_env)))do.call(Sys.setenv,as.list(old_env[!is.na(old_env)]))
    },add=TRUE)
    do.call(Sys.setenv,as.list(child_env))
    helper <- file.path(dirname(here),'R','input_data.R')
    results_helper <- file.path(dirname(here),'R','results.R')
    code <- c(file.path(here,c('run.R','build_normative_csv.R','grouped_folds.R','config.R',
      'rti_core.R','rti_cdf.R','run_analysis.R','cantab_design_audit.R')),helper,results_helper)
    code_hash <- setNames(vapply(code,hash,character(1)),basename(code))
    meta_file <- file.path(output_dir,'pipeline_manifest.rds')
    fold_file <- file.path(output_dir,'frozen_subject_folds.csv')
    analysis_input <- file.path(output_dir,'analysis_input.csv')
    if(mode=='resume') {
      meta <- readRDS(meta_file)
      if(!identical(meta$version,'RTI_reaction_R_normative_v1') || meta$mode!='fresh')
        stop('Only a fresh run created by this R launcher can be resumed.')
      if(!identical(meta$code_sha256,code_hash))stop('Resume rejected: R source changed.')
      if(!identical(meta$input_sha256,hash(input_file)) ||
         !identical(meta$input_sha256,hash(file.path(output_dir,'input.csv'))))
        stop('Resume rejected: input changed.')
      if(!identical(meta$analysis_sha256,hash(analysis_input)) || !identical(meta$fold_sha256,hash(fold_file)))
        stop('Resume rejected: prepared input or folds changed.')
      if(!identical(meta$input_audit_sha256,hash(file.path(output_dir,'INPUT_AUDIT_R.rds'))))
        stop('Resume rejected: input audit changed.')
      if(!is.null(reference_dir)) {
        prep <- new.env(parent=globalenv());sys.source(helper,envir=prep)
        prep$cantab_prepare_input(input_file,output_dir,reference_dir,'RTI_reaction')
        supplied <- unique(c(file.path(reference_dir,'input','frozen_subject_folds.csv'),
          file.path(reference_dir,'frozen_subject_folds.csv'),
          file.path(reference_dir,'RTI_reaction','input','frozen_subject_folds.csv')))
        supplied <- supplied[file.exists(supplied)]
        if(length(supplied)>1L || (length(supplied)==1L && hash(supplied)!=meta$fold_sha256))
          stop('Resume rejected: supplied reference folds changed.')
      }
      for(name in names(code_hash))
        if(hash(file.path(output_dir,'source',name))!=code_hash[[name]])stop('Archived source changed: ',name)
      if(as.integer(cores)!=meta$cores)stop('Saved execution settings differ. Use a new output directory.')
    } else {
      before <- hash(input_file)
      if(!file.copy(input_file,file.path(output_dir,'input.csv')) ||
         hash(file.path(output_dir,'input.csv'))!=before || hash(input_file)!=before)stop('Input changed during snapshot.')
      prep <- new.env(parent=globalenv());sys.source(helper,envir=prep)
      analysis_input <- prep$cantab_prepare_input(input_file,output_dir,reference_dir,'RTI_reaction')
      e <- new.env(parent=globalenv());source(file.path(here,'rti_core.R'),local=e)
      data <- e$rti_prepare(analysis_input)
      hc <- data$population[data$population$Diagnosis %in% 'HC',,drop=FALSE]
      supplied <- if(is.null(reference_dir))character() else unique(c(
        file.path(reference_dir,'input','frozen_subject_folds.csv'),
        file.path(reference_dir,'frozen_subject_folds.csv'),
        file.path(reference_dir,'RTI_reaction','input','frozen_subject_folds.csv')))
      supplied <- supplied[file.exists(supplied)]
      if(length(supplied)>1L)stop('Ambiguous reference: more than one frozen_subject_folds.csv.')
      if(length(supplied)) {
        if(!file.copy(supplied,fold_file))stop('Cannot snapshot supplied folds.')
        fold_origin <- list(method='supplied_frozen_participant_folds',source=supplied,
                            allocation_scope='supplied_input')
      } else {
        # Generate a participant-grouped partition using 500 candidates and seed 1334.
        fe <- new.env(parent=asNamespace('dplyr'))
        fe$cv_balance_folds <- TRUE;fe$cv_balance_candidates <- 500L
        fe$cv_balance_include_fsiq <- FALSE
        sys.source(file.path(here,'grouped_folds.R'),envir=fe)
        old_kind <- RNGkind(); has_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
        if(has_seed)old_seed <- get('.Random.seed',envir=.GlobalEnv)
        fold_result <- tryCatch({
          RNGkind('Mersenne-Twister','Inversion','Rejection')
          fe$make_grouped_folds(hc,5L,1334L,TRUE,500L)
        },finally={
          do.call(RNGkind,as.list(old_kind))
          if(has_seed)assign('.Random.seed',old_seed,envir=.GlobalEnv) else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
        })
        map <- unique(data.frame(Subject_ID=as.character(hc$Subject_ID),CV_fold=fold_result$folds))
        write.csv(map,fold_file,row.names=FALSE)
        write.csv(fold_result$audit,file.path(output_dir,'fold_balance_audit.csv'),row.names=FALSE)
        fold_origin <- list(method='new_cohort_balanced_grouped_5fold',seed=1334L,
          candidates=500L,selected_candidate=fold_result$selected_candidate,
          balance_score=fold_result$balance_score,include_fsiq=FALSE,
          RNGkind=c('Mersenne-Twister','Inversion','Rejection'),allocation_scope='supplied_input')
        message('Created a deterministic participant partition for this input.')
      }
      e$rti_prepare(analysis_input,fold_file)
      jsonlite::write_json(fold_origin,file.path(output_dir,'fold_provenance.json'),auto_unbox=TRUE,pretty=TRUE)
      dir.create(file.path(output_dir,'source'))
      if(!all(file.copy(code,file.path(output_dir,'source'))))stop('Cannot snapshot R source.')
      for(name in names(code_hash))
        if(hash(file.path(output_dir,'source',name))!=code_hash[[name]])stop('R source changed during snapshot: ',name)
      meta <- list(version='RTI_reaction_R_normative_v1',mode=mode,cores=as.integer(cores),
        model='normative',input_sha256=before,
        analysis_sha256=hash(analysis_input),fold_sha256=hash(fold_file),
        input_audit_sha256=hash(file.path(output_dir,'INPUT_AUDIT_R.rds')),
        code_sha256=code_hash,fold_provenance=fold_origin,created=as.character(Sys.time()))
      saveRDS(meta,meta_file)
      display_meta <- meta;display_meta$code_sha256 <- as.list(code_hash)
      jsonlite::write_json(display_meta,file.path(output_dir,'pipeline_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
    }
    engine <- file.path(here,'run_analysis.R'); model <- file.path(output_dir,'model')
    rscript <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
    args <- c('--vanilla',engine,'--input',analysis_input,'--folds',fold_file)
    products_file <- file.path(output_dir,'model_products.rds')
    if(file.exists(file.path(model,'manifest.rds'))) {
      args <- c(args,'--resume',model)
      if(startsWith(status(file.path(model,'STATUS.txt')),'NORMATIVE_COMPLETE') && file.exists(products_file)) {
        products <- readRDS(products_file)
        for(name in names(products))
          if(!file.exists(file.path(model,name)) || hash(file.path(model,name))!=products[[name]])
            stop('Saved model output changed: ',name)
        args <- c(args,'--validate-only')
      }
    } else args <- c(args,if(meta$mode=='check')'--check' else '--run','--output-dir',model)
    dir.create(file.path(output_dir,'logs'),showWarnings=FALSE)
    log <- tempfile(paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_model_'),
                    tmpdir=file.path(output_dir,'logs'),fileext='.log')
    progress_started <- TRUE
    writeLines('RUNNING_MODEL_STAGE',file.path(output_dir,'PIPELINE_STATUS.txt'))
    message('RTI reaction model log: ',log)
    code <- system2(rscript,args=vapply(args,shQuote,character(1)),stdout=log,stderr=log)
    if(code!=0L)stop('Model stage failed (exit ',code,'); inspect ',log)
    if(meta$mode=='check') {
      if(status(file.path(model,'STATUS.txt'))!='CHECK_PASSED_NO_MCMC')stop('Model check did not complete.')
      writeLines('CHECK_PASSED_NO_MCMC',file.path(output_dir,'PIPELINE_STATUS.txt'))
      completed <- TRUE
      return(invisible(output_dir))
    }
    if(!startsWith(status(file.path(model,'STATUS.txt')),'NORMATIVE_COMPLETE'))stop('Normative model stage requires review; inspect ',log)
    if(!file.exists(products_file)) {
      paths <- list.files(model,recursive=TRUE,full.names=FALSE)
      paths <- paths[grepl('\\.(csv|rds|R)$',paths)]
      saveRDS(setNames(vapply(file.path(model,paths),hash,character(1)),paths),products_file)
    }
    export <- new.env(parent=globalenv());sys.source(file.path(here,'build_normative_csv.R'),envir=export)
    export$rti_build_normative_csv(model,file.path(output_dir,'analysis'))
    sys.source(helper,envir=export);sys.source(results_helper,envir=export)
    result <- export$cantab_export_results('RTI_reaction',analysis_input,output_dir)
    jsonlite::write_json(list(model='normative',
      input_sha256=meta$input_sha256,analysis_manifest_sha256=hash(file.path(output_dir,'analysis','analysis_manifest.json')),
      completed=as.character(Sys.time())),file.path(output_dir,'final_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
    writeLines(c('NORMATIVE_COMPLETE','model=normative'),file.path(output_dir,'PIPELINE_STATUS.txt'))
    completed <- TRUE
    message('Normative scores ready: ',result)
    invisible(result)
  }
})
