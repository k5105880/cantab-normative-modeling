# Source this file from R or RStudio, then call run_outcome(). Sourcing does not run it.
run_outcome <- local({
  source_file<-tail(Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile)),1L)[[1L]]
  root<-dirname(normalizePath(source_file,winslash='/',mustWork=TRUE))
  function(input_file,output_dir,mode=c('check','fresh','resume'),reference_dir=NULL) {
    mode<-match.arg(mode)
    # One worker; all four sampling chains and iteration counts are preserved.
    cores <- 1L
    for(p in c('jsonlite','digest'))if(!requireNamespace(p,quietly=TRUE))stop('Install required R package: ',p)
    input_file<-normalizePath(input_file,winslash='/',mustWork=TRUE)
    if(!is.null(reference_dir))reference_dir<-normalizePath(reference_dir,winslash='/',mustWork=TRUE)
    if(length(output_dir)!=1L||is.na(output_dir)||!nzchar(trimws(output_dir)))stop('output_dir must be a nonempty path')
    local_env<-new.env(parent=environment())
    for(file in c(file.path(dirname(root),c('R/input_data.R','R/results.R')),file.path(root,'R/prepare_raw.R'),
      file.path(root,'R/portable_pipeline.R'),file.path(root,'R/config.R')))sys.source(file,envir=local_env)
    output_dir<-local_env$cantab_output_directory(output_dir,dirname(root))
    if(mode=='resume') {
      output_dir<-normalizePath(output_dir,winslash='/',mustWork=TRUE)
    } else {
      if(dir.exists(output_dir)&&length(list.files(output_dir,all.files=TRUE,no..=TRUE)))
        stop('Output directory is not empty; choose a new directory or use mode="resume"')
      dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
      output_dir<-normalizePath(output_dir,winslash='/',mustWork=TRUE)
    }
    lock<-file.path(output_dir,'.run.lock')
    if(!dir.create(lock,showWarnings=FALSE))stop('Output directory is locked by another run. If a prior R process crashed, remove its .run.lock directory after confirming it has stopped.')
    on.exit(unlink(lock,recursive=TRUE),add=TRUE)
    # Child R processes inherit the active R/RStudio package library, not a host-specific path.
    env<-c(CANTAB_JOB_CORES=as.character(as.integer(cores)),R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),
      OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
      VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',STAN_NUM_THREADS='1',
      MAKEFLAGS='-j1',CMAKE_BUILD_PARALLEL_LEVEL='1')
    old<-Sys.getenv(names(env),unset=NA_character_)
    on.exit({Sys.unsetenv(names(old)[is.na(old)]);do.call(Sys.setenv,as.list(old[!is.na(old)]))},add=TRUE)
    do.call(Sys.setenv,as.list(env))
    with(local_env,tryCatch({
      if(mode=='resume') {
        manifest<-jsonlite::fromJSON(file.path(output_dir,'identity.json'))
        if(!identical(manifest$profile,'fitting'))stop('Only a fitting run can be resumed')
        aprime_assert_identity(manifest,root,input_file,reference_dir)
        if(!identical(as.integer(manifest$cores),as.integer(cores)))stop('Saved execution settings differ. Use a new output directory.')
        if(!identical(trimws(readLines(file.path(output_dir,'identity.sha256'))),aprime_identity_token(manifest)))stop('Run identity corrupted')
        if(!identical(jsonlite::fromJSON(file.path(output_dir,'config.json')),manifest$configuration))stop('Saved configuration changed')
        aprime_verify_seals(output_dir)
      } else {
        cfg<-jsonlite::fromJSON(file.path(root,'config.json'))
        validate_config(cfg)
        cfg<-apply_runtime_budget(cfg)
        manifest<-list(schema_version=1L,profile=if(mode=='check')'check' else 'fitting',
          raw_path=input_file,raw_sha256=aprime_hash(input_file),
          cores=as.integer(cores),code=aprime_inventory(root),references=aprime_reference_inventory(reference_dir),configuration=cfg)
        aprime_json(manifest,file.path(output_dir,'identity.json'))
        # Use the same parsed representation on fresh and resumed runs.
        manifest<-jsonlite::fromJSON(file.path(output_dir,'identity.json'))
        writeLines(aprime_identity_token(manifest),file.path(output_dir,'identity.sha256'))
        aprime_json(cfg,file.path(output_dir,'config.json'))
      }
      aprime_status(output_dir,if(mode=='check')'RUNNING_CHECK' else 'RUNNING',mode=mode)
      # Optional participant selection and unit validation are shared across outcomes.
      analysis_input<-cantab_prepare_input(input_file,output_dir,reference_dir=reference_dir,outcome='RVP_A_prime')
      # Supplied participant folds may be at the reference root or in reference/.
      supplied_folds<-aprime_reference_folds(reference_dir)
      if(length(supplied_folds)) {
        fold_hash<-aprime_hash(supplied_folds)
        dir.create(file.path(output_dir,'reference'),showWarnings=FALSE)
        snapshot<-file.path(output_dir,'reference','primary_subject_folds.csv')
        if(!file.copy(supplied_folds,snapshot,overwrite=TRUE)||
           !identical(aprime_hash(snapshot),fold_hash)||
           !identical(aprime_hash(supplied_folds),fold_hash))
          stop('Unable to snapshot unchanged supplied participant folds')
      }
      aprime_prepare_raw(analysis_input,file.path(output_dir,'input'))
      aprime_run_stage('check',root,output_dir)
      aprime_run_stage('prepare',root,output_dir)
      aprime_seal(output_dir,'prepare',c(analysis_input,list.files(file.path(output_dir,'input'),full.names=TRUE)))
      aprime_assert_identity(manifest,root,input_file,reference_dir)
      if(mode=='check') {
        aprime_status(output_dir,'CHECK_PASSED_NO_MCMC')
        message('A-prime input, dependencies, six model designs and Stan syntax passed; no MCMC performed.')
        return(invisible(list(status='CHECK_PASSED_NO_MCMC',output_dir=output_dir)))
      }
      aprime_run_stage('fit',root,output_dir)
      aprime_run_stage('score',root,output_dir)
      aprime_seal(output_dir,'score',file.path(output_dir,'results/normative_score_precision.csv'))
      aprime_run_stage('prior',root,output_dir)
      aprime_seal(output_dir,'prior',list.files(file.path(output_dir,'results/prior_direct'),recursive=TRUE,full.names=TRUE))
      aprime_assemble(output_dir)
      aprime_seal(output_dir,'assemble',file.path(output_dir,'results',c('RVP_A_prime_scores.csv','score_availability.csv','score_dictionary.csv')))
      aprime_run_stage('calibrate',root,output_dir)
      aprime_seal(output_dir,'calibrate',file.path(output_dir,'results',c('OOF_calibration.csv','OOF_PIT_moments.csv','calibration_methods.txt')))
      aprime_run_stage('verify_models',root,output_dir)
      checks<-aprime_verify_artifacts(output_dir,manifest,root,input_file,reference_dir)
      write.csv(checks,file.path(output_dir,'results/run_verification.csv'),row.names=FALSE)
      result<-cantab_export_results('RVP_A_prime',analysis_input,output_dir)
      aprime_status(output_dir,'COMPLETE',checks_passed=nrow(checks),score_sha256=aprime_hash(result))
      message('Normative scores ready: ',result)
      invisible(list(status='COMPLETE',output_dir=output_dir,scores=result))
    },error=function(e){aprime_status(output_dir,if(mode=='check')'FAILED_CHECK' else 'FAILED',error=conditionMessage(e));stop(e)}))
  }
})
