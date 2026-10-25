#!/usr/bin/env Rscript
# Internal process entry point; sourcing only defines aprime_driver().
aprime_driver <- local({
  source_files <- Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
  script <- if(length(source_files))tail(source_files,1L)[[1L]] else
    sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L])
  root <- dirname(dirname(normalizePath(gsub('~+~',' ',script,fixed=TRUE),mustWork=TRUE)))
  function(args=commandArgs(TRUE)) {
    stopifnot(length(args)==2L);mode<-match.arg(args[1],c('check','prepare','fit','score','prior','calibrate','verify_models'));run_dir<-normalizePath(args[2],mustWork=TRUE)
    source(file.path(root,'R/prepare_data.R'),local=TRUE)
    source(file.path(root,'R/folds.R'),local=TRUE)
    if(!requireNamespace('jsonlite',quietly=TRUE))stop('Missing R package jsonlite')
    config<-jsonlite::fromJSON(file.path(run_dir,'config.json'),simplifyVector=TRUE)
    source(file.path(root,'R/config.R'),local=TRUE);validate_config(config)
    config<-apply_runtime_budget(config)
    jsonlite::write_json(aprime_resource_plan(config),file.path(run_dir,'runtime_resources.json'),pretty=TRUE,auto_unbox=TRUE)
    if(mode=='check') {
      pkgs<-c('brms','posterior','digest','cmdstanr','jsonlite','loo','mgcv')
      missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
      if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
      software<-c(list(R=paste(R.version$major,R.version$minor,sep='.'),cmdstan=as.character(cmdstanr::cmdstan_version())),
        as.list(setNames(vapply(pkgs,function(p)as.character(utils::packageVersion(p)),character(1)),pkgs)))
      if(!nzchar(software$cmdstan)||!dir.exists(cmdstanr::cmdstan_path()))stop('CmdStan is unavailable')
      if(!'xbeta'%in%getNamespaceExports('brms'))stop('brms must include xbeta')
      path<-file.path(run_dir,'software.json')
      if(file.exists(path) && !identical(jsonlite::fromJSON(path),software))stop('R/CmdStan/package versions changed; use a fresh run')
      jsonlite::write_json(software,path,pretty=TRUE,auto_unbox=TRUE)
      cat('R/CmdStan dependencies and fixed precision settings checked.\n');return(invisible(TRUE))
    }
    suppressPackageStartupMessages({library(brms);library(posterior)})
    source(file.path(root,'R/model_spec.R'),local=TRUE);source(file.path(root,'R/fit_models.R'),local=TRUE)
    if(mode=='prepare') {
      reference_dir<-if(dir.exists(file.path(run_dir,'reference')))file.path(run_dir,'reference') else NULL
      prepared<-prepare_data(file.path(run_dir,'input/aligned_raw.csv'),reference_dir,file.path(run_dir,'input'))
      tasks<-fit_tasks(prepared,config,run_dir)
      full_audit<-list()
      audit<-do.call(rbind,lapply(tasks,function(t) {
        model<-build_model(t$train);a<-model_design_audit(model,t$train)
        full_audit[[length(full_audit)+1L]]<<-cbind(kind=t$kind,fold=t$fold,attr(a,'full_audit'))
        if(t$fold==0L) {
          stan_dir<-file.path(run_dir,'stan');dir.create(stan_dir,showWarnings=FALSE)
          sf<-file.path(stan_dir,paste0(t$kind,'.stan'))
          writeLines(brms::make_stancode(model$formula,data=t$train,family=model$family,prior=model$prior),sf)
          cmdstanr::cmdstan_model(sf,compile=FALSE)$check_syntax(quiet=TRUE)
        }
        cbind(kind=t$kind,fold=t$fold,a)
      }))
      write.csv(do.call(rbind,full_audit),file.path(run_dir,'input/design_audit_full.csv'),row.names=FALSE)
      write.csv(audit,file.path(run_dir,'input/model_design_audit.csv'),row.names=FALSE)
      print(read.csv(file.path(run_dir,'input/prepared_counts.csv')),row.names=FALSE)
      cat('All ',length(tasks),' model designs validated before fitting.\n',sep='');return(invisible(TRUE))
    }
    identity<-trimws(readLines(file.path(run_dir,'identity.sha256')))
    prepared<-readRDS(file.path(run_dir,'input/prepared.rds'))
    previous_options<-options(cmdstanr_write_stan_file_dir=file.path(run_dir,'stan'))
    on.exit(options(previous_options),add=TRUE)
    dir.create(file.path(run_dir,'stan'),showWarnings=FALSE)
    source(file.path(root,'R/precision_helpers.R'),local=TRUE)
    source(file.path(root,'R/scoring.R'),local=TRUE);source(file.path(root,'R/calibration.R'),local=TRUE)
    if(mode=='fit') {
      fit_all_models(prepared,config,run_dir,identity)
    } else if(mode=='score') {
      score_all(prepared,config,run_dir,identity)
    } else if(mode=='prior') {
      source(file.path(root,'R/prior_direct.R'),local=TRUE)
      fit<-read_accepted_fit(file.path(run_dir,'models/primary_full.rds'),identity)
      run_prior_direct(fit,file.path(run_dir,'results/prior_direct'),draws=config$prior$draws,seed=config$seed,batch_size=config$prior$batch_size,training_data=prepared$hc)
    } else if(mode=='calibrate') {
      calibrate_scores(run_dir,config)
    } else if(mode=='verify_models') {
      # Rebuild covariates from raw rows; a saved prepared object alone is not evidence.
      temp<-tempfile('rvp_verify_prepared_');dir.create(temp)
      reference_dir<-if(dir.exists(file.path(run_dir,'reference')))file.path(run_dir,'reference') else NULL
      fresh<-prepare_data(file.path(run_dir,'input/aligned_raw.csv'),reference_dir,temp)
      if(!identical(prepared,fresh))stop('Prepared RDS differs from reconstruction from raw data')
      for(name in list.files(temp,pattern='[.]csv$'))
        if(!identical(readLines(file.path(temp,name)),readLines(file.path(run_dir,'input',name))))stop('Prepared CSV differs: ',name)
      result<-validate_all_models(fresh,config,run_dir,identity)
      tasks<-fit_tasks(fresh,config,run_dir)
      hashes<-as.list(setNames(vapply(tasks,function(t)digest::digest(file=t$path,algo='sha256'),character(1)),
        vapply(tasks,function(t)basename(t$path),character(1))))
      jsonlite::write_json(list(identity=identity,model_hashes=hashes),file.path(run_dir,'models/verification_stamp.json'),auto_unbox=TRUE,pretty=TRUE)
      write.csv(result,file.path(run_dir,'models/actual_verification_diagnostics.csv'),row.names=FALSE)
      unlink(temp,recursive=TRUE)
      cat('All ',length(tasks),' saved model objects and exact raw-data training/holdout identities verified.\n',sep='')
    } else stop('Unknown R stage: ',mode)
  }
})
if(sys.nframe()==0L)aprime_driver()
