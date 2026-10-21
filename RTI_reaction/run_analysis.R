#!/usr/bin/env Rscript
# Internal R model engine. The public entry point is run.R::run_outcome().
rti_verify_source <- function(out,signature) {
  hashes<-c(signature$code,'input.csv'=signature$input,'frozen_subject_folds.csv'=signature$folds)
  for(name in names(hashes))if(!identical(digest::digest(file=file.path(out,'source',name),algo='sha256'),unname(hashes[name])))
    stop('Archived source changed: ',name)
  invisible(TRUE)
}

rti_initialize_output <- function(out,manifest,codefiles,input,fold_file,before_commit=function()NULL) {
  if(file.exists(out))stop('Output already exists: ',out)
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile('.initializing_model_',tmpdir=dirname(out))
  if(!dir.create(tmp))stop('Cannot create initialization staging directory')
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  dir.create(file.path(tmp,'source'))
  if(!all(file.copy(codefiles,file.path(tmp,'source')))||
     !file.copy(input,file.path(tmp,'source/input.csv'))||
     !file.copy(fold_file,file.path(tmp,'source/frozen_subject_folds.csv')))stop('Snapshot copy failed')
  rti_verify_source(tmp,manifest$signature)
  before_commit()
  saveRDS(manifest,file.path(tmp,'manifest.rds'))
  jm<-manifest;jm$signature$code<-as.list(manifest$signature$code);jm$signature$versions<-as.list(manifest$signature$versions)
  jsonlite::write_json(jm,file.path(tmp,'manifest.json'),auto_unbox=TRUE,pretty=TRUE,digits=NA)
  if(file.exists(out)||!file.rename(tmp,out))stop('Cannot publish initialized output')
  invisible(out)
}

main <- function() {
  # Bound threads only while executing; reading config.R has no process side effects.
  thread_names<-c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','STAN_NUM_THREADS')
  previous_threads<-Sys.getenv(thread_names,unset=NA_character_)
  on.exit({
    Sys.unsetenv(names(previous_threads)[is.na(previous_threads)])
    if(any(!is.na(previous_threads)))do.call(Sys.setenv,as.list(previous_threads[!is.na(previous_threads)]))
  },add=TRUE)
  do.call(Sys.setenv,as.list(setNames(rep('1',length(thread_names)),thread_names)))
  here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE),value=TRUE)[1])))
  source(file.path(here,'config.R'));source(file.path(here,'rti_core.R'));source(file.path(here,'rti_cdf.R'))
  args<-commandArgs(TRUE);mode<-'check';resume<-NULL;requested_out<-NULL;validate_only<-FALSE;fold_file<-NULL;input<-NULL
  while(length(args)) {
    a<-args[1];args<-args[-1]
    if(a%in%c('--check','--run'))mode<-substring(a,3) else
    if(a=='--validate-only')validate_only<-TRUE else
    if(a%in%c('--input','--resume','--output-dir','--folds')) {
      if(!length(args))stop('Missing argument for ',a)
      if(a=='--input')input<-args[1]else if(a=='--resume')resume<-args[1]else if(a=='--folds')fold_file<-args[1]else requested_out<-args[1]
      args<-args[-1]
    } else stop('Unknown argument: ',a)
  }
  if(validate_only&&is.null(resume))stop('--validate-only requires --resume')
  if(!is.null(resume)&&!is.null(requested_out))stop('Cannot combine --resume and --output-dir')
  if(is.null(input)||is.null(fold_file))stop('Supply --input and --folds; use run_outcome().')
  if(is.null(resume)&&is.null(requested_out))stop('Supply --output-dir; use run_outcome().')
  packages<-c('brms','cmdstanr','posterior','digest','jsonlite','mgcv')
  absent<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent))stop('Missing R packages: ',paste(absent,collapse=', '))
  cmdstanr::check_cmdstan_toolchain(quiet=TRUE)
  input<-normalizePath(input,mustWork=TRUE);fold_file<-normalizePath(fold_file,mustWork=TRUE)
  if(!is.null(resume)) {
    out<-normalizePath(resume,mustWork=TRUE)
    old<-readRDS(file.path(out,'manifest.rds'));mode<-old$mode
    if(mode!='run')stop('Only analysis runs can resume')
  }
  cfg<-rti_config;d<-rti_prepare(input,fold_file);def<-rti_definition()
  hc<-d$population[d$population$Diagnosis%in%'HC',,drop=FALSE]
  checks<-list();full_checks<-list()
  for(fold in c(0L,cfg$folds)) {
    train<-if(fold==0L)hc else hc[hc$CV_fold!=fold,,drop=FALSE]
    test<-if(fold==0L)hc[FALSE,]else hc[hc$CV_fold==fold,,drop=FALSE]
    if(length(intersect(as.character(train$Subject_ID),as.character(test$Subject_ID))))stop('Fold leakage')
    for(v in c('Platform','Gender','Practice_bin'))
      if(!setequal(unique(as.character(train[[v]])),levels(hc[[v]])))stop('Training fold missing factor level: ',fold,' ',v)
    sd<-rti_standata(train,def);x<-cbind(sd$X,sd$Xs)
    full_checks[[length(full_checks)+1L]]<-cbind(fold,cantab_audit_standata(sd,paste('normative fold',fold)))
    checks[[length(checks)+1L]]<-data.frame(fold,n_train=nrow(train),n_train_subjects=length(unique(train$Subject_ID)),
      n_test=nrow(test),mean_rank=qr(x)$rank,mean_columns=ncol(x),scale_rank=qr(sd$X_sigma)$rank,scale_columns=ncol(sd$X_sigma))
  }
  stancode<-brms::make_stancode(def$formula,data=hc,family=def$family,prior=def$priors)
  versions<-c(R=as.character(getRversion()),vapply(packages,function(p)as.character(utils::packageVersion(p)),character(1)),CmdStan=as.character(cmdstanr::cmdstan_version()))
  codefiles<-file.path(here,c('config.R','rti_core.R','rti_cdf.R','run_analysis.R','cantab_design_audit.R'))
  signature<-list(mode=mode,input=digest::digest(file=input,algo='sha256'),folds=digest::digest(file=fold_file,algo='sha256'),
    code=setNames(vapply(codefiles,function(f)digest::digest(file=f,algo='sha256'),character(1)),basename(codefiles)),config=cfg,versions=versions)
  signature_hash<-digest::digest(signature,algo='sha256')
  if(is.null(resume)) {
    out<-requested_out
    manifest<-list(mode=mode,signature_hash=signature_hash,signature=signature,input_source=input,started=as.character(Sys.time()),output=out)
    rti_initialize_output(out,manifest,codefiles,input,fold_file)
  } else if(!identical(old$signature_hash,signature_hash))stop('Resume rejected: input, code, configuration, mode or software changed. Start a new run.')
  rti_verify_source(out,signature)
  if(validate_only) {cat('RESUME_VALIDATED_NO_MCMC\n');return(invisible(out))}
  lock<-file.path(out,'.running_lock')
  if(!dir.create(lock))stop('Run is locked. Check for an active process before removing ',lock)
  writeLines(as.character(Sys.getpid()),file.path(lock,'pid'))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  status<-file.path(out,'STATUS.txt');writeLines(paste('RUNNING',mode,format(Sys.time())),status)
  complete<-FALSE
  on.exit(if(!complete)writeLines(paste('INTERRUPTED_OR_FAILED',mode,'Inspect pipeline logs and diagnostic CSVs before resuming.'),status),add=TRUE)
  cat('Output directory:',out,'\nMode:',mode,'\n')
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  write.csv(do.call(rbind,checks),file.path(out,'design_checks.csv'),row.names=FALSE)
  write.csv(do.call(rbind,full_checks),file.path(out,'design_audit_full.csv'),row.names=FALSE)
  counts<-data.frame(data_set=c('all_input','eligible'),rows=c(nrow(d$all),nrow(d$population)),
    hc_rows=c(sum(d$all$Diagnosis%in%'HC'),nrow(hc)))
  write.csv(counts,file.path(out,'input_counts.csv'),row.names=FALSE)
  dir.create(file.path(out,'stan'),showWarnings=FALSE)
  sf<-file.path(out,'stan','normative.stan');writeLines(stancode,sf)
  # Syntax/type checking without C++ compilation or sampling.
  model<-cmdstanr::cmdstan_model(sf,compile=FALSE);model$check_syntax(quiet=TRUE)
  writeLines(capture.output(print(def$formula),print(def$priors)),file.path(out,'stan','normative_specification.txt'))
  if(mode=='check') {
    complete<-TRUE;writeLines('CHECK_PASSED_NO_MCMC',status)
    cat('Input, ',length(checks),' design checks, and model syntax passed; no MCMC run.\n',sep='');return(invisible(out))
  }
  dir.create(file.path(out,'fits'),showWarnings=FALSE);dir.create(file.path(out,'predictions'),showWarnings=FALSE)
  md<-d$population;parts<-list()
  for(fold in c(cfg$folds,0L)) {
    isfull<-fold==0L;tag<-paste('normative',if(isfull)'full_HC'else paste0('fold',fold),sep='_')
    seed<-rti_model_seed(cfg,fold)
    train<-if(isfull)hc else hc[hc$CV_fold!=fold,,drop=FALSE]
    nd<-if(isfull)md[!md$Diagnosis%in%'HC',,drop=FALSE]else hc[hc$CV_fold==fold,,drop=FALSE]
    message(format(Sys.time()),' ',tag,' training=',nrow(train),' scoring=',nrow(nd))
    fit<-rti_fit(train,def,cfg,file.path(out,'fits',tag),seed,full=isfull)
    # An all-HC cohort has no rows to score with the full-HC fit.
    if(nrow(nd)) {
      ppfile<-file.path(out,'predictions',paste0(tag,'.rds'))
      if(file.exists(ppfile))p<-readRDS(ppfile)else {
        p<-rti_predict_normative(fit,nd,cfg)
        p$source<-if(isfull)'patient_full_HC_new_subject'else'HC_subject_OOF'
        p$fold<-if(isfull)NA_integer_ else fold
        p$n_posterior_draws<-posterior::ndraws(fit)
        saveRDS(p,paste0(ppfile,'.tmp'));if(!file.rename(paste0(ppfile,'.tmp'),ppfile))stop('Cannot save prediction checkpoint')
      }
      if(!identical(p$row_id_global,nd$row_id_global))stop('Prediction checkpoint row mismatch')
      parts[[length(parts)+1L]]<-p
    }
    if(isfull&&isTRUE(cfg$run_prior_predictive))rti_prior_check(train,def,cfg,file.path(out,'fits',paste0(tag,'_prior')),seed+5L)
    rm(fit);gc(FALSE)
  }
  pred<-do.call(rbind,parts)
  if(anyDuplicated(pred$row_id_global))stop('Duplicate prediction rows')
  write.csv(pred,file.path(out,'normative_prediction_details.csv'),row.names=FALSE)
  scores<-d$all;m<-match(scores$row_id_global,pred$row_id_global)
  for(v in intersect(c('z','z_estimate','numerical_ok','source','score_ok',
    'posterior_mcse_z','posterior_mc_ok','posterior_mc_status','posterior_tail_ess',
    'posterior_n_chains','posterior_iterations_per_chain'),names(pred)))
    scores[[paste('normative',v,sep='_')]]<-pred[[v]][m]
  reason<-ifelse(is.na(m),'missing_response_or_covariate',ifelse(pred$numerical_ok[m],'ok','numerical_unresolved'))
  if('score_status'%in%names(pred))reason[!is.na(m)]<-pred$score_status[m[!is.na(m)]]
  scores$normative_status<-reason;scores$analysis_mode<-mode;scores$model<-'normative'
  # Preserve original response text in CSVs; numeric model data remain unchanged.
  raw_column<-'RTI_Median_five-choice_reaction_time'
  raw_response<-read.csv(file.path(out,'source/input.csv'),check.names=FALSE,
    colClasses='character',na.strings=NULL,fileEncoding='UTF-8-BOM')[[raw_column]]
  if(length(raw_response)!=nrow(scores))stop('Raw response row coverage changed')
  scores[[raw_column]]<-raw_response;scores$NM_raw<-raw_response
  write.csv(scores,file.path(out,'normative_all_rows_scores.csv'),row.names=FALSE,na='NA')
  calibration<-rti_calibration(pred[pred$source=='HC_subject_OOF',],md)
  write.csv(calibration,file.path(out,'HC_OOF_calibration.csv'),row.names=FALSE)
  q<-data.frame(score='normative',n=nrow(pred),numerical_unresolved=sum(!pred$numerical_ok),
    posterior_MC_unresolved=sum(!pred$posterior_mc_ok),max_abs_z_estimate=max(abs(pred$z_estimate)),
    max_quad_z_diff=max(pred$quad_z_diff),max_quad_log_score_diff=max(pred$quad_log_score_diff))
  write.csv(q,file.path(out,'numerical_quality.csv'),row.names=FALSE)
  dgfiles<-list.files(file.path(out,'fits'),pattern='_diagnostics_attempt[0-9]+\\.csv$',full.names=TRUE)
  dg<-lapply(dgfiles,function(f){x<-read.csv(f);x$fit_attempt<-basename(f);x})
  if(length(dg))write.csv(do.call(rbind,dg),file.path(out,'sampling_diagnostics.csv'),row.names=FALSE)
  unresolved<-sum(!pred$score_ok)
  finalstatus<-if(unresolved)'COMPLETE_WITH_UNRESOLVED_NUMERICAL_SCORES'else'NORMATIVE_COMPLETE'
  complete<-TRUE;writeLines(c(finalstatus,paste('Unresolved scores:',unresolved)),status)
  cat(finalstatus,'\nOutput:',out,'\n')
  invisible(out)
}

rti_prior_check <- function(train,def,cfg,path,seed) {
  if(file.exists(paste0(path,'_summary.csv')))return(invisible(NULL))
  rti_standata(train,def)
  fit<-brms::brm(def$formula,data=train,family=def$family,prior=def$priors,sample_prior='only',backend='cmdstanr',
    chains=2,cores=min(2,cfg$cores),iter=1500,warmup=750,seed=seed,init=0,
    control=list(adapt_delta=cfg$adapt_delta,max_treedepth=cfg$max_treedepth),refresh=500)
  # Prior predictive sanity check only; no arithmetic mean on exp(Student).
  nd<-train;nd$Subject_ID<-factor(paste0('prior_',nd$Subject_ID));set.seed(seed)
  pp<-brms::posterior_predict(fit,newdata=nd,allow_new_levels=TRUE,sample_new_levels='gaussian',ndraws=cfg$prior_draws)
  if(any(!is.finite(pp)))stop('Nonfinite prior predictive log draws')
  result<-lapply(split(seq_len(ncol(pp)),paste(nd$Platform,ifelse(nd$Timepoint==1,'baseline','followup'),sep='/')),function(i){
    v<-as.numeric(pp[,i,drop=FALSE]);q<-quantile(v,c(.005,.025,.5,.975,.995))
    data.frame(stratum=paste(nd$Platform[i[1]],ifelse(nd$Timepoint[i[1]]==1,'baseline','followup'),sep='/'),
      log1p_q005=q[1],log1p_q025=q[2],log1p_median=q[3],log1p_q975=q[4],log1p_q995=q[5],
      raw_nonpositive_fraction=mean(v<=0),raw_above_5000ms_fraction=mean(v>log1p(5000)),raw_overflow_fraction=mean(v>log(.Machine$double.xmax)))})
  dg<-rti_diagnose(fit,cfg);write.csv(dg,paste0(path,'_diagnostics.csv'),row.names=FALSE)
  x<-do.call(rbind,result);x$prior_sampler_pass<-dg$pass
  write.csv(x,paste0(path,'_summary.csv'),row.names=FALSE)
  invisible(NULL)
}

if(sys.nframe()==0L)main()
