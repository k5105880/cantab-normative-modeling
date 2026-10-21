swm_atomic_rds<-function(x,path) {
  tmp<-paste0(path,'.writing-',Sys.getpid());saveRDS(x,tmp)
  if(!file.rename(tmp,path))stop('Cannot atomically save ',path)
}

swm_fit_diagnostics<-function(fit,label,depth) {
  diag<-swm_diagnostics(fit,label,depth)
  n<-brms::nuts_params(fit);energy<-split(n$Value[n$Parameter=='energy__'],n$Chain[n$Parameter=='energy__'])
  diag$min_bfmi<-min(vapply(energy,function(e)mean(diff(e)^2)/var(e),numeric(1)))
  diag$pass<-with(diag,n_divergent==0 & treedepth_hits==0 & max_rhat<1.01 &
    min_ess_bulk>=400 & min_ess_tail>=400 & n_parameters_without_finite_diagnostics==0 & min_bfmi>=.3)
  diag
}

swm_retry_config<-function(cfg) {
  retry<-cfg
  retry$warmup<-max(cfg$warmup,2500L)
  retry$iter<-retry$warmup+(cfg$iter-cfg$warmup)
  retry$adapt_delta<-max(cfg$adapt_delta,.9999)
  if(retry$warmup==cfg$warmup&&retry$adapt_delta==cfg$adapt_delta)
    stop('Retry settings are already exhausted; investigate model geometry')
  retry
}

swm_divergence_only<-function(d) {
  required<-c('n_divergent','treedepth_hits','max_rhat','min_ess_bulk','min_ess_tail',
    'n_parameters_without_finite_diagnostics','min_bfmi')
  if(!is.data.frame(d)||nrow(d)!=1L||!all(required%in%names(d))||!identical(d$pass,FALSE))return(FALSE)
  if(!all(vapply(d[required],function(x)is.numeric(x)&&length(x)==1L&&is.finite(x),logical(1))))return(FALSE)
  with(d,n_divergent>0 & n_divergent==floor(n_divergent) & treedepth_hits==0 & max_rhat<1.01 &
    min_ess_bulk>=400 & min_ess_tail>=400 & n_parameters_without_finite_diagnostics==0 & min_bfmi>=.3)
}

swm_sample_fit<-function(ctx,def,train,fold,cfg,rawdir) {
  brms::brm(def$formula,data=train,family=def$family,prior=def$prior,stanvars=ctx$stanvars,
    backend='cmdstanr',chains=cfg$chains,iter=cfg$iter,warmup=cfg$warmup,cores=cfg$cores,
    seed=cfg$seed+if(is.null(fold))100L else as.integer(fold),init=cfg$init_radius,refresh=cfg$refresh,
    control=list(adapt_delta=cfg$adapt_delta,max_treedepth=cfg$max_treedepth),
    output_dir=rawdir,save_pars=brms::save_pars(all=TRUE))
}

swm_fit<-function(ctx,fold=NULL,attempt=1L,fit_allowed=TRUE) {
  stopifnot(length(attempt)==1L,attempt%in%c(1L,2L))
  train<-if(is.null(fold))ctx$hc else ctx$hc[ctx$hc$CV_fold!=fold,]
  tag<-paste0('main',if(is.null(fold))'_full_HC' else paste0('_fold',fold))
  cfg<-if(attempt==1L)ctx$cfg else swm_retry_config(ctx$cfg)
  if(attempt==2L)tag<-paste0(tag,'_retry1')
  path<-file.path(ctx$out,paste0(tag,'_fit.rds'));keypath<-file.path(ctx$out,paste0(tag,'_key.rds'))
  def<-swm_model_definition()
  sd<-brms::make_standata(def$formula,data=train,family=def$family,prior=def$prior,stanvars=ctx$stanvars)
  cantab_audit_standata(sd,tag,file.path(ctx$out,paste0(tag,'_design_ranks.csv')))
  key<-list(model='normative',fold=fold,rows=train$row_id_global,fingerprint=ctx$fingerprint,
    formula=capture.output(print(def$formula)),prior=as.data.frame(def$prior))
  if(attempt==2L){key$attempt<-2L;key$sampling<-cfg}
  if(xor(file.exists(path),file.exists(keypath)))stop('Incomplete fit/key pair: ',path)
  if(file.exists(path)) {
    if(!file.exists(keypath)||!identical(readRDS(keypath),key))stop('Stale fit: ',path)
    fit<-readRDS(path)
    swm_assert_training(fit,train,if(is.null(fold))NULL else ctx$hc[ctx$hc$CV_fold==fold,])
    message('Reusing ',tag);return(fit)
  }
  if(!fit_allowed)stop('Missing ',tag,'; run fit or all first')
  rawdir<-file.path(ctx$out,'stan_csv',tag)
  if(attempt==2L) {
    if(dir.exists(rawdir))stop('Incomplete retry evidence for ',tag,'; preserve files and investigate before any further sampling')
    # Persist the retry identity before Stan starts. An interrupted or failed
    # retry must not silently reset the one-retry budget on the next resume.
    swm_atomic_rds(key,keypath)
  }
  dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  message('Fitting ',tag,': ',nrow(train),' HC rows; ',cfg$chains,' chains')
  fit<-swm_sample_fit(ctx,def,train,fold,cfg,rawdir)
  swm_assert_training(fit,train,if(is.null(fold))NULL else ctx$hc[ctx$hc$CV_fold==fold,])
  swm_atomic_rds(fit,path);swm_atomic_rds(key,keypath);fit
}

swm_checked_fit<-function(ctx,fold=NULL,fit_allowed=FALSE) {
  tag<-paste0('main',if(is.null(fold))'_full_HC' else paste0('_fold',fold))
  initial<-file.path(ctx$out,paste0(tag,'_fit.rds'))
  retry_tag<-paste0(tag,'_retry1');retry_path<-file.path(ctx$out,paste0(retry_tag,'_fit.rds'))
  retry_evidence<-any(file.exists(c(retry_path,file.path(ctx$out,paste0(retry_tag,'_key.rds')))))||
    dir.exists(file.path(ctx$out,'stan_csv',retry_tag))
  if(retry_evidence&&!file.exists(initial))stop('Missing initial fit for existing retry: ',tag)
  fit<-swm_fit(ctx,fold,attempt=1L,fit_allowed=fit_allowed)
  diag<-swm_fit_diagnostics(fit,tag,ctx$cfg$max_treedepth)
  swm_write_csv(diag,ctx$out,paste0(tag,'_initial_diagnostics.csv'))
  swm_write_csv(diag,ctx$out,paste0(tag,'_diagnostics.csv'))
  sampling<-ctx$cfg;attempt<-1L;selected<-initial
  if(isTRUE(diag$pass)) {
    if(retry_evidence)stop('Unexpected retry for accepted initial fit: ',tag)
  } else {
    if(!swm_divergence_only(diag))stop('Sampling diagnostics failed for ',tag,'; not eligible for divergence-only retry. Scores withheld.')
    sampling<-swm_retry_config(ctx$cfg)
    rm(fit);invisible(gc())
    fit<-swm_fit(ctx,fold,attempt=2L,fit_allowed=fit_allowed)
    diag<-swm_fit_diagnostics(fit,tag,sampling$max_treedepth)
    swm_write_csv(diag,ctx$out,paste0(retry_tag,'_diagnostics.csv'))
    swm_write_csv(diag,ctx$out,paste0(tag,'_diagnostics.csv'))
    if(!isTRUE(diag$pass))stop('Sampling diagnostics failed after retry for ',tag,'; both fits preserved, scores withheld.')
    attempt<-2L;selected<-retry_path
  }
  diag$attempt<-attempt;diag$iter<-sampling$iter;diag$warmup<-sampling$warmup
  diag$adapt_delta<-sampling$adapt_delta;diag$max_treedepth<-sampling$max_treedepth
  swm_write_csv(diag,ctx$out,paste0(tag,'_diagnostics.csv'))
  list(fit=fit,diagnostics=diag,fit_path=selected)
}

# Prediction RNG streams are indexed by scoring role and batch number.
# These fixed offsets keep the diagnostic PMF draws reproducible independently
# of cache filenames. Normative z scores use deterministic integration.
swm_prediction_seed<-function(seed,role,index) {
  role<-match.arg(role,c('fold','other_participants'))
  if(length(index)!=1L||!is.finite(index)||index<1||index!=floor(index))stop('Invalid prediction batch index')
  digits<-utf8ToInt(as.character(index))
  offset<-if(role=='fold')2766L else 16259L
  start<-if(role=='fold')8L else 18L
  seed+offset+sum(digits*seq.int(start,length.out=length(digits)))
}

swm_prediction_cache<-function(ctx,fit,d,role,index,fit_path) {
  role<-match.arg(role,c('fold','other_participants'))
  seed<-swm_prediction_seed(ctx$cfg$seed,role,index)
  tag<-paste0('main_',if(role=='fold')'fold'else'other_participants_chunk',index)
  path<-file.path(ctx$out,paste0(tag,'_predictions.rds'))
  fit_md5<-tools::md5sum(fit_path)
  if(anyNA(fit_md5))stop('Prediction source fit is missing')
  key<-list(fingerprint=ctx$fingerprint,fit_md5=fit_md5,row_ids=d$row_id_global,
    model='normative',role=role,index=index,seed=seed)
  if(file.exists(path)) {
    x<-readRDS(path);if(!identical(x$key,key))stop('Stale predictions: ',path)
    return(x$pred)
  }
  message('Predicting ',tag,': ',nrow(d),' visits')
  p<-swm_predict(fit,d,ctx$cfg$ndraws,seed,
    diagnostics_file=file.path(ctx$out,paste0(tag,'_state_precision.csv')))
  swm_atomic_rds(list(key=key,pred=p),path);p
}

swm_combine<-function(predictions,d) {
  ids<-unlist(lapply(predictions,`[[`,'row_ids'),use.names=FALSE)
  if(anyDuplicated(ids)||!setequal(ids,d$row_id_global))stop('Missing/duplicated prediction rows')
  ix<-match(d$row_id_global,ids)
  list(model='normative',row_ids=d$row_id_global,
    pmf=unlist(lapply(predictions,`[[`,'pmf'),recursive=FALSE)[ix],
    rows=do.call(rbind,lapply(predictions,`[[`,'rows'))[ix,,drop=FALSE])
}
