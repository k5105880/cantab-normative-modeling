local({
  paths <- vapply(sys.frames(), function(f) if(is.null(f$ofile)) '' else f$ofile, character(1))
  source(file.path(dirname(normalizePath(tail(paths[nzchar(paths)],1))), 'cantab_design_audit.R'), local=parent.env(environment()))
})
# Fixed normative reaction-time model, validation, and prediction.
rti_model_seed <- function(cfg,fold) {
  # Fixed offset preserves the model's reproducible fit seeds.
  cfg$seed+31000L+fold*10L
}

rti_predict_normative <- function(fit,nd,cfg) {
  p<-rti_predict_marginal(fit,nd,quad_nodes=cfg$quad_nodes,check_nodes=cfg$check_nodes,
    chunk_size=cfg$cdf_chunk_size,z_tolerance=cfg$quad_z_tolerance,
    log_score_tolerance=cfg$quad_log_score_tolerance,max_nodes=cfg$quad_max_nodes,
    posterior_mcse_limit=cfg$posterior_mcse_z_limit)
  rti_apply_marginal_precision(p)
}

rti_num <- function(x,column) {
  missing<-is.na(x)|trimws(as.character(x))=='';v<-suppressWarnings(as.numeric(x))
  if(any(!missing & !is.finite(v))) stop('Invalid numeric value: ',column)
  v[missing]<-NA_real_;v
}

rti_prepare <- function(input,fold_file=NULL) {
  d<-read.csv(input,check.names=FALSE,colClasses='character',na.strings=c('','NA','NaN'),fileEncoding='UTF-8-BOM')
  names(d)<-sub('^\ufeff','',names(d))
  required<-c('Subject_ID','Diagnosis','Timepoint','Gender','Platform','FSIQ','RTI_Median_five-choice_reaction_time')
  if(!all(required %in% names(d)))stop('Missing columns: ',paste(setdiff(required,names(d)),collapse=', '))
  for(v in c('Timepoint','FSIQ','RTI_Median_five-choice_reaction_time'))d[[v]]<-rti_num(d[[v]],v)
  for(stem in c('Age_baseline','Time_since_baseline','Age')) {
    y<-paste0(stem,'_years')
    # Current age is retained for reporting but is not a model predictor.
    if(stem=='Age' && !any(c(stem,y)%in%names(d)))
      d[[y]]<-d$Age_baseline_years+d$Time_since_baseline_years
    if(!any(c(stem,y)%in%names(d)))stop('Missing time/age field: ',y)
    if(stem%in%names(d))d[[stem]]<-rti_num(d[[stem]],stem)
    if(y%in%names(d))d[[y]]<-rti_num(d[[y]],y) else d[[y]]<-d[[stem]]/365.25
    if(stem%in%names(d)) {
      ok<-is.finite(d[[stem]])&is.finite(d[[y]])
      if(any(abs(d[[stem]][ok]/365.25-d[[y]][ok])>1e-6))stop('Conflicting units: ',stem)
    }
  }
  if(anyNA(d[c('Subject_ID','Timepoint','Diagnosis','Platform')]) || any(!nzchar(d$Subject_ID)))stop('Missing participant/visit/group/platform')
  if(anyDuplicated(d[c('Subject_ID','Timepoint')]))stop('Duplicate participant/timepoint')
  if(any(!d$Platform%in%c('RS6','Connect')))stop('Unknown platform')
  if(any(!d$Diagnosis%in%c('HC','BD','MDD')))stop('Unknown diagnosis label')
  if(any(!d$Gender%in%c('F','M') & !is.na(d$Gender)))stop('Unknown Gender coding')
  if(any(d$Timepoint<1|d$Timepoint!=floor(d$Timepoint)))stop('Invalid Timepoint')
  d$NM_raw<-d[['RTI_Median_five-choice_reaction_time']]
  if(any(d$NM_raw<=0,na.rm=TRUE))stop('Reaction time must be positive')
  if(any(d$Time_since_baseline_years<0,na.rm=TRUE))stop('Negative follow-up time')
  if(any(d$Timepoint==1 & (is.na(d$Time_since_baseline_years)|abs(d$Time_since_baseline_years)>1e-6)))stop('Baseline must be time zero')
  if(any(d$Timepoint!=1 & abs(d$Time_since_baseline_years)<1e-6,na.rm=TRUE))stop('Time zero must be baseline')
  d$row_id_global<-seq_len(nrow(d));d$NM_response<-log1p(d$NM_raw)
  d$Age_baseline_10y<-(d$Age_baseline_years-50)/10
  d$Time_since_baseline_1y<-d$Time_since_baseline_years;d$FSIQ_15<-(d$FSIQ-100)/15
  d$Practice_bin<-factor(ifelse(d$Timepoint==1,'0_baseline',ifelse(d$Timepoint==2,'1_second',ifelse(d$Timepoint==3,'2_third','3_fourth_or_later'))),
    levels=c('0_baseline','1_second','2_third','3_fourth_or_later'))
  # Fixed treatment contrasts define the coefficient priors.
  d$Platform<-factor(d$Platform,levels=c('RS6','Connect'));d$Gender<-factor(d$Gender,levels=c('F','M'));d$Subject_ID<-factor(d$Subject_ID)
  for(a in split(d,d$Subject_ID,drop=TRUE)) {
    if(length(unique(a$Diagnosis%in%'HC'))!=1)stop('HC membership changes within a participant')
    if(length(unique(na.omit(a$Age_baseline_years)))>1)stop('Baseline age changes within a participant')
  }
  vars<-c('NM_raw','Age_baseline_10y','Time_since_baseline_1y','FSIQ_15','Gender','Platform','Practice_bin')
  complete<-complete.cases(d[vars]);md<-droplevels(d[complete,])
  # The portable launcher first prepares eligibility before generating a new
  # participant partition. Supplied frozen partitions use the original checks.
  if(is.null(fold_file)) {
    md$CV_fold<-NA_integer_
    return(list(all=d,population=md))
  }
  map<-read.csv(fold_file,stringsAsFactors=FALSE,colClasses=c(Subject_ID='character',CV_fold='integer'))
  if(!all(c('Subject_ID','CV_fold')%in%names(map))||anyNA(map)||anyDuplicated(map$Subject_ID)||!setequal(map$CV_fold,1:5)||!all(map$CV_fold%in%1:5))stop('Invalid frozen fold map')
  md$CV_fold<-map$CV_fold[match(as.character(md$Subject_ID),map$Subject_ID)]
  if(!setequal(as.character(md$Subject_ID[md$Diagnosis%in%'HC']),map$Subject_ID))stop('HC cohort changed; review and freeze a new fold map')
  if(anyNA(md$CV_fold[md$Diagnosis%in%'HC']))stop('HC absent from frozen fold map; create a reviewed new partition for changed cohort')
  list(all=d,population=md)
}

rti_definition <- function() {
  fixed<-c('s(Age_baseline_10y, k=5)',"s(Time_since_baseline_1y, k=4, bs='tp', m=c(2,0))",'Gender','FSIQ_15','Practice_bin','Platform','Platform:Time_since_baseline_1y','Platform:Practice_bin')
  re<-'(0 + Platform + Time_since_baseline_1y | Subject_ID)'
  scale<-'Platform * Practice_bin + Platform * Age_baseline_10y'
  f<-brms::bf(as.formula(paste('NM_response ~',paste(c(fixed,re),collapse=' + '))),as.formula(paste('sigma ~',scale)))
  pr<-c(brms::set_prior('normal(0, 0.20)',class='b'),brms::set_prior('exponential(5)',class='sd'),
    brms::set_prior('exponential(5)',class='sds'),brms::set_prior('lkj(2)',class='cor'),
    brms::set_prior('normal(5.75, 0.35)',class='Intercept'),
    brms::set_prior('normal(-1.71, 0.50)',class='Intercept',dpar='sigma'),
    brms::set_prior('normal(0, 0.25)',class='b',dpar='sigma'),brms::set_prior('gamma(5, 0.5)',class='nu'))
  list(formula=f,priors=pr,family=brms::student(),model='normative')
}

rti_standata <- function(d,def) {
  sd<-brms::make_standata(def$formula,data=d,family=def$family,prior=def$priors)
  cantab_audit_standata(sd,'normative')
  sd
}

rti_diagnose <- function(fit,cfg) {
  dr<-posterior::as_draws_array(fit);dr<-posterior::subset_draws(dr,variable=setdiff(posterior::variables(dr),'lp__'))
  # Cholesky upper triangles and L[1,1] are structural constants, not sampled
  # parameters. Other constant/nonfinite parameters must still fail the gate.
  vn<-posterior::variables(dr);lc<-grep('^L_[0-9]+\\[[0-9]+,[0-9]+\\]$',vn,value=TRUE)
  structural<-vapply(lc,function(v){ij<-as.integer(strsplit(gsub('.*\\[|\\]','',v),',')[[1]]);ij[1]<ij[2]||all(ij==1L)},logical(1))
  if(any(structural))dr<-posterior::subset_draws(dr,variable=setdiff(vn,lc[structural]))
  s<-posterior::summarise_draws(dr,'rhat','ess_bulk','ess_tail')
  np<-brms::nuts_params(fit)
  energy<-split(np$Value[np$Parameter=='energy__'],np$Chain[np$Parameter=='energy__'])
  bfmi<-vapply(energy,function(e)mean(diff(e)^2)/var(e),numeric(1))
  out<-data.frame(n_divergent=sum(np$Parameter=='divergent__'&np$Value==1),
    n_max_treedepth=sum(np$Parameter=='treedepth__'&np$Value>=cfg$max_treedepth),
    max_rhat=max(s$rhat),min_ess_bulk=min(s$ess_bulk),min_ess_tail=min(s$ess_tail),min_bfmi=min(bfmi))
  out$pass<-all(is.finite(as.matrix(out)))&&with(out,n_divergent==0&&n_max_treedepth==0&&max_rhat<cfg$max_rhat&&min_ess_bulk>=cfg$min_ess&&min_ess_tail>=cfg$min_ess&&min_bfmi>=cfg$min_bfmi)
  out
}

rti_fit <- function(train,def,cfg,path,seed,full=FALSE) {
  if(any(!train$Diagnosis%in%'HC'))stop('Non-HC training row')
  rti_standata(train,def)
  if(file.exists(paste0(path,'.rds'))) {
    cp<-readRDS(paste0(path,'.rds'))
    if(!identical(cp$training_rows,train$row_id_global)||!identical(cp$model,def$model))stop('Fit checkpoint metadata mismatch')
    if(!isTRUE(cp$diagnostics$pass))stop('Cached fit did not pass convergence gate')
    return(cp$fit)
  }
  it<-if(full)cfg$full_iter else cfg$cv_iter;wa<-if(full)cfg$full_warmup else cfg$cv_warmup
  ad<-cfg$adapt_delta
  for(attempt in seq_len(2L)) {
    fit<-brms::brm(def$formula,data=train,family=def$family,prior=def$priors,backend='cmdstanr',
      chains=cfg$chains,cores=cfg$cores,iter=it,warmup=wa,seed=seed+attempt-1L,init=0,
      save_pars=brms::save_pars(all=TRUE),refresh=200,
      control=list(adapt_delta=ad,max_treedepth=cfg$max_treedepth))
    dg<-rti_diagnose(fit,cfg);dg$attempt<-attempt
    write.csv(dg,paste0(path,'_diagnostics_attempt',attempt,'.csv'),row.names=FALSE)
    if(isTRUE(dg$pass))break
    ad<-cfg$retry_adapt_delta;it<-it*cfg$retry_iteration_multiplier;wa<-wa*cfg$retry_iteration_multiplier
  }
  if(!isTRUE(dg$pass))stop('Convergence gate failed: ',path)
  cp<-list(fit=fit,training_rows=train$row_id_global,model=def$model,diagnostics=dg)
  saveRDS(cp,paste0(path,'.tmp'));if(!file.rename(paste0(path,'.tmp'),paste0(path,'.rds')))stop('Cannot save fit checkpoint')
  fit
}

rti_calibration <- function(pred,md) {
  d<-md[match(pred$row_id_global,md$row_id_global),];if(anyNA(d$row_id_global))stop('Calibration row mapping failed')
  d$z_estimate<-if('z_estimate'%in%names(pred))pred$z_estimate else pred$z
  d$valid<-pred$numerical_ok
  d$overall<-'all';d$visit<-ifelse(d$Timepoint==1,'baseline','followup')
  d$age<-cut(d$Age_years,c(-Inf,40,60,Inf),labels=c('<=40','40-60','>60'))
  d$platform_visit<-paste(d$Platform,d$visit,sep='/');d$platform_age<-paste(d$Platform,d$age,sep='/')
  groups<-c('overall','Platform','platform_visit','platform_age')
  out<-list()
  for(g in groups)for(a in split(d,d[[g]],drop=TRUE))for(weighting in c('row_equal','subject_equal')) {
    valid<-a$valid&is.finite(a$z_estimate);z<-a$z_estimate[valid];aa<-a[valid,]
    w<-if(weighting=='row_equal')rep(1,length(z))else 1/as.numeric(table(as.character(aa$Subject_ID))[as.character(aa$Subject_ID)])
    m<-if(length(z))weighted.mean(z,w)else NA_real_
    out[[length(out)+1]]<-data.frame(model='normative',score_kind='normative',stratum=g,level=as.character(a[[g]][1]),weighting,
      n_total=nrow(a),n_valid=sum(valid),n_numerical_unresolved=sum(!valid),n_subjects=length(unique(a$Subject_ID)),
      mean_z=m,sd_z=if(length(z))sqrt(weighted.mean((z-m)^2,w))else NA_real_,
      lower_tail=if(length(z))weighted.mean(z<qnorm(.025),w)else NA_real_,
      upper_tail=if(length(z))weighted.mean(z>qnorm(.975),w)else NA_real_,
      central95=if(length(z))weighted.mean(abs(z)<=qnorm(.975),w)else NA_real_)
  }
  do.call(rbind,out)
}
