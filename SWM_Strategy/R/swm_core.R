# Minimum-inflated discretized beta model and normative predictions.
swm_encode_all<-function(x) {
  req<-c('row_id_global','Subject_ID','Diagnosis','Platform','SWM_Strategy',
    'SWM_Strategy_QC_reason','Age_baseline_10y','Age_baseline_years',
    'Time_since_baseline_1y','Timepoint','Practice_bin','Gender','FSIQ_15')
  if(!all(req%in%names(x))||anyNA(x$row_id_global)||anyDuplicated(x$row_id_global)||
     anyNA(x$Subject_ID)||any(!x$Platform%in%c('RS6','Connect')))stop('Invalid input identifiers/platform/columns')
  x$score<-NA_real_;x$index<-NA_integer_;x$grid<-ifelse(x$Platform=='RS6',48L,12L)
  ok<-x$SWM_Strategy_QC_reason=='ok' & is.finite(x$SWM_Strategy)
  e<-swm_encode(x$SWM_Strategy[ok],as.character(x$Platform[ok]))
  x[ok,c('score','index','grid')]<-e
  x$Platform<-factor(x$Platform,levels=c('RS6','Connect'))
  x$Practice_bin<-factor(x$Practice_bin,levels=c('0_baseline','1_second','2_third','3_fourth_or_later'))
  x$Gender<-factor(x$Gender);x$Subject_ID<-factor(x$Subject_ID)
  x$eligible<-ok & complete.cases(x[,req]) & apply(x[,c('Age_baseline_10y','Age_baseline_years',
    'Time_since_baseline_1y','FSIQ_15')],1,function(y)all(is.finite(y)))
  x$Age_band<-cut(x$Age_baseline_years,c(-Inf,40,65,Inf),right=FALSE,
    labels=c('under40','40_to_under65','65plus'))
  x
}

swm_validate_baselines<-function(d) {
  for(ix in split(seq_len(nrow(d)),as.character(d$Subject_ID),drop=TRUE)) {
    base<-ix[which(d$Timepoint[ix]==1 & abs(d$Time_since_baseline_1y[ix])<1e-8)]
    if(length(base)>1)stop('Multiple baseline rows for one participant')
  }
  invisible(TRUE)
}

swm_model_definition<-function() {
  main<-index | vint(grid) ~ s(Age_baseline_10y,k=5)+s(Time_since_baseline_1y,k=4,bs='tp',m=c(2,0))+
    Gender+FSIQ_15+Practice_bin+Platform+Platform:Time_since_baseline_1y+
    Platform:Practice_bin+(1+Time_since_baseline_1y|Subject_ID)
  precision<-phi~Platform+Practice_bin+Age_baseline_10y
  formula<-brms::bf(main,precision,zi~Platform+Practice_bin)
  prior<-c(brms::set_prior('normal(0, 1.5)',class='Intercept'),
    brms::set_prior('normal(0, 0.5)',class='b'),brms::set_prior('exponential(1)',class='sd'),
    brms::set_prior('exponential(1)',class='sds'),brms::set_prior('lkj(2)',class='cor'),
    brms::set_prior('normal(2.7, 1)',class='Intercept',dpar='phi'),
    brms::set_prior('normal(0, 0.4)',class='b',dpar='phi'),
    brms::set_prior('normal(-4, 2)',class='Intercept',dpar='zi'),
    brms::set_prior('normal(0, 1)',class='b',dpar='zi'))
  list(formula=formula,prior=prior,family=swm_family())
}

swm_pit_moments<-function(lo,hi) {
  a<-qnorm(pmax(0,lo));b<-qnorm(pmin(1,hi));pa<-dnorm(a);pb<-dnorm(b)
  width<-hi-lo
  ez<-(pa-pb)/width
  ez2<-1+(ifelse(is.finite(a),a*pa,0)-ifelse(is.finite(b),b*pb,0))/width
  same<-width==0 & lo>0 & lo<1;ez[same]<-a[same];ez2[same]<-a[same]^2
  tail<-width==0 & (lo==0|lo==1);ez[tail]<-ez2[tail]<-NA_real_
  data.frame(ez=ez,ez2=ez2)
}

swm_predict_sampled<-function(fit,d,ndraws=8000,seed=20260921) {
  if(any(!d$eligible))stop('Prediction input contains excluded rows')
  original_ids<-as.character(d$Subject_ID);nd<-d
  nd$Subject_ID<-factor(paste0('new_person_',original_ids))
  set.seed(seed)
  available<-posterior::ndraws(fit)
  draw_ids<-sample(seq_len(available),ndraws,replace=ndraws>available)
  # mgcv's smooth prediction needs a matrix even for a single scoring visit.
  model_nd<-if(nrow(nd)==1L)nd[c(1L,1L),,drop=FALSE]else nd
  prep<-brms::prepare_predictions(fit,newdata=model_nd,allow_new_levels=TRUE,
    sample_new_levels='gaussian',re_formula=NULL,draw_ids=draw_ids)
  pars<-lapply(seq_len(nrow(d)),function(i)swm_parameters(i,prep))
  pmf<-logpmf<-vector('list',nrow(d))
  for(i in seq_len(nrow(d))) {
    p<-pars[[i]]
    logP<-vapply(0:p$grid,function(k)swm_logpmf(k,p$mu,p$phi,p$grid,p$zi),numeric(ndraws))
    P<-exp(logP)
    if(any(!is.finite(P))||any(abs(rowSums(P)-1)>1e-6))stop('Predictive mass invalid')
    P<-P/rowSums(P);pmf[[i]]<-colMeans(P)
    logpmf[[i]]<-swm_average_logpmf(logP)
  }
  set.seed(20260917);uniforms<-runif(nrow(d))
  rows<-swm_rows(pmf,d,'normative',uniforms,logpmf)
  rows$Time_since_baseline_1y<-d$Time_since_baseline_1y
  lo<-mapply(function(p,k)if(k==0)0 else sum(p[seq_len(k)]),pmf,d$index)
  hi<-lo+mapply(function(p,k)p[k+1],pmf,d$index)
  rows<-cbind(rows,pit_lo=lo,pit_hi=hi,swm_pit_moments(lo,hi))
  list(model='normative',row_ids=d$row_id_global,pmf=pmf,
    rows=rows,prediction_seed=seed,ndraws=ndraws)
}

# Visit state uses all posterior draws and deterministic Gaussian integration.
swm_predict<-function(fit,d,ndraws=8000,seed=20260921,
 diagnostics_file=tempfile('swm_state_precision_',fileext='.csv')) {
 if(any(!d$eligible))stop('Prediction input contains excluded rows')
 state<-swm_state_predict_release(fit,d,'beta',diagnostics_file)
 result<-swm_predict_sampled(fit,d,ndraws,seed)
 ix<-match(result$row_ids,state$row_id_global);stopifnot(!anyNA(ix))
 state<-state[ix,,drop=FALSE]
 result$rows$sampled_pmf_draws<-ndraws
 result$rows$pmf_diagnostic_method<-'sampled_new_person_predictive_distribution'
 for(n in c('mid_p','log_mid_p_lower','log_mid_p_upper','numerical_ok',
  'posterior_precision_ok','posterior_mcse_z','posterior_draws','posterior_chains',
  'quad_nodes_used','quad_z_diff','quad_relative_tail_diff','scoring_method'))result$rows[[n]]<-state[[n]]
 result$rows$z_mid<-state$z
 result$rows$numeric_tail_status<-'ok'
 result$state_diagnostics<-state
 result
}
