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

swm_pit_moments<-function(lo,hi,log_lower=log(lo),log_mass=log(hi-lo),log_upper=log1p(-hi)) {
  ez<-ez2<-rep(NA_real_,length(lo))
  nodes<-(c(-.906179845938664,-.538469310105683,0,.538469310105683,.906179845938664)+1)/2
  weights<-c(.236926885056189,.478628670499366,.568888888888889,.478628670499366,.236926885056189)/2
  for(i in which(is.finite(log_mass))) {
    ll<-log_lower[i];lm<-log_mass[i];lu<-log_upper[i]
    if(lm<min(ll,lu)+log(1e-3)) {
      # Integrate a narrow probability interval directly, preserving its small
      # tail instead of subtracting nearly identical endpoint densities.
      q<-if(ll<=lu)qnorm(swm_state_logadd(ll,lm+log(nodes)),log.p=TRUE)else
        qnorm(swm_state_logadd(lu,lm+log1p(-nodes)),lower.tail=FALSE,log.p=TRUE)
      ez[i]<-sum(weights*q);ez2[i]<-sum(weights*q*q)
    } else {
      a<-if(ll<=log(.5))qnorm(ll,log.p=TRUE)else qnorm(swm_state_logadd(lu,lm),lower.tail=FALSE,log.p=TRUE)
      b<-if(lu<=log(.5))qnorm(lu,lower.tail=FALSE,log.p=TRUE)else qnorm(swm_state_logadd(ll,lm),log.p=TRUE)
      pa<-if(is.finite(a))exp(dnorm(a,log=TRUE)-lm)else 0
      pb<-if(is.finite(b))exp(dnorm(b,log=TRUE)-lm)else 0
      ez[i]<-pa-pb
      ez2[i]<-1+(if(is.finite(a))a*pa else 0)-(if(is.finite(b))b*pb else 0)
    }
    ez2[i]<-max(ez2[i],ez[i]^2)
  }
  same<-!is.finite(log_mass)&lo==hi&lo>0&lo<1
  ez[same]<-qnorm(lo[same]);ez2[same]<-ez[same]^2
  data.frame(ez=ez,ez2=ez2)
}

swm_new_person_ids<-function(fit,ids) {
  prefix<-'new_person_'
  known<-union(as.character(fit$data$Subject_ID),levels(fit$data$Subject_ID))
  while(any(paste0(prefix,ids)%in%known))prefix<-paste0(prefix,'_')
  factor(paste0(prefix,ids))
}

swm_predict_sampled<-function(fit,d,ndraws=8000,seed=20260921) {
  if(any(!d$eligible))stop('Prediction input contains excluded rows')
  original_ids<-as.character(d$Subject_ID);nd<-d
  nd$Subject_ID<-swm_new_person_ids(fit,original_ids)
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
    # Keep the linear predictors: inverse-logit rounding can erase a positive
    # beta shape or mixture tail even though the original logits are finite.
    eta<-matrix(brms::get_dpar(prep,'mu',i=i,inv_link=FALSE),ncol=1L)
    zero_eta<-matrix(brms::get_dpar(prep,'zi',i=i,inv_link=FALSE),ncol=1L)
    phi<-matrix(p$phi,ncol=1L)
    logP<-vapply(0:p$grid,function(k)as.numeric(swm_state_log_probabilities(
      k,eta,phi,zero_eta,'beta',p$grid)$mass),numeric(ndraws))
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
  ll<-mapply(function(p,k)swm_logsumexp(p[seq_len(k)]),logpmf,d$index)
  lm<-mapply(function(p,k)p[k+1L],logpmf,d$index)
  lu<-mapply(function(p,k)swm_logsumexp(if(k<length(p)-1L)p[seq.int(k+2L,length(p))]else numeric()),logpmf,d$index)
  rows<-cbind(rows,pit_lo=lo,pit_hi=hi,swm_pit_moments(lo,hi,ll,lm,lu))
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
