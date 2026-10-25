# Public runtime regression tests use synthetic data and never sample a model.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/load.R'))
err<-function(x)inherits(tryCatch(force(x),error=identity),'error')
cfg<-swm_config(root);retry<-swm_retry_config(cfg)
stopifnot(cfg$chains==4L,cfg$iter==4000L,cfg$warmup==1500L,cfg$ndraws==8000L,
  cfg$seed==20260921L,cfg$fold_seed==1334L,cfg$fold_candidates==500L,
  retry$iter==5000L,retry$warmup==2500L,retry$adapt_delta==.9999,
  retry$iter-retry$warmup==cfg$iter-cfg$warmup,
  identical(retry[setdiff(names(cfg),c('iter','warmup','adapt_delta'))],cfg[setdiff(names(cfg),c('iter','warmup','adapt_delta'))]),
  err(swm_retry_config(retry)),
  identical(swm_prediction_seed(cfg$seed,'fold',1),20264079L),
  identical(swm_prediction_seed(cfg$seed,'other_participants',1),20278062L))
good<-data.frame(n_divergent=1,treedepth_hits=0,max_rhat=1,min_ess_bulk=500,min_ess_tail=500,
  n_parameters_without_finite_diagnostics=0,min_bfmi=.4,pass=FALSE)
stopifnot(swm_divergence_only(good))
for(n in names(good)) {
  bad<-good
  bad[[n]]<-switch(n,n_divergent=0,treedepth_hits=1,max_rhat=1.01,min_ess_bulk=399,
    min_ess_tail=399,n_parameters_without_finite_diagnostics=1,min_bfmi=.29,pass=TRUE)
  stopifnot(!swm_divergence_only(bad))
}
raw<-data.frame(Subject_ID=c('synthetic_control','synthetic_other','synthetic_missing'),
  Diagnosis=c('HC','Other','HC'),Platform=c('RS6','Connect','RS6'),Timepoint=1,
  Gender=c('F','M','F'),FSIQ=c(100,105,100),Age_baseline_years=c(45,55,60),
  Time_since_baseline_years=0,SWM_Strategy=c(32,2,NA))
x<-swm_adapt_raw(raw);x$CV_fold<-c(1L,NA_integer_,2L)
pred<-list(row_ids=1:2,rows=data.frame(z_mid=c(.4,-1.2),log_mid_p_lower=log(pnorm(c(.4,-1.2))),
  log_mid_p_upper=log(pnorm(c(.4,-1.2),lower.tail=FALSE)),posterior_precision_ok=TRUE,
  numerical_ok=TRUE,posterior_mcse_z=c(.01,.02)))
out<-tempfile('strategy_export_');dir.create(out)
input<-file.path(out,'input.csv');write.csv(raw,input,row.names=FALSE,na='')
ctx<-list(all=x,raw=raw,out=out,input=input)
result<-swm_export(ctx,pred)
stopifnot(identical(result[names(raw)],raw),identical(result$normative_z,c(.4,-1.2,NA_real_)),
  identical(result$normative_model,rep('normative',3)),
  identical(result$normative_score_source,c('participant_OOF','full_HC_new_person','excluded')),
  identical(result$score_status,c('ok','ok','missing_raw')),
  identical(result$precision_status,c('passed','passed','not_scored')),
  file.exists(file.path(out,'SWM_Strategy_scores.csv')),
  err(swm_adapt_raw(result)))
bad<-pred;bad$rows$numerical_ok[1]<-FALSE;stopifnot(err(swm_export(ctx,bad)))
bad<-pred;bad$rows$posterior_precision_ok[1]<-FALSE;stopifnot(err(swm_export(ctx,bad)))
bad<-pred;bad$row_ids<-c(1L,1L);stopifnot(err(swm_export(ctx,bad)))
# Fold-only predictions export correctly when there are no other diagnoses.
controls<-raw;controls$Diagnosis<-'HC';cx<-swm_adapt_raw(controls);cx$CV_fold<-1:3
blocks<-lapply(1:2,function(i)list(row_ids=i,pmf=list(c(.5,.5)),rows=pred$rows[i,,drop=FALSE]))
combined<-swm_combine(blocks,cx[cx$eligible,]);control_out<-tempfile('strategy_controls_');dir.create(control_out)
control_input<-file.path(control_out,'input.csv');write.csv(controls,control_input,row.names=FALSE,na='')
scored<-swm_export(list(all=cx,raw=controls,out=control_out,input=control_input),combined)
stopifnot(identical(scored$normative_z,c(.4,-1.2,NA_real_)),
  all(scored$normative_score_source[1:2]=='participant_OOF'))
cat('PASS: sampling/retry settings, prediction streams, divergence gate, complete raw-row exports and precision refusal.\n')

# Exercise the real brms prediction adapter with invented posterior draws;
# a synthetic new-person label must never resolve to a training participant.
invisible(loadNamespace('rstan'))
training<-data.frame(index=c(1L,3L),grid=12L,
  Subject_ID=factor(c('new_person_heldout','other')))
fixture<-brms::brm(brms::bf(index|vint(grid)~1+(1|Subject_ID),phi~1,zi~1),
  data=training,family=swm_family(),empty=TRUE)
draw<-list(b_Intercept=rep(qlogis(.3),400),b_phi_Intercept=rep(log(30),400),
  b_zi_Intercept=rep(-4,400),sd_Subject_ID__Intercept=rep(.2,400))
draw[['r_Subject_ID[new_person_heldout,Intercept]']]<-rep(5,400)
draw[['r_Subject_ID[other,Intercept]']]<-rep(0,400)
fixture$fit<-methods::new('stanfit',sim=list(samples=list(draw),warmup2=0L,n_save=400L,chains=1L),mode=0L)
new_person<-swm_adapt_raw(data.frame(Subject_ID='heldout',Diagnosis='HC',
  Platform='Connect',Timepoint=1,Gender='F',FSIQ=100,Age_baseline_years=50,
  Time_since_baseline_years=0,SWM_Strategy=6))
new_person$CV_fold<-1L
renamed<-new_person;renamed$Subject_ID<-factor('noncolliding')
collision<-swm_predict_sampled(fixture,new_person,100L,1L)
reference<-swm_predict_sampled(fixture,renamed,100L,1L)
stopifnot(identical(collision$pmf,reference$pmf),
  identical(collision$rows$mean_score,reference$rows$mean_score))
cat('PASS: sampled Strategy predictions are invariant to training-label collisions.\n')

# A real supported ceiling category has finite randomized-PIT moments even
# when its ordinary cumulative probabilities both round to one.
ceiling_person<-renamed;ceiling_person$SWM_Strategy<-ceiling_person$score<-14
ceiling_person$index<-12L
ceiling<-swm_predict_sampled(fixture,ceiling_person,100L,1L)
log_mass<-ceiling$rows$log_mid_p_upper+log(2)
cutoff<-qnorm(log_mass,lower.tail=FALSE,log.p=TRUE)
expected_mean<-exp(dnorm(cutoff,log=TRUE)-log_mass)
stopifnot(is.finite(ceiling$rows$ez),is.finite(ceiling$rows$ez2),
 abs(ceiling$rows$ez-expected_mean)<1e-8,
 abs(ceiling$rows$ez2-(1+cutoff*expected_mean))<1e-8,
 is.finite(swm_calibration(ceiling)$sampled_randomized_pit_mean_z[1]))
cat('PASS: upper-tail Strategy calibration moments remain finite.\n')

# Large but finite logits still define valid beta shapes. Auxiliary sampled
# diagnostics must not reject them after inverse-logit rounding to one.
saturated<-fixture
saturated_draw<-draw;saturated_draw$b_Intercept[]<-40;saturated_draw$lp__<-rep(0,400)
saturated$fit<-methods::new('stanfit',sim=list(samples=rep(list(saturated_draw),4L),
  warmup2=rep(0L,4L),n_save=rep(400L,4L),chains=4L),mode=0L)
saturated_result<-swm_predict(saturated,ceiling_person,100L,1L,
  diagnostics_file=tempfile('strategy_saturated_',fileext='.csv'))
stopifnot(all(is.finite(unlist(saturated_result$pmf))),
 abs(sum(saturated_result$pmf[[1]])-1)<1e-12,
 saturated_result$rows$numerical_ok,saturated_result$rows$posterior_precision_ok,
 abs(saturated_result$rows$z_mid-qnorm(.5+.5*plogis(-4)))<1e-10)
cat('PASS: sampled Strategy probabilities retain finite extreme logits.\n')
