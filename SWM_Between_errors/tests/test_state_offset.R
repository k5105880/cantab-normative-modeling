# The normative model includes exposure offsets in both mu and sigma. No MCMC.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/bootstrap.R'))
runtime<-load_runtime(root)
sys.source(file.path(root,'tests/synthetic_data.R'),envir=runtime)
with(runtime,{
err<-function(x)inherits(tryCatch(force(x),error=identity),'error')
p<-prepare_raw_data(synthetic_raw_data(),CONFIG)
p$hc$CV_fold<-rep(1:5,length.out=nrow(p$hc))
for(id in CONFIG$models) {
 d<-p$hc;def<-fixed_definition(id,d)
 fit<-brms::brm(def$formula,data=d,family=def$family,prior=def$priors,stanvars=def$stanvars,empty=TRUE)
 tt<-brms::brmsterms(fit$formula)
 stopifnot('offset'%in%names(tt$dpars$mu),'offset'%in%names(tt$dpars$sigma))
 swm_state_validate_fit(fit,'hdg')
 bad<-fit;bad$formula<-brms::bf(def$formula$formula,sigma~1+(1|Subject_ID),hu~1,family=def$family)
 failure<-tryCatch(swm_state_validate_fit(bad,'hdg'),error=identity)
 stopifnot(inherits(failure,'error'),grepl('Random auxiliary effects',conditionMessage(failure),fixed=TRUE))
}
cat('PASS: the normative model accepts exposure offsets and still rejects auxiliary random effects. No sampling.\n')

# Real brms prediction with invented draws, avoiding MCMC: both coverage and
# prior diagnostics must use a new participant despite adversarial input IDs.
loadNamespace('rstan')
training<-data.frame(NM_response=c(1L,3L,4L),
 Subject_ID=factor(c('new_subject_heldout','other','prior_new_heldout')))
fixture<-brms::brm(brms::bf(NM_response~1+(1|Subject_ID),sigma~1,hu~1),
 data=training,family=hdg_family(),empty=TRUE)
draw<-list(b_Intercept=rep(log(4),400),b_sigma_Intercept=rep(log(2),400),
 b_hu_Intercept=rep(qlogis(.1),400),sd_Subject_ID__Intercept=rep(.2,400))
for(id in levels(training$Subject_ID))
 draw[[paste0('r_Subject_ID[',id,',Intercept]')]]<-rep(if(id=='other')0 else 5,400)
fixture$fit<-methods::new('stanfit',sim=list(samples=list(draw),warmup2=0L,n_save=400L,chains=1L),mode=0L)
new_person<-data.frame(row_id_global=1:2,NM_response=5L,Subject_ID=factor('heldout'),
 Platform=c('RS6','Connect'))
renamed<-new_person;renamed$Subject_ID<-factor('noncolliding')
collision<-predict_normative_coverage(fixture,new_person,1L,100L,TRUE)
reference<-predict_normative_coverage(fixture,renamed,1L,100L,TRUE)
stopifnot(identical(collision,reference),
 identical(summarize_prior_predictions(fixture,new_person,1L),
           summarize_prior_predictions(fixture,renamed,1L)))
cat('PASS: sampled Between-errors predictions are invariant to training-label collisions.\n')
})
