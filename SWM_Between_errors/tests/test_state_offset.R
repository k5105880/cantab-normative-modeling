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
})
