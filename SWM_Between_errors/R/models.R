fixed_definition <- function(id, data) {
 if(!identical(id, 'normative'))stop('Only the normative model is available.')
 swm_final_model(data)
}
family_is_hdg<-function(fit)identical(fit$family$name,'hurdle_discrete_gaussian')
family_probabilities<-function(fit,y,mu,scale,hu){
 if(!family_is_hdg(fit))stop('Only the normative hurdle discretized Gaussian family is supported')
 hdg_probabilities(y,mu,scale,hu)
}
check_design <- function(def,data,label=def$id,output_file=NULL) {
 s<-brms::make_standata(def$formula,data=data,family=def$family,prior=def$priors,stanvars=def$stanvars)
 cantab_audit_standata(s,label,output_file)
}
atomic_rds<-function(object,file){
 dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
 tmp<-paste0(file,'.tmp-',Sys.getpid());saveRDS(object,tmp)
 if(!file.rename(tmp,file))stop('Cannot finalize cache ',file)
}
cache_load<-function(path,signature){
 if(!file.exists(path))return(NULL)
 x<-readRDS(path)
 if(!identical(x$signature,signature))stop('Cache signature mismatch: ',path)
 x$value
}
cache_save<-function(value,path,signature)atomic_rds(list(signature=signature,value=value),path)
fit_fixed <- function(id,data,tag,seed,ctx,prior_only=FALSE) {
 def<-fixed_definition(id,data)
 for(attempt in 1:2){
  check_design(def,data,paste0(tag,'_a',attempt),file.path(ctx$run,'aggregate',paste0(tag,'_design_ranks.csv')))
  settings<-list(chains=ctx$config$chains,
   iter=if(prior_only)1200L else if(attempt==1)ctx$config$iter else ctx$config$retry_iter,
   warmup=if(prior_only)600L else if(attempt==1)ctx$config$warmup else ctx$config$retry_warmup,
   adapt_delta=if(attempt==1)ctx$config$adapt_delta else ctx$config$retry_adapt_delta,
   max_treedepth=if(attempt==1)ctx$config$max_treedepth else ctx$config$retry_max_treedepth)
  sig<-digest::digest(list(ctx$signature,id,data,settings,seed,prior_only),algo='sha256')
  file<-file.path(ctx$run,'models',paste0(tag,'_a',attempt,'.rds'))
  obj<-cache_load(file,sig)
  if(is.null(obj)){
   message(format(Sys.time()),' fitting ',tag,' attempt ',attempt)
   fit<-brms::brm(def$formula,data=data,family=def$family,prior=def$priors,stanvars=def$stanvars,
    backend='cmdstanr',chains=settings$chains,cores=min(ctx$config$cores,settings$chains),iter=settings$iter,warmup=settings$warmup,
    seed=seed,init=0,save_pars=brms::save_pars(all=TRUE),sample_prior=if(prior_only)'only'else'no',refresh=500,
    control=list(adapt_delta=settings$adapt_delta,max_treedepth=settings$max_treedepth))
   obj<-list(fit=fit,diagnostics=if(prior_only)NULL else diagnostics(fit,settings$max_treedepth),settings=settings,fit_signature=sig)
   cache_save(obj,file,sig)
  }
  if(prior_only||diagnostic_pass(obj$diagnostics)){obj$file<-file;return(obj)}
  message('Diagnostic gate did not pass; retrying ',tag)
 }
 stop('Diagnostics failed after retry: ',tag,'. Outputs will not be marked complete.')
}
add_prediction_metadata<-function(p,data,id){
 if(!nrow(p))return(p)
 ix<-match(p$row_id_global,data$row_id_global);stopifnot(!anyNA(ix),!anyDuplicated(p$row_id_global))
 vars<-intersect(c('Subject_ID','Platform','Practice_bin','Age_baseline_years','Time_since_baseline_years','NM_response'),names(data))
 p<-cbind(p,data[ix,vars,drop=FALSE]);p$model<-id;p
}
predict_state<-function(fit,data,id,seed,draws,coverage=TRUE,diagnostics_file=tempfile('swm_state_precision_',fileext='.csv')){
 if(!nrow(data))return(data.frame())
 prediction<-predict_normative(fit,data,seed,ndraws=draws,coverage=coverage,diagnostics_file=diagnostics_file)
 add_prediction_metadata(prediction,data,id)
}
