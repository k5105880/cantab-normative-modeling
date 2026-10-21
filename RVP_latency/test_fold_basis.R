# Regression against the actual production model constructors; no MCMC.
args<-commandArgs(TRUE)
here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE)))
scratch<-tempfile('normative_regression_');dir.create(scratch)
source(file.path(here,'core.R'));source(file.path(here,'synthetic_input.R'))
input<-if(length(args))normalizePath(args[1],mustWork=TRUE)else file.path(scratch,'synthetic.csv')
if(!length(args))write.csv(synthetic_normative_input(),input,row.names=FALSE)
e<-normative_context(input,scratch)
setwd(scratch)
spec<-e$outcome_specs[[1]];family<-spec$primary_family
md<-e$prepare_model_data(e$data,spec,family)
folds<-e$make_grouped_folds(md,e$cv_K,e$seed_main+100L)$folds
map<-unique(data.frame(id=as.character(md$Subject_ID),fold=folds))
for(role in 'normative') {
  d<-md
  f<-e$build_model_definition(d,family,outcome_prefix=spec$outcome_prefix)
  x<-e$cantab_brm(f$formula,data=d,family=f$family,prior=f$priors,empty=TRUE)
  fs<-map$fold[match(as.character(d$Subject_ID),map$id)]
  for(k in 1:5) {
    train<-x$data[fs!=k,,drop=FALSE]
    # Empty fits have no sampler metadata. A different unused algorithm label
    # skips that metadata lookup; testmode updates only the same formula/data/bases.
    updated<-update(x,newdata=train,recompile=FALSE,testmode=TRUE,algorithm='meanfield')
    actual<-brms::standata(updated)
    checked<-brms::make_standata(x$formula,data=train,family=x$family,prior=x$prior,
      knots=attr(x$data,'knots',exact=TRUE),drop_unused_levels=attr(x$data,'drop_unused_levels',exact=TRUE),stanvars=x$stanvars)
    keys<-names(actual)[grepl('^(X($|_)|Xs($|_)|Zs_)',names(actual)) & vapply(actual,is.matrix,logical(1))]
    comparison<-all.equal(actual[keys],checked[keys],tolerance=1e-9)
    if(!isTRUE(comparison)) print(comparison)
    stopifnot(length(keys)>0,isTRUE(comparison))
    e$cantab_audit_standata(actual,paste(role,'actual update fold',k))
  }
  cat('ACTUAL_UPDATE_FIVE_FOLD_BASIS_PASS',role,'\n')
}
