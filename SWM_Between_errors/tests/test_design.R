#!/usr/bin/env Rscript
# Regression: an audit that drops Zs columns silently misses duplicates.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/bootstrap.R'))
runtime<-load_runtime(root)
sys.source(file.path(root,'tests/synthetic_data.R'),envir=runtime)
with(runtime,{
err<-function(expr) inherits(tryCatch({force(expr);NULL},error=identity),'error')
x<-prepare_raw_data(synthetic_raw_data(),CONFIG)
x$hc$CV_fold<-make_grouped_folds(x$hc,CONFIG$folds,CONFIG$fold_seed,TRUE,CONFIG$fold_candidates)$folds

# Interleaved unused non-control levels must not leave participants unassigned.
clean<-make_grouped_folds(x$hc,5L,1334L,n_candidates=5L)
with_unused<-x$hc;lev<-levels(with_unused$Subject_ID)
with_unused$Subject_ID<-factor(with_unused$Subject_ID,levels=as.vector(rbind(lev,paste0('unused_',seq_along(lev)))))
repaired<-make_grouped_folds(with_unused,5L,1334L,n_candidates=5L)
stopifnot(!anyNA(repaired$folds),identical(clean,repaired),setequal(repaired$folds,1:5))
ranks<-list()
for(id in CONFIG$models) {
 d<-x$hc
 for(k in 0:5) {
  train<-if(k==0)d else d[d$CV_fold!=k,]
  def<-fixed_definition(id,train)
  sd<-brms::make_standata(def$formula,data=train,family=def$family,prior=def$priors,stanvars=def$stanvars)
  r<-check_design(def,train,paste0(id,'_fold',k));ranks[[length(ranks)+1L]]<-r
  stopifnot(all(r$full_rank),setequal(r$parameter,c('mu','sigma','hu')),
    r$columns[r$parameter=='mu']>ncol(sd$X)+ncol(sd$Xs),any(grepl('Zs_',r$matrices)))
  sd$Zs_sigma_99_1<-sd$X_sigma[,1,drop=FALSE]
  stopifnot(err(cantab_audit_standata(sd,'sigma_penalty_duplicate')))
 }
}
print(do.call(rbind,ranks),row.names=FALSE)
cat('PASS: normative model x full HC/five folds; all distribution parameters audited; inserted penalty duplicates rejected.\n')
})
