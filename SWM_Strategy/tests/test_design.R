#!/usr/bin/env Rscript
# Synthetic public fixture: no participant data or external input is required.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/load.R'));suppressPackageStartupMessages(library(brms))
err<-function(expr) inherits(tryCatch({force(expr);NULL},error=identity),'error')
def<-swm_model_definition()
expected<-read.csv(file.path(root,'tests/expected_priors.csv'),check.names=FALSE,stringsAsFactors=FALSE,na.strings='NA',colClasses='character')
stopifnot(isTRUE(all.equal(as.data.frame(def$prior),expected,check.attributes=FALSE)))
set.seed(20260921)
n<-120L;visits<-5L;id<-rep(seq_len(n),each=visits);tp<-rep(seq_len(visits),n)
platform<-sample(c('RS6','Connect'),length(id),replace=TRUE)
raw<-data.frame(Subject_ID=sprintf('synthetic_%03d',id),Diagnosis='HC',Platform=platform,
  Timepoint=tp,Gender=rep(sample(c('F','M'),n,replace=TRUE),each=visits),
  FSIQ=rep(round(rnorm(n,100,15)),each=visits),
  Age_baseline_years=rep(runif(n,20,85),each=visits),
  Time_since_baseline_years=ifelse(tp==1,0,(tp-1)*runif(length(id),.7,1.3)),
  SWM_Strategy=ifelse(platform=='RS6',sample(8:56,length(id),replace=TRUE),sample(2:14,length(id),replace=TRUE)))
x<-swm_adapt_raw(raw);hc<-droplevels(x[x$eligible&x$Diagnosis%in%'HC',])
hc<-swm_folds(hc,swm_config(root))$data
stanvars<-brms::stanvar(scode=paste(readLines(file.path(root,'stan/swm_discrete_functions.stan')),collapse='\n'),block='functions')

# Interleaved unused non-control levels preserve the allocation of active IDs.
clean<-make_grouped_folds(hc,5L,1334L,n_candidates=5L)
with_unused<-hc;lev<-levels(with_unused$Subject_ID)
with_unused$Subject_ID<-factor(with_unused$Subject_ID,levels=as.vector(rbind(lev,paste0('unused_',seq_along(lev)))))
repaired<-make_grouped_folds(with_unused,5L,1334L,n_candidates=5L)
stopifnot(!anyNA(repaired$folds),identical(clean,repaired),setequal(repaired$folds,1:5))
rows<-list()
for(k in 0:5) {
 train<-if(k==0)hc else hc[hc$CV_fold!=k,]
 sd<-brms::make_standata(def$formula,data=train,family=def$family,prior=def$prior,stanvars=stanvars)
 r<-cantab_audit_standata(sd,paste0('main_fold',k));rows[[k+1L]]<-r
 stopifnot(all(r$full_rank),setequal(r$parameter,c('mu','phi','zi')),
   r$columns[r$parameter=='mu']>ncol(sd$X)+ncol(sd$Xs),any(grepl('Zs_',r$matrices)))
 stopifnot(identical(as.integer(sd$Y),train$index),identical(as.integer(sd$vint1),train$grid))
 # A duplicate hidden solely in a penalized smooth must be rejected.
 sd$Zs_99_1<-sd$X[,1,drop=FALSE]
 stopifnot(err(cantab_audit_standata(sd,'penalty_duplicate')))
 if(k>0)swm_assert_training(list(data=train),train,hc[hc$CV_fold==k,])
}
print(do.call(rbind,rows),row.names=FALSE)
cat('PASS: full-control and all five folds have full-rank distributional designs, unchanged priors and correct score grids.\n')
