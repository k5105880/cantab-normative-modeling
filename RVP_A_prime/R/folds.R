# Outcome-blind grouped five-fold search: fixed seed 1334, 500 candidates.
# FSIQ and the response are deliberately absent from the balance objective.
aprime_generate_folds <- function(hc, seed=1334L, candidates=500L) {
  if(!requireNamespace('loo',quietly=TRUE))stop('Missing R package loo')
  # loo encodes factors as integers: patient-only levels must not leave gaps.
  hc<-droplevels(hc)
  ids<-as.character(hc$Subject_ID)
  groups<-split(seq_len(nrow(hc)),ids)
  if(length(groups)<5L)stop('At least five eligible HC participants are required for five-fold validation')
  profile<-do.call(rbind,lapply(groups,function(ix)data.frame(
    Subject_ID=ids[ix[1]],n_rows=length(ix),n_visits=length(unique(hc$Timepoint[ix])),
    longitudinal=length(unique(hc$Timepoint[ix]))>=2L,
    Age_baseline_10y=median(hc$Age_baseline_10y[ix]),
    Gender=as.character(hc$Gender[ix[1]]),
    platform_pattern=paste(sort(unique(as.character(hc$Platform[ix]))),collapse='+'),
    stringsAsFactors=FALSE)))
  balance<-function(folds) {
    f<-folds[match(profile$Subject_ID,ids)]
    numeric_dispersion<-function(x) {
      ok<-is.finite(x);x<-x[ok];sx<-sd(x)
      if(length(x)<2L||!is.finite(sx)||sx<1e-8)return(0)
      mean(((tapply(x,f[ok],mean)-mean(x))/sx)^2)
    }
    categorical_dispersion<-function(x) {
      x<-as.character(x);lev<-sort(unique(x[!is.na(x)]))
      if(length(lev)<=1L)return(0)
      sum(vapply(lev,function(lv) {
        p0<-mean(x==lv,na.rm=TRUE);pf<-tapply(x==lv,f,mean,na.rm=TRUE)
        mean((pf-p0)^2,na.rm=TRUE)/max(p0*(1-p0),.05)
      },numeric(1)))
    }
    n<-as.numeric(table(factor(f,levels=1:5)))
    rows<-as.numeric(tapply(profile$n_rows,factor(f,levels=1:5),sum))
    2*mean(((n-mean(n))/max(mean(n),1))^2)+
      mean(((rows-mean(rows))/max(mean(rows),1))^2)+
      numeric_dispersion(profile$Age_baseline_10y)+.5*numeric_dispersion(profile$n_visits)+
      categorical_dispersion(profile$Gender)+categorical_dispersion(profile$platform_pattern)+
      categorical_dispersion(profile$longitudinal)
  }
  had_seed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed)old_seed<-get('.Random.seed',envir=.GlobalEnv)
  old_kind<-RNGkind()
  on.exit({
    do.call(RNGkind,as.list(old_kind))
    if(had_seed)assign('.Random.seed',old_seed,envir=.GlobalEnv) else
      if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  },add=TRUE)
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  best<-NULL;best_score<-Inf;selected<-NA_integer_
  for(i in seq_len(candidates)) {
    set.seed(seed+i-1L)
    f<-loo::kfold_split_grouped(K=5L,x=hc$Subject_ID)
    score<-balance(f)
    if(is.finite(score)&&score<best_score){best<-f;best_score<-score;selected<-i}
  }
  if(is.null(best))stop('No valid grouped five-fold partition')
  list(assignment=unique(data.frame(Subject_ID=ids,CV_fold=as.integer(best))),
    audit=list(source='generated_for_input',
      method='Outcome-blind grouped five-fold covariate balance; FSIQ excluded',
      seed=seed,candidates=candidates,selected_candidate=selected,balance_score=best_score))
}
