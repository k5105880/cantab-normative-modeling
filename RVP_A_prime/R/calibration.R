pit_moments <- function(a,b,log_lower=log(a),log_mass=log(b-a),log_upper=log1p(-b)) {
  stopifnot(all(is.finite(a)&is.finite(b)&a>=0&b<=1&a<=b))
  stopifnot(length(a)==length(b),all(lengths(list(log_lower,log_mass,log_upper))==length(a)))
  logs<-c(log_lower,log_mass,log_upper)
  stopifnot(!anyNA(logs),all(logs<=0),
    all(abs(aprime_logadd(aprime_logadd(log_lower,log_mass),log_upper))<1e-10))
  z<-aprime_normal_quantile(aprime_logadd(log_lower,log_mass-log(2)),
    aprime_logadd(log_upper,log_mass-log(2)))
  m<-z;m2<-z^2
  low<-as.numeric(log_lower<log(.025));high<-as.numeric(log_upper<log(.025))
  nodes<-(c(-.906179845938664,-.538469310105683,0,.538469310105683,.906179845938664)+1)/2
  weights<-c(.236926885056189,.478628670499366,.568888888888889,.478628670499366,.236926885056189)/2
  for(i in which(is.finite(log_mass))) {
    ll<-log_lower[i];lm<-log_mass[i];lu<-log_upper[i]
    if(lm<min(ll,lu)+log(1e-3)) {
      # Avoid cancellation between almost equal endpoint densities.
      q<-aprime_normal_quantile(aprime_logadd(ll,lm+log(nodes)),
        aprime_logadd(lu,lm+log1p(-nodes)))
      m[i]<-sum(weights*q);m2[i]<-sum(weights*q*q)
    } else {
      za<-aprime_normal_quantile(ll,aprime_logadd(lm,lu))
      zb<-aprime_normal_quantile(aprime_logadd(ll,lm),lu)
      da<-if(is.finite(za))exp(dnorm(za,log=TRUE)-lm)else 0
      db<-if(is.finite(zb))exp(dnorm(zb,log=TRUE)-lm)else 0
      m[i]<-da-db
      m2[i]<-1+(if(is.finite(za))za*da else 0)-(if(is.finite(zb))zb*db else 0)
    }
    m2[i]<-max(m2[i],m[i]^2)
    # Fraction of the PIT interval below/above the 2.5% tail thresholds.
    fraction<-function(log_start,log_width) {
      cutoff<-log(.025)
      if(log_start>=cutoff)return(0)
      if(aprime_logadd(log_start,log_width)<=cutoff)return(1)
      exp(cutoff+log(-expm1(log_start-cutoff))-log_width)
    }
    low[i]<-fraction(ll,lm);high[i]<-fraction(lu,lm)
  }
  if(any(!is.finite(c(m,m2,low,high))))stop('Nonfinite randomized-PIT moments')
  data.frame(m=m,m2=m2,low=low,high=high)
}
calibrate_scores <- function(run_dir,config) {
  out<-file.path(run_dir,'results');d<-read.csv(file.path(out,'RVP_A_prime_scores.csv'),
    colClasses=c(Subject_ID='character'))
  s<-read.csv(file.path(out,'normative_score_precision.csv'))
  sets<-list(normative=s[grepl('^normative_OOF',s$role),])
  summaries<-moments<-list()
  for(metric in names(sets)) {
    v<-sets[[metric]];a<-d[match(v$row_id_global,d$row_id_global),]
    stopifnot(nrow(v)>0,all(a$Diagnosis=='HC'),!anyDuplicated(a$row_id_global))
    required<-c('log_lower','log_mass','log_upper')
    if(!all(required%in%names(v)))stop('Calibration requires log probabilities from the current scorer')
    a<-cbind(a,pit_moments(v$p_lower,v$p_upper,v$log_lower,v$log_mass,v$log_upper));a$metric<-metric
    moments[[metric]]<-a[c('row_id_global','Subject_ID','Platform','Practice_bin','metric','m','m2','low','high')]
    groups<-list(overall=rep('all',nrow(a)),platform=a$Platform,platform_practice=paste(a$Platform,a$Practice_bin),
      followup_window=ifelse(a$Time_since_baseline_years==0,'baseline',ifelse(a$Time_since_baseline_years<=5,'0_to_5_years','over_5_years')))
    for(group in names(groups))for(level in unique(groups[[group]])) {
      aa<-a[groups[[group]]==level,];by<-split(aa,aa$Subject_ID)
      sums<-do.call(rbind,lapply(by,function(x)c(n=nrow(x),m=sum(x$m),m2=sum(x$m2),low=sum(x$low),high=sum(x$high))))
      calc<-function(x){v<-colSums(x);c(mean=v['m']/v['n'],sd=sqrt(max(0,v['m2']/v['n']-(v['m']/v['n'])^2)),low=v['low']/v['n'],high=v['high']/v['n'])}
      est<-calc(sums);set.seed(config$seed)
      boot<-replicate(config$calibration_bootstraps,calc(sums[sample.int(nrow(sums),nrow(sums),TRUE),,drop=FALSE]))
      ci<-apply(boot,1,quantile,c(.025,.975))
      summaries[[length(summaries)+1L]]<-data.frame(metric=metric,stratum=group,level=level,n_rows=nrow(aa),n_subjects=nrow(sums),
        mean=est[1],sd=est[2],low_tail=est[3],high_tail=est[4],PIT_central95=1-est[3]-est[4],
        mean_CI_low=ci[1,1],mean_CI_high=ci[2,1],sd_CI_low=ci[1,2],sd_CI_high=ci[2,2])
    }
  }
  write.csv(do.call(rbind,summaries),file.path(out,'OOF_calibration.csv'),row.names=FALSE)
  write.csv(do.call(rbind,moments),file.path(out,'OOF_PIT_moments.csv'),row.names=FALSE)
  writeLines(c('Observation-weighted integrated randomized-PIT moments, with subject-cluster bootstrap.',
    'Central95 is PIT central mass, not raw-score predictive interval coverage.',
    'Bootstrap conditions on fitted models; it is not external validation and cannot manufacture unseen tails.',
    'Calibration summarizes held-out healthy controls.'),file.path(out,'calibration_methods.txt'))
}
