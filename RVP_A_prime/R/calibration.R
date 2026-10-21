pit_moments <- function(a,b) {
  stopifnot(all(is.finite(a)&is.finite(b)&a>=0&b<=1&a<=b))
  z<-qnorm((a+b)/2);m<-z;m2<-z^2;low<-as.numeric(a<.025);high<-as.numeric(a>.975)
  ix<-which(b>a);za<-qnorm(a[ix]);zb<-qnorm(b[ix]);w<-b[ix]-a[ix]
  ap<-ifelse(is.finite(za),za*dnorm(za),0);bp<-ifelse(is.finite(zb),zb*dnorm(zb),0)
  m[ix]<-(dnorm(za)-dnorm(zb))/w;m2[ix]<-1+(ap-bp)/w
  low[ix]<-pmax(0,pmin(b[ix],.025)-a[ix])/w
  high[ix]<-pmax(0,b[ix]-pmax(a[ix],.975))/w
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
    a<-cbind(a,pit_moments(v$p_lower,v$p_upper));a$metric<-metric
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
