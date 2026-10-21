# Deterministic, invented assessments for public model-design tests.
synthetic_raw_data <- function() {
 set.seed(20260921)
 n<-120L;visits<-5L;id<-rep(seq_len(n),each=visits);tp<-rep(seq_len(visits),n)
 platform<-sample(c('RS6','Connect'),length(id),replace=TRUE)
 data.frame(Subject_ID=sprintf('synthetic_%03d',id),Diagnosis='HC',Platform=platform,
  Timepoint=tp,Gender=rep(sample(c('F','M'),n,replace=TRUE),each=visits),
  FSIQ=rep(round(rnorm(n,100,15)),each=visits),
  Age_baseline_years=rep(runif(n,20,85),each=visits),
  Time_since_baseline_years=unlist(lapply(seq_len(n),function(i)c(0,cumsum(runif(visits-1,.7,1.3))))),
  SWM_Between_errors=sample(0:50,length(id),replace=TRUE))
}
