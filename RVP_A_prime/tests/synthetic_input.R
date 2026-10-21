# Synthetic records only; no original participants or measurements.
aprime_test_input <- function() {
  id<-rep(seq_len(62),each=4);visit<-rep(1:4,62)
  data.frame(Subject_ID=sprintf('synthetic_%03d',id),Timepoint=visit,
    Diagnosis=ifelse(id<=60,'HC',ifelse(id==61,'BD','MDD')),
    Platform=ifelse((id+visit)%%3==0,'Connect','RS6'),
    Gender=ifelse(id%%3==0,'F','M'),Age_baseline_years=23+(id*7)%%48,
    Time_since_baseline_years=ifelse(visit==1,0,(visit-1)*1.3+(id%%11)/25),
    FSIQ=85+(id*11)%%35,RVP_A_prime=.65+(id%%15)/100+(visit%%3)/100,
    stringsAsFactors=FALSE,check.names=FALSE)
}
