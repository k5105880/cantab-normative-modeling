# Synthetic participant records for executable regression tests.
synthetic_normative_input <- function() {
  id<-rep(seq_len(62),each=4);visit<-rep(1:4,62)
  d<-data.frame(Subject_ID=sprintf('synthetic_%03d',id),Diagnosis=ifelse(id<=60,'HC',ifelse(id==61,'BD','MDD')),
    Timepoint=visit,Platform=ifelse((id+visit)%%3==0,'Connect','RS6'),Gender=ifelse(id%%3==0,'F','M'),
    FSIQ=85+(id*11)%%35,Age_baseline_years=23+(id*7)%%48,
    Time_since_baseline_years=ifelse(visit==1,0,(visit-1)*1.3+(id%%11)/25),
    SWM_Between_errors=5+id%%30,SWM_Strategy=ifelse((id+visit)%%3==0,8,32),
    RVP_A_prime=.8+id%%15/100,RVP_Median_latency=350+id%%25*8,check.names=FALSE)
  d[['RTI_Median_five-choice_reaction_time']]<-350+id%%25*4
  d[['RTI_Median_five-choice_movement_time']]<-250+id%%25*5
  d
}
