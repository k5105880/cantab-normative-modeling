required_raw <- c('Subject_ID','Diagnosis','Platform','Timepoint','Gender','FSIQ','SWM_Between_errors')
prepare_raw_data <- function(raw,config) {
 missing<-setdiff(required_raw,names(raw));if(length(missing))stop('Missing raw columns: ',paste(missing,collapse=', '))
 d<-raw;d$row_id_global<-seq_len(nrow(d))
 to_number<-function(x,nm) {
  out<-suppressWarnings(as.numeric(as.character(x)))
  bad<-!is.na(x)&nzchar(trimws(as.character(x)))&!as.character(x)%in%c('NA','NaN')&is.na(out)
  if(any(bad))stop('Nonnumeric values in ',nm)
  if(any(!is.finite(out)&!is.na(out)))stop('Infinite values in ',nm)
  out
 }
 for(n in c('Timepoint','FSIQ','SWM_Between_errors'))d[[n]]<-to_number(d[[n]],n)
 days<-all(c('Age_baseline','Time_since_baseline')%in%names(d))
 years<-all(c('Age_baseline_years','Time_since_baseline_years')%in%names(d))
 if(!days&&!years)stop('Need Age_baseline and Time_since_baseline in DAYS, or explicitly named *_years columns.')
 if(days){
  d$Age_baseline_years<-to_number(d$Age_baseline,'Age_baseline')/365.25
  d$Time_since_baseline_years<-to_number(d$Time_since_baseline,'Time_since_baseline')/365.25
  if(years && (any(abs(d$Age_baseline_years-to_number(raw$Age_baseline_years,'Age_baseline_years'))>1e-5,na.rm=TRUE)||any(abs(d$Time_since_baseline_years-to_number(raw$Time_since_baseline_years,'Time_since_baseline_years'))>1e-5,na.rm=TRUE)))stop('Conflicting day and year columns')
 } else {
  d$Age_baseline_years<-to_number(d$Age_baseline_years,'Age_baseline_years')
  d$Time_since_baseline_years<-to_number(d$Time_since_baseline_years,'Time_since_baseline_years')
 }
 if(any(d$Age_baseline_years<5|d$Age_baseline_years>110,na.rm=TRUE))stop('Age outside supported 5-110 years; check day/year units.')
 if(any(d$Time_since_baseline_years<0|d$Time_since_baseline_years>50,na.rm=TRUE))stop('Invalid elapsed time; check units.')
 if(any(d$SWM_Between_errors<0|d$SWM_Between_errors!=round(d$SWM_Between_errors),na.rm=TRUE))stop('SWM errors must be nonnegative integers.')
 if(any(d$Timepoint<1|d$Timepoint!=round(d$Timepoint),na.rm=TRUE))stop('Timepoint must be positive integer.')
 for(n in c('Subject_ID','Diagnosis','Platform','Gender'))d[[n]]<-trimws(as.character(d[[n]]))
 if(anyNA(d$Subject_ID)||any(!nzchar(d$Subject_ID)))stop('Missing Subject_ID')
 if(anyNA(d$Diagnosis)||any(!d$Diagnosis%in%config$groups))stop('Diagnosis must be one of configured groups; no relabeling is automatic.')
 if(anyNA(d$Platform)||any(!d$Platform%in%names(config$trials)))stop('Unknown platform')
 if(anyDuplicated(paste(d$Subject_ID,d$Timepoint,sep='::')))stop('Duplicate subject/timepoint')
 for(ix in split(seq_len(nrow(d)),d$Subject_ID)){
  if(length(unique(d$Diagnosis[ix]))>1)stop('Diagnosis changes within a subject; resolve grouping before analysis.')
  ag<-d$Age_baseline_years[ix];ag<-ag[is.finite(ag)]
  if(length(ag)>1&&diff(range(ag))>1e-5)stop('Inconsistent baseline age within subject')
  tt<-d$Time_since_baseline_years[ix];tp<-d$Timepoint[ix]
  if(any(abs(tt[tp==1])>1e-6,na.rm=TRUE)||any(tt[tp>1]<=0,na.rm=TRUE))stop('Baseline/timepoint mismatch')
  ok<-is.finite(tt)&is.finite(tp)
  if(sum(ok)>1&&any(diff(tt[ok][order(tp[ok])])<=0))stop('Visit times must increase with Timepoint')
 }
 d$Age_baseline_10y<-(d$Age_baseline_years-50)/10
 d$Time_since_baseline_1y<-d$Time_since_baseline_years
 d$FSIQ_15<-(d$FSIQ-100)/15
 labels<-c('0_baseline','1_second','2_third','3_fourth_or_later')
 d$Practice_bin<-factor(labels[pmin(d$Timepoint,4)],levels=labels)
 d$Platform<-factor(d$Platform,levels=c('RS6','Connect'))
 d$Gender<-factor(d$Gender);d$Subject_ID<-factor(d$Subject_ID)
 d$Diagnosis<-factor(d$Diagnosis,levels=config$groups)
 d$trial_exposure<-unname(config$trials[as.character(d$Platform)])
 d$log_trial_exposure<-log(d$trial_exposure)
 d$Connect_followup<-as.numeric(d$Platform=='Connect'&d$Practice_bin!='0_baseline')
 d$NM_response<-d$NM_raw<-d$SWM_Between_errors
 vars<-c('NM_response','Age_baseline_10y','Time_since_baseline_1y','Gender','FSIQ_15','Platform','Practice_bin','Subject_ID','Timepoint','Diagnosis')
 ok<-complete.cases(d[vars]);d$score_eligible<-ok
 d$score_exclusion<-vapply(seq_len(nrow(d)),function(i)paste(vars[is.na(d[i,vars])],collapse=';'),character(1))
 valid<-droplevels(d[ok,,drop=FALSE]);hc<-droplevels(valid[valid$Diagnosis==config$control_label,,drop=FALSE])
 if(!nrow(hc))stop('No complete-case HC rows')
 list(raw=raw,all=d,valid=valid,hc=hc,units=if(days)'days converted to years'else'explicit years')
}
align_scores<-function(target,pred){
 if(anyDuplicated(target$row_id_global)||anyNA(target$row_id_global))stop('Invalid target row keys')
 if(!nrow(pred))return(rep(NA_real_,nrow(target)))
 if(anyDuplicated(pred$row_id_global)||anyNA(pred$row_id_global)||any(!pred$row_id_global%in%target$row_id_global))stop('Invalid prediction row keys')
 if(any(!is.finite(pred$z)))stop('Nonfinite prediction scores')
 ix<-match(pred$row_id_global,target$row_id_global)
 if('Subject_ID'%in%names(pred)&&any(as.character(pred$Subject_ID)!=as.character(target$Subject_ID[ix])))stop('Subject mismatch')
 out<-rep(NA_real_,nrow(target));out[ix]<-pred$z;out
}
choose_analysis_score<-function(diagnosis,oof,full,control='HC')ifelse(diagnosis==control,oof,full)
