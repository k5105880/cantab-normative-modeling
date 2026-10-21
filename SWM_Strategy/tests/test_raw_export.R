# Raw score cell text must survive CSV output, including excluded assessments.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/load.R'))
source(file.path(dirname(root),'R/input_data.R'))
source(file.path(dirname(root),'R/results.R'))
out<-tempfile('strategy_raw_export_');dir.create(out)
input<-file.path(out,'input.csv')
cells<-data.frame(Subject_ID=c('00001','00002'),Diagnosis='HC',Platform='RS6',Timepoint=1,
 Gender=c('F',NA),FSIQ=100,Age_baseline_years=c(45,55),Time_since_baseline_years=0,
 SWM_Strategy=c('32','14.000000000000002'))
other_outcomes<-c(SWM_Between_errors=4,SWM_Strategy=32,RVP_A_prime=.9,RVP_Median_latency=300,
  'RTI_Median_five-choice_reaction_time'=400,'RTI_Median_five-choice_movement_time'=300)
for(n in setdiff(names(other_outcomes),names(cells)))cells[[n]]<-unname(other_outcomes[n])
write.csv(cells,input,row.names=FALSE,na='')
raw<-read.csv(input,colClasses=c(Subject_ID='character'),check.names=FALSE,na.strings=c('','NA','NaN'))
x<-swm_adapt_raw(raw);x$CV_fold<-c(1L,NA_integer_)
stopifnot(identical(x$eligible,c(TRUE,FALSE)),raw$SWM_Strategy[2]!=14)
pred<-list(row_ids=1L,rows=data.frame(z_mid=.4,log_mid_p_lower=log(pnorm(.4)),
 log_mid_p_upper=log(pnorm(.4,lower.tail=FALSE)),posterior_precision_ok=TRUE,
 numerical_ok=TRUE,posterior_mcse_z=.01))
result<-swm_export(list(all=x,raw=raw,input=input,out=out),pred)
stopifnot(identical(result[names(raw)],raw),identical(result$normative_z,c(.4,NA_real_)))
native<-cantab_read_characters(file.path(out,'SWM_Strategy_scores.csv'))
stopifnot(identical(native$SWM_Strategy,cells$SWM_Strategy),identical(native$Subject_ID,cells$Subject_ID))
public<-cantab_export_results('SWM_Strategy',input,out)
exported<-cantab_read_characters(public)
stopifnot(identical(exported$raw_score,cells$SWM_Strategy),identical(exported$Subject_ID,cells$Subject_ID),
 identical(exported$score_status,c('available','not_available')))
# Even a tiny numeric change in the source must still fail the identity check.
cells$SWM_Strategy[2]<-'14';write.csv(cells,input,row.names=FALSE,na='')
failure<-tryCatch(swm_export(list(all=x,raw=raw,input=input,out=out),pred),error=identity)
stopifnot(inherits(failure,'error'))
cat('PASS: exact raw-score digits and leading-zero IDs survive native and strict shared exports; numeric model results unchanged.\n')
