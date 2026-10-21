# Raw count spelling and participant identifiers survive native/shared exports.
a<-commandArgs(FALSE);root<-dirname(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',a[grepl('^--file=',a)][1]),fixed=TRUE))))
source(file.path(root,'R/pipeline.R'))
source(file.path(dirname(root),'R/input_data.R'))
source(file.path(dirname(root),'R/results.R'))
out<-tempfile('between_raw_export_');dir.create(file.path(out,'private'),recursive=TRUE)
input<-file.path(out,'input.csv')
cells<-data.frame(Subject_ID=c('00001','00002'),Diagnosis='HC',Platform='RS6',Timepoint=1,
 Gender=c('F',NA),FSIQ=100,Age_baseline_years=c(45,55),Time_since_baseline_years=0,
 SWM_Between_errors=c('0004','1.4000000000000000e1'))
other_outcomes<-c(SWM_Between_errors=4,SWM_Strategy=32,RVP_A_prime=.9,RVP_Median_latency=300,
  'RTI_Median_five-choice_reaction_time'=400,'RTI_Median_five-choice_movement_time'=300)
for(n in setdiff(names(other_outcomes),names(cells)))cells[[n]]<-unname(other_outcomes[n])
write.csv(cells,input,row.names=FALSE,na='')
raw<-read.csv(input,colClasses=c(Subject_ID='character'),check.names=FALSE,na.strings=c('','NA','NaN'))
scores<-raw;scores$normative_z<-c(.4,NA_real_)
scores$normative_score_source<-'HC_subject_held_out'
scores$posterior_precision_ok<-scores$numerical_ok<-c(TRUE,NA)
scores$posterior_mcse_z<-c(.01,NA_real_)
stopifnot(identical(write_score_table(scores,list(input=input,run=out)),scores),
 identical(scores[names(raw)],raw))
native<-cantab_read_characters(file.path(out,'private/SWM_scores.csv'))
stopifnot(identical(native$SWM_Between_errors,cells$SWM_Between_errors),
 identical(native$Subject_ID,cells$Subject_ID))
public<-cantab_export_results('SWM_Between_errors',input,out)
exported<-cantab_read_characters(public)
stopifnot(identical(exported$raw_score,cells$SWM_Between_errors),identical(exported$Subject_ID,cells$Subject_ID),
 identical(exported$score_status,c('available','not_available')))
# Changed values cannot be disguised by replacing the CSV spelling.
cells$SWM_Between_errors[1]<-'5';write.csv(cells,input,row.names=FALSE,na='')
failure<-tryCatch(write_score_table(scores,list(input=input,run=out)),error=identity)
stopifnot(inherits(failure,'error'))
cat('PASS: raw count cells and leading-zero IDs survive strict native/shared export; changed raw input is refused.\n')
