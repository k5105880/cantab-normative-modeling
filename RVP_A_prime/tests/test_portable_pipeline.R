args<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
for(f in c('prepare_raw.R','prepare_data.R','folds.R','portable_pipeline.R','precision_helpers.R','calibration.R'))source(file.path(root,'R',f))
source(file.path(root,'tests/synthetic_input.R'))
# Sourcing only defines the public entry point and does not create a run or mutate options.
before<-list.files(root,all.files=TRUE);env<-Sys.getenv()
source(file.path(root,'run.R'))
stopifnot(is.function(run_outcome),identical(before,list.files(root,all.files=TRUE)),identical(env,Sys.getenv()))
stopifnot(identical(aprime_cdf_map(as.list(1:3),function(x)x^2,3L,os_type='windows'),as.list(c(1,4,9))))
tmp<-tempfile('aprime_assembly_');dir.create(tmp);dir.create(file.path(tmp,'results'))
raw<-aprime_test_input();raw$Subject_ID<-paste0('00',raw$Subject_ID);raw$RVP_A_prime[1]<-NA_real_
raw$row_id_global<-1000L+seq_len(nrow(raw))
write.csv(raw,file.path(tmp,'raw.csv'),row.names=FALSE,na='')
aprime_prepare_raw(file.path(tmp,'raw.csv'),file.path(tmp,'input'))
p<-prepare_data(file.path(tmp,'input/aligned_raw.csv'),NULL,file.path(tmp,'input'))
idx<-which(p$data$model_eligible)
scores<-data.frame(row_id_global=idx,role=ifelse(p$data$Diagnosis[idx]=='HC',paste0('normative_OOF_fold',p$data$CV_fold[idx]),'normative_patient'),
 z=seq_along(idx)/1000,precision_status='ok',quadrature_z_difference=.001,posterior_mcse_z=.01)
write.csv(scores,file.path(tmp,'results/normative_score_precision.csv'),row.names=FALSE)
assembled<-aprime_assemble(tmp)
stopifnot(identical(assembled$Subject_ID,raw$Subject_ID),assembled$normative_z[1]=='',
  identical(assembled$row_id_global,as.character(seq_len(nrow(raw)))),
  identical(assembled$precision_status,c('not_scored',rep('ok',247))),
  identical(aprime_read_cells(file.path(tmp,'results/RVP_A_prime_scores.csv')),aprime_assembled_rows(tmp)))
scores$role[1]<-'normative_patient';write.csv(scores,file.path(tmp,'results/normative_score_precision.csv'),row.names=FALSE)
stopifnot(inherits(try(aprime_assemble(tmp),silent=TRUE),'try-error'))
writeLines('test-identity',file.path(tmp,'identity.sha256'))
aprime_seal(tmp,'prepare',file.path(tmp,'raw.csv'));aprime_verify_seals(tmp,'prepare')
cat('\n',file=file.path(tmp,'raw.csv'),append=TRUE)
stopifnot(inherits(try(aprime_verify_seals(tmp,'prepare'),silent=TRUE),'try-error'))
# An ineligible row is retained even when its platform has no eligible HC support.
unsupported<-file.path(tmp,'unsupported platform');dir.create(file.path(unsupported,'results'),recursive=TRUE)
single_platform<-aprime_test_input();single_platform$Platform<-'RS6'
last<-nrow(single_platform)
single_platform$Platform[last]<-'Connect';single_platform$RVP_A_prime[last]<-NA_real_
write.csv(single_platform,file.path(unsupported,'raw.csv'),row.names=FALSE,na='')
aprime_prepare_raw(file.path(unsupported,'raw.csv'),file.path(unsupported,'input'))
single<-prepare_data(file.path(unsupported,'input/aligned_raw.csv'),NULL,file.path(unsupported,'input'))
idx<-which(single$data$model_eligible)
accepted<-data.frame(row_id_global=idx,
  role=ifelse(single$data$Diagnosis[idx]=='HC',
    paste0('normative_OOF_fold',single$data$CV_fold[idx]),'normative_patient'),
  z=seq_along(idx)/1000,precision_status='ok',quadrature_z_difference=.001,posterior_mcse_z=.01)
write.csv(accepted,file.path(unsupported,'results/normative_score_precision.csv'),row.names=FALSE)
retained<-aprime_assemble(unsupported)
stopifnot(nrow(retained)==248L,retained$normative_z[last]=='',
  retained$precision_status[last]=='not_scored',retained$outside_hc_support[last]=='TRUE',
  all(retained$outside_hc_support[retained$Diagnosis=='HC']=='FALSE'),
  identical(as.numeric(retained$normative_z[idx]),accepted$z))
# A scored row still requires HC covariate support on its own platform.
single_platform$RVP_A_prime[last]<-.8
write.csv(single_platform,file.path(unsupported,'raw.csv'),row.names=FALSE,na='')
aprime_prepare_raw(file.path(unsupported,'raw.csv'),file.path(unsupported,'input'))
single<-prepare_data(file.path(unsupported,'input/aligned_raw.csv'),NULL,file.path(unsupported,'input'))
accepted<-rbind(accepted,data.frame(row_id_global=last,role='normative_patient',
  z=.2,precision_status='ok',quadrature_z_difference=.001,posterior_mcse_z=.01))
write.csv(accepted,file.path(unsupported,'results/normative_score_precision.csv'),row.names=FALSE)
err<-tryCatch(aprime_assemble(unsupported),error=identity)
stopifnot(inherits(err,'error'),
  grepl('No eligible HC covariate support for platform Connect',conditionMessage(err),fixed=TRUE))
# Numeric-looking participant identifiers stay distinct bootstrap clusters.
cal<-file.path(tmp,'calibration');dir.create(file.path(cal,'results'),recursive=TRUE)
write.csv(data.frame(row_id_global=1:4,Subject_ID=rep(c('001','1'),each=2),
  Diagnosis='HC',Platform='RS6',Practice_bin=rep(c('0_baseline','1_second'),2),
  Time_since_baseline_years=rep(c(0,1),2)),
  file.path(cal,'results/RVP_A_prime_scores.csv'),row.names=FALSE)
write.csv(data.frame(row_id_global=1:4,role='normative_OOF_fold1',
  p_lower=c(.2,.3,.7,.8),p_upper=c(.2,.3,.7,.8),
  log_lower=log(c(.2,.3,.7,.8)),log_mass=-Inf,log_upper=log(c(.8,.7,.3,.2))),
  file.path(cal,'results/normative_score_precision.csv'),row.names=FALSE)
calibrate_scores(cal,list(seed=20260922L,calibration_bootstraps=32L))
summary<-read.csv(file.path(cal,'results/OOF_calibration.csv'))
stopifnot(summary$n_subjects[summary$stratum=='overall']==2L,all(summary$n_subjects==2L))
unlink(tmp,recursive=TRUE)
cat('PASS: side-effect-free source, Windows serial fallback, raw-cell/score assembly, unsupported-platform missing-row retention and scored-row rejection, role and artifact tamper rejection, distinct numeric-looking participant calibration clusters.\n')
