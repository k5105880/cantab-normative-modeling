# Synthetic fixture only. No fitting, private input or Python is required.
args <- gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE)
here <- dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
source(file.path(here,'build_normative_csv.R'))
root <- tempfile('reaction_export_test_');dir.create(root)
run <- file.path(root,'model');dir.create(file.path(run,'source'),recursive=TRUE)
raw <- expand.grid(Timepoint=1:2,subject=1:12)
raw$Subject_ID <- sprintf('synthetic_%02d',raw$subject)
raw$Diagnosis <- ifelse(raw$subject<=10,'HC','MDD')
raw$Platform <- ifelse(raw$subject<=5 | raw$subject==11,'RS6','Connect')
raw$Gender <- 'F';raw$FSIQ <- '100'
raw$Age_baseline_years <- ifelse(raw$subject<=10,'50','80')
raw$Time_since_baseline_years <- as.character(raw$Timepoint-1)
raw[['RTI_Median_five-choice_reaction_time']] <- as.character(300+raw$Timepoint)
raw[['RTI_Median_five-choice_reaction_time']][1] <- '350.12345678901234'
raw[['RTI_Median_five-choice_movement_time']] <- '250'
raw$RVP_A_prime <- '0.8';raw$RVP_Median_latency <- '400'
raw$SWM_Between_errors <- '4';raw$SWM_Strategy <- '5'
raw$subject <- NULL
raw[nrow(raw),'RTI_Median_five-choice_reaction_time'] <- 'NA'
write.csv(raw,file.path(run,'source/input.csv'),row.names=FALSE)
folds <- data.frame(Subject_ID=sprintf('synthetic_%02d',1:10),CV_fold=rep(1:5,2))
write.csv(folds,file.path(run,'source/frozen_subject_folds.csv'),row.names=FALSE)
base <- raw
base$row_id_global <- as.character(seq_len(nrow(base)))
base$NM_raw <- base[['RTI_Median_five-choice_reaction_time']]
base$model <- 'normative';base$analysis_mode <- 'run'
base$normative_z <- '0.25000000000000001';base$normative_z_estimate <- base$normative_z
base$normative_source <- ifelse(base$Diagnosis=='HC','HC_subject_OOF','patient_full_HC_new_subject')
base$normative_numerical_ok <- 'TRUE';base$normative_posterior_mc_ok <- 'TRUE'
base$normative_posterior_mcse_z <- '0.01'
base$normative_posterior_mc_status <- 'ok';base$normative_status <- 'ok'
for (nm in grep('^normative_',names(base),value=TRUE))base[nrow(base),nm] <- 'NA'
base[nrow(base),'normative_status'] <- 'missing_response_or_covariate'
write.csv(base,file.path(run,'normative_all_rows_scores.csv'),row.names=FALSE)
details <- data.frame(row_id_global=base$row_id_global[-nrow(base)],
  source=base$normative_source[-nrow(base)],
  fold=folds$CV_fold[match(base$Subject_ID[-nrow(base)],folds$Subject_ID)])
write.csv(details,file.path(run,'normative_prediction_details.csv'),row.names=FALSE)
hash <- function(path)digest::digest(file=path,algo='sha256')
jsonlite::write_json(list(mode='run',signature=list(config=list(version='reaction_normative_v1'),
  input=hash(file.path(run,'source/input.csv')),folds=hash(file.path(run,'source/frozen_subject_folds.csv')))),
  file.path(run,'manifest.json'),auto_unbox=TRUE)
writeLines('NORMATIVE_COMPLETE',file.path(run,'STATUS.txt'))
out <- file.path(root,'analysis')
rti_build_normative_csv(run,out)
result <- read.csv(file.path(out,'RTI_analysis_ready.csv'),colClasses='character',na.strings=NULL,check.names=FALSE)
stopifnot(identical(result$normative_z,base$normative_z),
  identical(result$NM_raw,base$NM_raw),
  identical(result[['RTI_Median_five-choice_reaction_time']],
            raw[['RTI_Median_five-choice_reaction_time']]),
  all(result$normative_usable[-nrow(result)]=='TRUE'),
  result$normative_usable[nrow(result)]=='FALSE',
  result$normative_extrapolation_any_global[21]=='TRUE',
  result$normative_within_reference_range_global[21]=='FALSE',
  result$HC_oof_fold[21]=='NA')
source(file.path(dirname(here),'R','input_data.R'))
source(file.path(dirname(here),'R','results.R'))
cantab_export_results('RTI_reaction',file.path(run,'source','input.csv'),root)
public <- read.csv(file.path(root,'normative_scores.csv'),colClasses='character',
                   na.strings=NULL,check.names=FALSE)
stopifnot(identical(public$normative_z,
  ifelse(base$normative_z=='NA','',base$normative_z)),
  identical(public$raw_score,ifelse(base$NM_raw=='NA','',base$NM_raw)),
  identical(public$Subject_ID,base$Subject_ID))
# Every visit can be ineligible for one HC participant. Only eligible HC
# participants belong in the frozen folds; the unscored row must still export.
extra_raw <- raw[nrow(raw),,drop=FALSE]
extra_raw$Subject_ID <- 'synthetic_ineligible_HC';extra_raw$Diagnosis <- 'HC'
extra_base <- base[nrow(base),,drop=FALSE]
extra_base$Subject_ID <- extra_raw$Subject_ID;extra_base$Diagnosis <- 'HC'
extra_base$row_id_global <- as.character(nrow(base)+1L)
raw_with_ineligible <- rbind(raw,extra_raw)
base_with_ineligible <- rbind(base,extra_base)
refresh_fixture <- function(input_rows,score_rows) {
  write.csv(input_rows,file.path(run,'source/input.csv'),row.names=FALSE)
  write.csv(score_rows,file.path(run,'normative_all_rows_scores.csv'),row.names=FALSE)
  jsonlite::write_json(list(mode='run',signature=list(config=list(version='reaction_normative_v1'),
    input=hash(file.path(run,'source/input.csv')),folds=hash(file.path(run,'source/frozen_subject_folds.csv')))),
    file.path(run,'manifest.json'),auto_unbox=TRUE)
}
refresh_fixture(raw_with_ineligible,base_with_ineligible)
ineligible_out <- file.path(root,'ineligible_hc_analysis')
rti_build_normative_csv(run,ineligible_out)
ineligible_result <- read.csv(file.path(ineligible_out,'RTI_analysis_ready.csv'),
  colClasses='character',na.strings=NULL,check.names=FALSE)
stopifnot(nrow(ineligible_result)==25L,
  ineligible_result$normative_z[25]=='NA',ineligible_result$normative_usable[25]=='FALSE',
  ineligible_result$HC_oof_fold[25]=='NA',ineligible_result$covariate_reference_set[25]=='NA',
  identical(ineligible_result$normative_z[1:24],base$normative_z))
# All documented missing markers for Gender identify the same excluded value.
for(marker in c('', 'NA', 'NaN')) {
  missing_raw <- raw;missing_base <- base
  missing_raw$Gender[nrow(raw)] <- marker;missing_base$Gender[nrow(base)] <- 'NA'
  refresh_fixture(missing_raw,missing_base)
  missing_out <- tempfile('missing_gender_',tmpdir=root)
  rti_build_normative_csv(run,missing_out)
  missing_result <- read.csv(file.path(missing_out,'RTI_analysis_ready.csv'),
    colClasses='character',na.strings=NULL,check.names=FALSE)
  stopifnot(identical(missing_result$normative_z,base$normative_z),
    missing_result$Gender[nrow(base)]=='NA',missing_result$normative_usable[nrow(base)]=='FALSE')
}
refresh_fixture(raw,base)
base$normative_posterior_mcse_z[1] <- 'NA'
write.csv(base,file.path(run,'normative_all_rows_scores.csv'),row.names=FALSE)
bad <- file.path(root,'bad_analysis')
err <- tryCatch(rti_build_normative_csv(run,bad),error=identity)
stopifnot(inherits(err,'error'),grepl('MCSE',conditionMessage(err)),
  readLines(file.path(bad,'STATUS.txt'))!='NORMATIVE_COMPLETE')
unlink(root,recursive=TRUE)
cat('PASS: exact intermediate/public score and long-decimal raw strings, participant sources, holdout support, ineligible rows, MCSE rejection; no MCMC.\n')
