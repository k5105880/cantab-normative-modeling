# Base-R preprocessing. Input is the raw-label-preserving CSV produced by
# prepare_raw.R. Numeric derivations use supplied years; raw source is untouched.
prepare_data <- function(input_csv, reference_dir, outdir) {
  d <- read.csv(input_csv, check.names=FALSE, stringsAsFactors=FALSE,
                colClasses='character', na.strings=c('', 'NA', 'NaN'))
  required <- c('Subject_ID','Diagnosis','Gender','Timepoint','Platform',
    'Age_baseline_years','Time_since_baseline_years','FSIQ','RVP_A_prime')
  if (!all(required %in% names(d))) stop('Missing required prepared input columns')
  if (!nrow(d) || anyDuplicated(names(d))) stop('Empty input or duplicate column names')
  key <- c('Subject_ID','Timepoint','Platform')
  if (anyNA(d[key]) || any(!nzchar(d$Subject_ID)) || anyDuplicated(d[key]))
    stop('Missing or duplicated identity key')
  if (anyNA(d$Diagnosis) || any(!d$Diagnosis %in% c('HC','BD','MDD')) ||
      any(!is.na(d$Gender) & !d$Gender %in% c('F','M')) ||
      any(!d$Platform %in% c('RS6','Connect'))) stop('Unrecognized input category')
  numeric_fields<-c('Timepoint','Age_baseline_years','Time_since_baseline_years','FSIQ','RVP_A_prime')
  if('Age_years'%in%names(d))numeric_fields<-c(numeric_fields,'Age_years')
  for (name in numeric_fields) {
    before <- d[[name]]
    d[[name]] <- suppressWarnings(as.numeric(before))
    if (any(!is.na(before) & !is.finite(d[[name]]))) stop('Invalid numeric field: ',name)
  }
  if (any(!is.finite(d$Timepoint) | d$Timepoint<1 | d$Timepoint!=floor(d$Timepoint)))
    stop('Timepoint must be a positive integer')
  if (any(d$RVP_A_prime<0 | d$RVP_A_prime>1,na.rm=TRUE)) stop('RVP_A_prime outside [0,1]')
  if (any(d$Time_since_baseline_years<0,na.rm=TRUE)) stop('Negative follow-up time')
  d$row_id_global <- seq_len(nrow(d))
  # Preserve the specified predictor arithmetic, including the year/day round-trip.
  if('Age_years'%in%names(d))d$Age <- d$Age_years*365.25
  d$Age_baseline <- d$Age_baseline_years*365.25
  d$Time_since_baseline <- d$Time_since_baseline_years*365.25
  if('Age_years'%in%names(d))d$Age_years <- d$Age/365.25
  d$Age_baseline_years <- d$Age_baseline/365.25
  d$Time_since_baseline_years <- d$Time_since_baseline/365.25
  d$Age_baseline_10y <- (d$Age_baseline_years-50)/10
  d$Time_since_baseline_1y <- d$Time_since_baseline_years
  d$FSIQ_15 <- (d$FSIQ-100)/15
  d$Platform <- factor(d$Platform,levels=c('RS6','Connect'))
  d$PlatformConnect_indicator <- as.numeric(d$Platform=='Connect')
  d$Gender <- factor(d$Gender,levels=c('F','M'))
  d$Subject_ID <- factor(d$Subject_ID)
  d$Diagnosis <- factor(d$Diagnosis)
  d$timepoint_f <- factor(d$Timepoint)
  d$Practice_order <- d$Timepoint-1
  bins <- c('0_baseline','1_second','2_third','3_fourth_or_later')
  d$Practice_bin <- factor(bins[pmin(d$Timepoint,4)],levels=bins)
  d$NM_raw <- d$RVP_A_prime
  d$NM_response <- d$RVP_A_prime
  eligible_cols <- c('NM_response','Age_baseline_10y','Time_since_baseline_1y',
    'Gender','FSIQ_15','Platform','Practice_bin','Subject_ID','Timepoint',
    'Time_since_baseline_years','Diagnosis')
  missing <- is.na(d[eligible_cols])
  d$model_eligible <- rowSums(missing)==0L
  d$model_exclusion_reason <- apply(missing,1,function(x)paste(eligible_cols[x],collapse=';'))

  hc_idx <- which(d$model_eligible & d$Diagnosis=='HC')
  fold_file<-if(!is.null(reference_dir))file.path(reference_dir,'primary_subject_folds.csv') else ''
  frozen<-nzchar(fold_file)&&file.exists(fold_file)
  if(frozen) {
    map <- read.csv(fold_file,stringsAsFactors=FALSE,colClasses='character')
    fold_audit<-list(source='supplied_frozen_reference',sha256=digest::digest(file=fold_file,algo='sha256'))
  } else {
    generated<-aprime_generate_folds(d[hc_idx,,drop=FALSE])
    map<-generated$assignment
    fold_audit<-generated$audit
  }
  if (!all(c('Subject_ID','CV_fold') %in% names(map)) || anyNA(map$Subject_ID) ||
      anyDuplicated(map$Subject_ID)) stop('Invalid frozen subject fold keys')
  original_fold <- map$CV_fold
  map$CV_fold <- suppressWarnings(as.integer(original_fold))
  if (anyNA(map$CV_fold) || any(map$CV_fold!=original_fold) ||
      !setequal(unique(map$CV_fold),1:5)) stop('Invalid frozen fold labels')
  if (!setequal(as.character(d$Subject_ID[hc_idx]),map$Subject_ID))
    stop('Frozen folds do not match the eligible HC subject cohort')
  d$CV_fold <- NA_integer_
  h <- which(d$Diagnosis=='HC')
  d$CV_fold[h] <- map$CV_fold[match(as.character(d$Subject_ID[h]),map$Subject_ID)]
  hc <- droplevels(d[hc_idx,,drop=FALSE])
  folds <- hc[c('row_id_global','Subject_ID','CV_fold')]
  if (anyNA(folds$CV_fold)) stop('An eligible HC row lacks a fold')

  result <- list(data=d,hc=hc,folds=folds)
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(map,file.path(outdir,'primary_subject_folds.csv'),row.names=FALSE)
  jsonlite::write_json(fold_audit,file.path(outdir,'fold_provenance.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  for(name in c('data','folds'))
    write.csv(result[[name]],file.path(outdir,paste0('prepared_',name,'.csv')),row.names=FALSE,na='')
  counts <- data.frame(metric=c('raw_rows','raw_subjects','eligible_rows','eligible_subjects',
    'eligible_HC_rows','eligible_HC_subjects'),
    value=c(nrow(d),nlevels(d$Subject_ID),sum(d$model_eligible),
      length(unique(d$Subject_ID[d$model_eligible])),nrow(hc),nlevels(hc$Subject_ID)))
  write.csv(counts,file.path(outdir,'prepared_counts.csv'),row.names=FALSE)
  tmp <- tempfile('prepared_',tmpdir=outdir,fileext='.rds')
  saveRDS(result,tmp)
  if (!file.rename(tmp,file.path(outdir,'prepared.rds'))) stop('Unable to save prepared data atomically')
  result
}
