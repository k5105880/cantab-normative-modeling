args<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
source(file.path(root,'R/prepare_raw.R'));source(file.path(root,'R/prepare_data.R'))
source(file.path(root,'R/folds.R'));source(file.path(root,'tests/synthetic_input.R'))
tmp<-tempfile('aprime_prepare_');dir.create(tmp)
raw<-aprime_test_input();raw$RVP_A_prime[1]<-NA_real_;raw$FSIQ[5]<-NA_real_
raw$Subject_ID<-paste0('00',raw$Subject_ID)
raw$Diagnosis2<-rep(c('subtype A','','NA','NaN'),length.out=nrow(raw))
raw$source_note<-rep(c('comma, retained','quote " retained','  padded  ','NA'),length.out=nrow(raw))
path<-file.path(tmp,'raw.csv');write.csv(raw,path,row.names=FALSE,na='')
cells<-aprime_read_cells(path)
raw<-read.csv(path,stringsAsFactors=FALSE)
aligned<-aprime_prepare_raw(path,file.path(tmp,'input'))
stopifnot(identical(aligned,cells),
  identical(aprime_read_cells(file.path(tmp,'input/aligned_raw.csv')),cells),
  setequal(list.files(file.path(tmp,'input')),c('aligned_raw.csv','source_manifest.json')))
manifest<-jsonlite::fromJSON(file.path(tmp,'input/source_manifest.json'))
stopifnot(identical(manifest$input_sha256,digest::digest(file=path,algo='sha256')),
  manifest$raw_rows==248L,manifest$raw_subjects==62L,
  identical(manifest$diagnosis_policy,'preserve_input_labels'))
p<-prepare_data(file.path(tmp,'input/aligned_raw.csv'),NULL,file.path(tmp,'prepared'))
stopifnot(nrow(p$data)==248L,sum(p$data$model_eligible)==246L,nrow(p$hc)==238L,
  !p$data$model_eligible[1],!p$data$model_eligible[5],
  identical(p$data$row_id_global,seq_len(248L)),
  identical(p$data$NM_response,raw$RVP_A_prime),
  all(p$data$Age_baseline_10y==(raw$Age_baseline_years*365.25/365.25-50)/10),
  all(p$data$Time_since_baseline_1y==raw$Time_since_baseline_years*365.25/365.25),
  !any(c('Age_years','Onset_age','illness_duration','HRSD17')%in%names(p$data)),
  !anyDuplicated(unique(p$folds[c('Subject_ID','CV_fold')])$Subject_ID))
frozen<-file.path(tmp,'reference');dir.create(frozen)
file.copy(file.path(tmp,'prepared/primary_subject_folds.csv'),file.path(frozen,'primary_subject_folds.csv'))
q<-prepare_data(file.path(tmp,'input/aligned_raw.csv'),frozen,file.path(tmp,'frozen'))
stopifnot(identical(p,q))
set.seed(29);before<-.Random.seed
fold1<-aprime_generate_folds(p$hc);stopifnot(identical(before,.Random.seed))
fold2<-aprime_generate_folds(p$hc);stopifnot(identical(fold1,fold2),
  fold1$audit$source=='generated_for_input')
with_unused_levels<-p$hc
with_unused_levels$Subject_ID<-factor(as.character(p$hc$Subject_ID),
  levels=c('A_patient_only',levels(p$hc$Subject_ID),'Z_patient_only'))
stopifnot(identical(aprime_generate_folds(with_unused_levels),fold1))
reject<-function(x){f<-file.path(tmp,'invalid.csv');write.csv(x,f,row.names=FALSE,na='');
  stopifnot(inherits(try(aprime_prepare_raw(f,file.path(tmp,'invalid')),silent=TRUE),'try-error'))}
bad<-raw;bad$Subject_ID[2]<-bad$Subject_ID[1];bad$Timepoint[2]<-bad$Timepoint[1];bad$Platform[2]<-bad$Platform[1];reject(bad)
bad<-raw;bad$RVP_A_prime[2]<-1.1;reject(bad)
bad<-raw;bad$Age_baseline_years[2]<-Inf;reject(bad)
bad<-raw;bad$Gender[2]<-'unknown';reject(bad)
bad<-raw;bad$FSIQ<-NULL;reject(bad)
bad<-raw;bad$Subject_ID[1]<-' ';reject(bad)
bad<-raw;bad$Timepoint[2]<-1.5;reject(bad)
# A frozen map may never silently assign a changed HC cohort.
map<-read.csv(file.path(frozen,'primary_subject_folds.csv'));map$Subject_ID[1]<-'unknown'
write.csv(map,file.path(frozen,'primary_subject_folds.csv'),row.names=FALSE)
stopifnot(inherits(try(prepare_data(file.path(tmp,'input/aligned_raw.csv'),frozen,file.path(tmp,'bad_folds')),silent=TRUE),'try-error'))
unlink(tmp,recursive=TRUE)
cat('PASS: minimal normative input, missingness, grouped fold reproducibility, frozen fold identity, and malformed input rejection.\n')
