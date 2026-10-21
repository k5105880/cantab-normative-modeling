aprime_read_cells <- function(path) {
  # A parser warning can otherwise leave a plausible but incomplete table.
  withCallingHandlers(
    invisible(scan(path,what=character(),sep=',',quote='"',comment.char='',
      blank.lines.skip=FALSE,na.strings=NULL,fileEncoding='UTF-8-BOM',quiet=TRUE)),
    warning=function(w)stop('Malformed CSV: ',conditionMessage(w),call.=FALSE))
  counts<-count.fields(path,sep=',',quote='"',comment.char='',blank.lines.skip=FALSE)
  counts<-counts[!is.na(counts)]
  if(!length(counts)||length(unique(counts))!=1L)stop('Malformed CSV row width')
  d<-read.csv(path,colClasses='character',check.names=FALSE,na.strings=NULL,
    stringsAsFactors=FALSE,blank.lines.skip=FALSE,fileEncoding='UTF-8-BOM')
  if(!nrow(d)||anyDuplicated(names(d))||any(!nzchar(names(d))))stop('Empty CSV or duplicate/missing column names')
  d
}

aprime_validate_raw <- function(d) {
  required<-c('Subject_ID','Timepoint','Platform','Diagnosis','Gender',
    'Age_baseline_years','Time_since_baseline_years','FSIQ','RVP_A_prime')
  if(!all(required%in%names(d)))stop('Missing required raw columns: ',paste(setdiff(required,names(d)),collapse=', '))
  missing<-function(x)is.na(x)|x%in%c('','NA','NaN')
  key<-c('Subject_ID','Timepoint','Platform')
  if(any(vapply(d[key],function(x)any(missing(x)|!nzchar(trimws(x))),logical(1)))||anyDuplicated(d[key]))
    stop('Missing or duplicated subject/visit/platform identity key')
  if(any(!d$Diagnosis%in%c('HC','BD','MDD'))||any(!d$Platform%in%c('RS6','Connect'))||
     any(!missing(d$Gender)&!d$Gender%in%c('F','M')))stop('Unrecognized diagnosis, gender or platform')
  numeric_fields<-c('Timepoint','Age_baseline_years','Time_since_baseline_years','FSIQ','RVP_A_prime')
  if('Age_years'%in%names(d))numeric_fields<-c(numeric_fields,'Age_years')
  numbers<-lapply(d[numeric_fields],function(x){x[missing(x)]<-NA_character_;suppressWarnings(as.numeric(x))})
  for(nm in numeric_fields)if(any(!missing(d[[nm]])&!is.finite(numbers[[nm]])))stop('Invalid numeric input in ',nm)
  tp<-numbers$Timepoint
  if(any(!is.finite(tp)|tp<1|tp!=floor(tp)))stop('Timepoint must be a positive integer')
  # A numeric duplicate such as "1" and "1.0" is also one visit.
  numeric_key<-d[key];numeric_key$Timepoint<-tp
  if(anyDuplicated(numeric_key))stop('Duplicated numeric subject/visit/platform identity key')
  if(any(numbers$RVP_A_prime<0|numbers$RVP_A_prime>1,na.rm=TRUE))stop('RVP_A_prime outside [0,1]')
  for(nm in intersect(c('Age_years','Age_baseline_years','Time_since_baseline_years'),names(numbers)))
    if(any(numbers[[nm]]<0,na.rm=TRUE))stop('Negative age or follow-up time')
  for(day in c('Age','Age_baseline','Time_since_baseline')) {
    year<-paste0(day,'_years')
    if(!day%in%names(d)||!year%in%names(numbers))next
    x<-d[[day]];x[missing(x)]<-NA_character_;x<-suppressWarnings(as.numeric(x))
    if(any(!missing(d[[day]])&!is.finite(x))||
       any(is.na(x)!=is.na(numbers[[year]]))||
       any(abs(x-numbers[[year]]*365.25)>.01,na.rm=TRUE))stop('Conflicting day/year unit fields: ',day)
  }
  invisible(TRUE)
}

aprime_prepare_raw <- function(input_file,outdir) {
  raw<-aprime_read_cells(input_file);aprime_validate_raw(raw)
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(raw,file.path(outdir,'aligned_raw.csv'),row.names=FALSE,na='')
  manifest<-list(schema_version=4L,diagnosis_policy='preserve_input_labels',
    input_sha256=digest::digest(file=input_file,algo='sha256'),
    raw_rows=nrow(raw),raw_subjects=length(unique(raw$Subject_ID)),
    row_order='original raw input order',
    age_time_input_unit='years')
  jsonlite::write_json(manifest,file.path(outdir,'source_manifest.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  invisible(raw)
}
