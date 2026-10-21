# Reject malformed CSVs before a parser can silently lose rows or invent headers.
local({
  source('R/input_data.R', local=TRUE)
  source('RVP_A_prime/R/prepare_raw.R', local=TRUE)
  work <- tempfile('csv structure ');dir.create(work)
  on.exit(unlink(work,recursive=TRUE),add=TRUE)
  raw <- data.frame(Subject_ID=sprintf('000%02d',1:12),Diagnosis='HC',
    Timepoint=1,Platform='Connect',Gender='F',FSIQ=100,
    Age_baseline_years=40,Time_since_baseline_years=0,
    SWM_Between_errors=4,SWM_Strategy=5,RVP_A_prime='0.90000000000001234',
    RVP_Median_latency=400,check.names=FALSE)
  raw[['RTI_Median_five-choice_reaction_time']] <- 320
  raw[['RTI_Median_five-choice_movement_time']] <- 250
  raw$Notes <- 'ordinary note'
  plain <- file.path(work,'plain.csv')
  write.csv(raw,plain,quote=FALSE,row.names=FALSE)
  lines <- readLines(plain)
  bad <- lines;bad[7] <- sub('ordinary note$','"unterminated note',bad[7])
  writeLines(bad,file.path(work,'unclosed_quote.csv'))
  bad <- lines;bad[-1] <- paste0('unlabelled_value,',bad[-1])
  writeLines(bad,file.path(work,'extra_fields.csv'))
  bad <- lines;bad[7] <- sub(',ordinary note$','',bad[7])
  writeLines(bad,file.path(work,'missing_field.csv'))
  writeLines(append(lines,'',after=6),file.path(work,'blank_record.csv'))
  readers <- list(shared=cantab_read_characters,validated=read_cantab_input,
                  aprime=aprime_read_cells)
  for(name in names(readers)) {
    for(case in c('unclosed_quote','extra_fields','missing_field','blank_record')) {
      result <- tryCatch(suppressWarnings(readers[[name]](file.path(work,paste0(case,'.csv')))),error=identity)
      if(!inherits(result,'error'))stop('Malformed CSV accepted by ',name,': ',case,
        '; returned ',nrow(result),' of 12 source rows')
    }
  }
  # Legal quoted commas, doubled quotes, newlines, Unicode and backslashes survive.
  raw$Notes <- rep(c('comma, retained','quote " retained','line one\nline two',
    '日本語の記録','literal \\ path',''),2)
  canonical <- file.path(work,'quoted.csv')
  write.csv(raw,canonical,row.names=FALSE,na='',fileEncoding='UTF-8')
  bytes <- readBin(canonical,'raw',n=file.info(canonical)$size)
  writeBin(c(as.raw(c(239,187,191)),bytes),file.path(work,'bom.csv'))
  text <- rawToChar(bytes)
  writeBin(charToRaw(gsub('\n','\r\n',text,fixed=TRUE)),file.path(work,'crlf.csv'))
  writeBin(bytes[-length(bytes)],file.path(work,'no_final_newline.csv'))
  for(name in names(readers))for(file in c('quoted.csv','bom.csv','crlf.csv','no_final_newline.csv')) {
    result <- readers[[name]](file.path(work,file))
    stopifnot(nrow(result)==12L,identical(result$Subject_ID,raw$Subject_ID),
      identical(result$Notes[nzchar(raw$Notes)],raw$Notes[nzchar(raw$Notes)]))
    if(name!='validated')stopifnot(identical(result$RVP_A_prime,raw$RVP_A_prime))
  }
  # Input-table construction and analysis preparation must fail before writing outputs.
  p <- data.frame(Subject_ID='00001',Diagnosis='HC',Gender='F',FSIQ=100,
    Age_baseline_years=40,Baseline_date='2020-01-01')
  a <- raw[1,c('Subject_ID','Timepoint','Platform',cantab_outcome_columns)]
  a$Assessment_date <- '2020-01-01'
  master <- file.path(work,'participants.csv');write.csv(p,master,row.names=FALSE)
  bad <- readLines(master);bad[2] <- paste0('extra,',bad[2]);writeLines(bad,master)
  dest <- file.path(work,'created.csv')
  result <- tryCatch(create_cantab_input(master,a,dest),error=identity)
  stopifnot(inherits(result,'error'),!file.exists(dest))
  out <- file.path(work,'analysis')
  result <- tryCatch(cantab_prepare_input(file.path(work,'unclosed_quote.csv'),out,
    outcome='RVP_A_prime'),error=identity)
  stopifnot(inherits(result,'error'),!file.exists(file.path(out,'analysis_input.csv')))
  cat('PASS: malformed CSV rejection, no partial input creation, and valid quoting/UTF-8/line-ending preservation.\n')
})
