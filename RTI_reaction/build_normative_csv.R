# Export normative results while preserving input score strings without rounding.
rti_build_normative_csv <- function(run, out) {
  need <- function(ok, text) if (!isTRUE(ok)) stop(text, call.=FALSE)
  hash <- function(path) digest::digest(file=path, algo='sha256')
  read <- function(path) read.csv(path, colClasses='character', na.strings=NULL,
    check.names=FALSE, fileEncoding='UTF-8-BOM', stringsAsFactors=FALSE)
  num <- function(x) suppressWarnings(as.numeric(x))
  flag <- function(x) ifelse(x, 'TRUE', 'FALSE')
  write <- function(x, path) {
    need(nrow(x)>0L, paste('Cannot export an empty table:', path))
    write.csv(x, path, row.names=FALSE, na='NA', fileEncoding='UTF-8')
  }
  run <- normalizePath(run, winslash='/', mustWork=TRUE)
  out <- normalizePath(out, winslash='/', mustWork=FALSE)
  need(out != run && !startsWith(out, paste0(run, '/')),
       'Output must be a separate analysis directory.')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  writeLines('BUILDING_NORMATIVE_CSV', file.path(out,'STATUS.txt'))
  manifest <- jsonlite::read_json(file.path(run,'manifest.json'), simplifyVector=TRUE)
  need(manifest$mode=='run', 'Only completed analysis output can be exported')
  need(startsWith(readLines(file.path(run,'STATUS.txt'),warn=FALSE)[1], 'NORMATIVE_COMPLETE'),
       'Normative model stage incomplete')
  paths <- c(file.path(run,c('normative_all_rows_scores.csv','normative_prediction_details.csv',
    'manifest.json','STATUS.txt','source/input.csv','source/frozen_subject_folds.csv')),
    list.files(file.path(run,'source'), pattern='\\.R$', full.names=TRUE))
  hashes <- setNames(vapply(paths,hash,character(1)),paths)
  need(hash(file.path(run,'source/input.csv'))==manifest$signature$input, 'Frozen input changed')
  need(hash(file.path(run,'source/frozen_subject_folds.csv'))==manifest$signature$folds, 'Frozen folds changed')
  for (name in names(manifest$signature$code))
    need(hash(file.path(run,'source',name))==manifest$signature$code[[name]], paste('Archived source changed:',name))
  raw <- read(file.path(run,'source/input.csv'))
  base <- read(file.path(run,'normative_all_rows_scores.csv'))
  details <- read(file.path(run,'normative_prediction_details.csv'))
  folds <- read(file.path(run,'source/frozen_subject_folds.csv'))
  need(!anyDuplicated(folds$Subject_ID) && setequal(num(folds$CV_fold),1:5),
       'Expected frozen five-fold assignment')
  order <- base$row_id_global
  need(identical(order, as.character(seq_len(nrow(raw)))), 'Raw row order or coverage changed')
  need(!anyDuplicated(base[c('Subject_ID','Timepoint')]), 'Duplicate subject visit')
  eligible <- details$row_id_global
  need(!anyDuplicated(eligible) && all(eligible %in% order), 'Normative prediction coverage invalid')
  need(all(base$model=='normative' & base$analysis_mode=='run'), 'Unexpected normative model')
  for (field in c('Subject_ID','Diagnosis','Platform','Gender'))
    need(identical(base[[field]], ifelse(raw[[field]] %in% c('','NA','NaN'),'NA',raw[[field]])),
         paste('Raw input identity mismatch:',field))
  need(all(num(base$Timepoint)==num(raw$Timepoint)), 'Visit identity mismatch')
  a <- num(base$NM_raw); b <- num(raw[['RTI_Median_five-choice_reaction_time']])
  need(all((is.finite(a)&is.finite(b)&a==b) | (!is.finite(a)&!is.finite(b))), 'Raw reaction time changed')
  fold_for <- function(ids) num(folds$CV_fold[match(ids,folds$Subject_ID)])
  index <- match(eligible,order)
  is_hc <- base$Diagnosis %in% 'HC'
  expected <- ifelse(is_hc[index],'HC_subject_OOF','patient_full_HC_new_subject')
  need(all(details$source==expected & base$normative_source[index]==expected), 'HC/patient score source mismatch')
  hc_detail <- is_hc[index]
  need(all(num(details$fold[hc_detail])==fold_for(base$Subject_ID[index[hc_detail]])), 'HC holdout fold mismatch')
  pool <- base[is_hc & order %in% eligible,,drop=FALSE]
  variables <- c('Age_baseline_years','Time_since_baseline_years','FSIQ')
  ranges <- list(); lookup <- list()
  key <- function(ref,platform,variable) paste(ref,platform,variable,sep='/')
  for (fold in 0:5) {
    ref <- if(fold==0L)'full_HC' else paste0('HC_fold_',fold,'_training')
    train <- if(fold==0L)pool else pool[fold_for(pool$Subject_ID)!=fold,,drop=FALSE]
    for (platform in c('All','RS6','Connect')) {
      selected <- if(platform=='All')train else train[train$Platform==platform,,drop=FALSE]
      need(nrow(selected)>0L, 'Empty HC reference training set')
      for (variable in variables) {
        values <- num(selected[[variable]])
        need(all(is.finite(values)), 'Nonfinite HC covariate')
        lo <- min(values); hi <- max(values); lookup[[key(ref,platform,variable)]] <- c(lo,hi)
        ranges[[length(ranges)+1L]] <- data.frame(model_role='normative', reference_set=ref,
          platform=platform,variable=variable,hc_min=lo,hc_max=hi,hc_rows=nrow(selected),
          hc_subjects=length(unique(selected$Subject_ID)), stringsAsFactors=FALSE)
      }
    }
  }
  rows <- base[,names(base)!='model',drop=FALSE]
  fold <- ifelse(is_hc,fold_for(base$Subject_ID),0L)
  need(!anyNA(fold[order %in% eligible]), 'HC absent from frozen fold map')
  ref <- ifelse(is.na(fold),'NA',ifelse(is_hc,paste0('HC_fold_',fold,'_training'),'full_HC'))
  rows$is_HC <- flag(is_hc); rows$is_baseline <- flag(num(base$Timepoint)==1)
  rows$HC_oof_fold <- ifelse(is_hc & !is.na(fold),as.character(fold),'NA')
  rows$covariate_reference_set <- ref; rows$normative_model <- 'normative'
  finite <- is.finite(num(base$normative_z))
  numeric_ok <- base$normative_numerical_ok=='TRUE'; mc_ok <- base$normative_posterior_mc_ok=='TRUE'
  scored <- order %in% eligible
  need(all(base$normative_numerical_ok[scored] %in% c('TRUE','FALSE') &
           base$normative_posterior_mc_ok[scored] %in% c('TRUE','FALSE')), 'Eligible normative lacks a precision flag')
  need(all(finite[scored]==(numeric_ok[scored]&mc_ok[scored])), 'Normative score/precision flag mismatch')
  need(all(numeric_ok[scored]&mc_ok[scored]), 'Unresolved normative precision prevents completion')
  mcse <- num(base$normative_posterior_mcse_z[scored])
  need(all(is.finite(mcse)&mcse>=0&mcse<=.03), 'Nonfinite or excessive normative MCSE cannot pass')
  need(all(base$normative_posterior_mc_status[scored]=='ok' & base$normative_status[scored]=='ok'),
       'Normative precision status mismatch')
  need(all(!finite[!scored] & base$normative_numerical_ok[!scored]=='NA'), 'Ineligible normative has a score')
  usable <- finite & numeric_ok & mc_ok; rows$normative_usable <- flag(usable)
  for (scope in c('global','platform')) {
    bad <- rep('NA',nrow(base))
    for (i in which(usable)) {
      platform <- if(scope=='platform')base$Platform[i] else 'All'
      outside <- character()
      for (variable in variables) {
        value <- num(base[[variable]][i]); limits <- lookup[[key(ref[i],platform,variable)]]
        need(is.finite(value), 'Nonfinite score covariate')
        if(value<limits[1])outside <- c(outside,paste0(variable,':below'))
        else if(value>limits[2])outside <- c(outside,paste0(variable,':above'))
      }
      bad[i] <- if(length(outside))paste(outside,collapse=';') else 'none'
    }
    rows[[paste0('normative_extrapolation_any_',scope)]] <- ifelse(usable,flag(bad!='none'),'NA')
    rows[[paste0('normative_extrapolation_variables_',scope)]] <- bad
    rows[[paste0('normative_within_reference_range_',scope)]] <- flag(usable & bad=='none')
  }
  write(rows,file.path(out,'RTI_analysis_ready.csv'))
  write(do.call(rbind,ranges),file.path(out,'HC_covariate_reference_ranges.csv'))
  write(data.frame(column=names(rows),source='Normative model or training-fold covariate flags'),
        file.path(out,'RTI_column_dictionary.csv'))
  write(data.frame(model='normative',score='normative',eligible_rows=length(eligible),usable_rows=length(eligible),
    numerical_unresolved=0,posterior_MC_unresolved=0,unresolved_scores=0),file.path(out,'final_score_quality.csv'))
  need(identical(rows,read(file.path(out,'RTI_analysis_ready.csv'))), 'CSV roundtrip changed scores')
  need(identical(hashes,setNames(vapply(paths,hash,character(1)),paths)), 'Source changed during export')
  outputs <- list.files(out,pattern='\\.csv$',full.names=TRUE)
  jsonlite::write_json(list(model='normative',run=basename(run),
    rows=nrow(rows),columns=ncol(rows),score_values_preserved_verbatim=TRUE,no_fit_or_rescoring=TRUE,
    input_sha256=as.list(hashes),output_sha256=as.list(setNames(vapply(outputs,hash,character(1)),basename(outputs))),
    R_version=as.character(getRversion())),file.path(out,'analysis_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
  writeLines('PASS: Normative, frozen input, HC subject folds, precision, support flags and CSV roundtrip.',
             file.path(out,'VALIDATION.txt'))
  writeLines('NORMATIVE_COMPLETE',file.path(out,'STATUS.txt'))
  invisible(file.path(out,'RTI_analysis_ready.csv'))
}
