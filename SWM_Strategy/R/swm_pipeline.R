swm_write_csv <- function(x,out,name) write.csv(x,file.path(out,name),row.names=FALSE,na='')

swm_adapt_raw <- function(raw) {
  required <- c('Subject_ID','Diagnosis','Platform','Timepoint','Gender','FSIQ','SWM_Strategy')
  if (!all(required %in% names(raw))) stop('Missing raw columns: ',paste(setdiff(required,names(raw)),collapse=', '))
  if (anyDuplicated(names(raw)) || any(c('normative_z','normative_model','normative_score_source')%in%names(raw))) stop('Use raw input, without duplicate or final-score columns')
  x <- raw
  if (!'row_id_global' %in% names(x)) x$row_id_global <- seq_len(nrow(x))
  if(anyNA(x$row_id_global)||anyDuplicated(x$row_id_global)||anyNA(x$Subject_ID)||any(!nzchar(as.character(x$Subject_ID))))stop('Invalid row/participant identifiers')
  if(anyNA(x$Platform)||any(!x$Platform %in% c('RS6','Connect')))stop('Unknown platform')
  if(anyNA(x$Diagnosis)||any(!nzchar(as.character(x$Diagnosis))))stop('Missing diagnosis label')
  if(any(vapply(split(x$Diagnosis,x$Subject_ID),function(a)length(unique(a))!=1L,logical(1))))stop('Diagnosis changes within participant')
  to_num <- function(a) suppressWarnings(as.numeric(as.character(a)))
  years <- function(year_col,day_col) {
    if (!year_col %in% names(raw) && !day_col %in% names(raw)) stop('Missing ',year_col,' (or ',day_col,')')
    y <- if(year_col %in% names(raw)) to_num(raw[[year_col]]) else to_num(raw[[day_col]])/365.25
    if(all(c(year_col,day_col)%in%names(raw))) {
      days_as_years<-to_num(raw[[day_col]])/365.25
      if(any(is.na(y)!=is.na(days_as_years))||any(abs(y-days_as_years)>1e-6,na.rm=TRUE))stop('Conflicting day/year columns: ',year_col)
    }
    y
  }
  x$Age_baseline_years <- years('Age_baseline_years','Age_baseline')
  x$Time_since_baseline_years <- years('Time_since_baseline_years','Time_since_baseline')
  x$Age_baseline_10y <- (x$Age_baseline_years-50)/10
  x$Time_since_baseline_1y <- x$Time_since_baseline_years
  for(n in c('Timepoint','FSIQ','SWM_Strategy'))x[[n]]<-to_num(x[[n]])
  x$FSIQ_15 <- (x$FSIQ-100)/15
  x$Practice_order <- x$Timepoint-1
  x$Practice_bin <- ifelse(x$Timepoint==1,'0_baseline',ifelse(x$Timepoint==2,'1_second',
    ifelse(x$Timepoint==3,'2_third',ifelse(x$Timepoint>=4,'3_fourth_or_later',NA))))
  invalid_timepoint<-!is.finite(x$Timepoint)|x$Timepoint!=round(x$Timepoint)|x$Timepoint<1
  x$Practice_bin[invalid_timepoint]<-NA_character_
  score <- x$SWM_Strategy/ifelse(x$Platform=='RS6',4,1)
  reason <- rep('ok',nrow(x))
  reason[is.finite(x$SWM_Strategy)&abs(x$SWM_Strategy-round(x$SWM_Strategy))>1e-8]<-'noninteger_raw'
  reason[is.finite(score)&(score<2-1e-8|score>14+1e-8)]<-'outside_harmonized_2_14'
  reason[is.finite(x$SWM_Strategy)&x$SWM_Strategy<=0]<-'nonpositive_placeholder'
  reason[!is.finite(x$SWM_Strategy)]<-'missing_raw'
  x$SWM_Strategy_QC_reason<-reason
  x$Gender[!nzchar(as.character(x$Gender))]<-NA
  x <- swm_encode_all(x)
  x$eligible <- x$eligible & !invalid_timepoint & x$Time_since_baseline_1y>=0
  x$eligible[is.na(x$eligible)] <- FALSE
  swm_validate_baselines(x) # reject duplicate original baseline rows before subsetting
  x
}

swm_folds <- function(hc,cfg) {
  suppressPackageStartupMessages(library(dplyr))
  result <- make_grouped_folds(hc,K=5L,seed=cfg$fold_seed,n_candidates=cfg$fold_candidates)
  hc$CV_fold <- result$folds
  swm_check_folds(droplevels(hc))
  if(!identical(sort(unique(hc$CV_fold)),1:5))stop('All five folds are required')
  list(data=hc,audit=result$audit,balance_score=result$balance_score,selected_candidate=result$selected_candidate)
}

swm_assert_key <- function(saved,current) {
  if(!identical(saved,current))stop('Cache fingerprint differs; use a new output_dir with mode="fresh". Existing results are preserved.')
  invisible(TRUE)
}

swm_assert_training <- function(fit,train,test=NULL) {
  # brms keeps only model variables; original data.frame row names preserve the
  # input row identity even when row_id_global is omitted by brms.
  if(!is.data.frame(fit$data)||nrow(fit$data)!=nrow(train)||!identical(rownames(fit$data),rownames(train)))stop('Saved fit training rows differ')
  if('row_id_global'%in%names(fit$data)&&!identical(as.character(fit$data$row_id_global),as.character(train$row_id_global)))stop('Saved fit row IDs differ')
  required<-c('Subject_ID','Platform','Gender','Practice_bin','index','grid',
    'Age_baseline_10y','Time_since_baseline_1y','FSIQ_15')
  for(n in intersect(names(train),names(fit$data))) {
    a<-fit$data[[n]];b<-train[[n]]
    if(is.factor(a)||is.factor(b)){a<-as.character(a);b<-as.character(b)}
    if(!isTRUE(all.equal(a,b,check.attributes=FALSE)))stop('Saved fit training data differ: ',n)
  }
  if(!all(required%in%names(fit$data)))stop('Saved fit lacks required training covariates')
  if(!is.null(test)&&length(intersect(as.character(fit$data$Subject_ID),as.character(test$Subject_ID))))stop('OOF participant leakage')
  invisible(TRUE)
}

swm_with_lock <- function(out,code) {
  lock<-file.path(out,'.run_lock')
  if(!dir.create(lock,showWarnings=FALSE))stop('Output is locked by another process: ',lock,
    '. If a process was killed, verify it has stopped before manually removing the lock directory.')
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(c(paste('PID',Sys.getpid()),paste('Host',Sys.info()[['nodename']]),as.character(Sys.time())),file.path(lock,'owner.txt'))
  force(code)
}

swm_config <- function(root) {
  e<-new.env(parent=baseenv());sys.source(file.path(root,'config.R'),e);cfg<-e$NORMATIVE_CONFIG
  cfg$cores<-e$cantab_worker_cores(cfg$cores)
  numeric_cfg<-cfg
  if(any(!vapply(numeric_cfg,function(x)length(x)==1&&is.numeric(x)&&is.finite(x)&&x>0,logical(1))))stop('Invalid configuration')
  ints<-setdiff(names(numeric_cfg),c('adapt_delta','init_radius'))
  if(any(vapply(cfg[ints],function(x)x!=round(x),logical(1)))||cfg$chains!=4||cfg$cores<1||cfg$cores>4||cfg$iter<=cfg$warmup||cfg$ndraws<2000||cfg$adapt_delta>=1)stop('Invalid production configuration')
  cfg
}

swm_context <- function(root,out) {
  suppressPackageStartupMessages(library(brms))
  cfg<-swm_config(root)
  input<-Sys.getenv('NM_INPUT_FILE','')
  if(!nzchar(input))stop('NM_INPUT_FILE is required; use run_outcome() to supply the input CSV.')
  input<-normalizePath(input,mustWork=TRUE)
  raw<-read.csv(input,check.names=FALSE,fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE,
    colClasses=c(Subject_ID='character'),na.strings=c('','NA','NaN'))
  all<-swm_adapt_raw(raw)
  hc<-droplevels(all[all$eligible & all$Diagnosis%in%'HC',])
  folds<-swm_folds(hc,cfg);hc<-folds$data
  # Any HC row excluded for outcome/covariates keeps its participant's fold for audit.
  subj<-!duplicated(as.character(hc$Subject_ID))
  all$CV_fold<-hc$CV_fold[which(subj)][match(as.character(all$Subject_ID),as.character(hc$Subject_ID[subj]))]
  code<-sort(c(file.path(root,'run_model.R'),file.path(root,'config.R'),file.path(root,'cantab_design_audit.R'),file.path(root,'run.R'),list.files(file.path(root,'R'),full.names=TRUE),list.files(file.path(root,'stan'),full.names=TRUE)))
  hashes<-vapply(code,function(p)digest::digest(file=p,algo='sha256'),character(1));names(hashes)<-substring(code,nchar(root)+2L)
  pkgs<-c('brms','cmdstanr','posterior','loo','dplyr','digest','mgcv','rstan','StanHeaders')
  fingerprint<-list(schema='cantab-strategy-normative',input_sha256=digest::digest(file=input,algo='sha256'),
    code_sha256=hashes,config=cfg,R=R.version.string,RNG=RNGkind(),platform=R.version$platform,
    packages=vapply(pkgs,function(p)as.character(utils::packageVersion(p)),character(1)),
    cmdstan=as.character(cmdstanr::cmdstan_version()),
    folds=hc[,c('row_id_global','Subject_ID','CV_fold')],compiler=Sys.getenv(c('CXX','CXXFLAGS','STAN_THREADS','MAKEFLAGS')))
  manifest<-file.path(out,'run_manifest.rds')
  if(file.exists(manifest))swm_assert_key(readRDS(manifest)$fingerprint,fingerprint)
  else swm_atomic_rds(list(created=Sys.time(),input_path=input,fingerprint=fingerprint),manifest)
  stanvars<-brms::stanvar(scode=paste(readLines(file.path(root,'stan/swm_discrete_functions.stan')),collapse='\n'),block='functions')
  list(root=root,out=out,input=input,raw=raw,all=all,hc=hc,cfg=cfg,fingerprint=fingerprint,fold_audit=folds$audit,stanvars=stanvars)
}

# Restore raw-score cell text only in the CSV copy; numeric model data stay intact.
swm_write_score_table <- function(result,ctx) {
 cells<-read.csv(ctx$input,check.names=FALSE,colClasses='character',stringsAsFactors=FALSE,
   fileEncoding='UTF-8-BOM',na.strings=c('','NA','NaN'))
 if(nrow(cells)!=nrow(result) || !identical(cells$Subject_ID,as.character(result$Subject_ID)) ||
    !identical(suppressWarnings(as.numeric(cells$SWM_Strategy)),as.numeric(result$SWM_Strategy)))
  stop('Input rows or raw score values changed before export')
 export<-result;export$SWM_Strategy<-cells$SWM_Strategy
 write.csv(export,file.path(ctx$out,'SWM_Strategy_scores.csv'),row.names=FALSE,na='')
 invisible(result)
}

swm_export <- function(ctx,pred) {
  x<-ctx$all;original<-ctx$raw;result<-original;ix<-match(pred$row_ids,x$row_id_global)
  if(anyNA(ix)||anyDuplicated(ix)||!all(x$eligible[ix])||sum(x$eligible)!=length(ix))stop('Export prediction rows differ')
  # Canonical derived fields are appended only; raw fields and raw order stay intact.
  for(n in c('row_id_global','Age_baseline_years','Time_since_baseline_years','Time_since_baseline_1y',
             'Age_baseline_10y','FSIQ_15','Practice_order','Practice_bin','SWM_Strategy_QC_reason'))
    if(!n%in%names(result))result[[n]]<-x[[n]]
  result$normative_model<-'normative';result$harmonized_score<-x$score;result$score_grid<-x$grid;result$cv_fold<-x$CV_fold
  result$normative_score_source<-ifelse(!x$eligible,'excluded',ifelse(x$Diagnosis%in%'HC','participant_OOF','full_HC_new_person'))
  result$score_status<-ifelse(x$eligible,'ok',ifelse(x$SWM_Strategy_QC_reason=='ok','missing_or_invalid_covariate',x$SWM_Strategy_QC_reason))
  if(!all(c('z_mid','log_mid_p_lower','log_mid_p_upper')%in%names(pred$rows)))
    stop('State predictions lack stable log-tail scores; regenerate predictions')
  result$normative_z<-NA_real_;result$normative_z[ix]<-pred$rows$z_mid
  for(n in c('log_mid_p_lower','log_mid_p_upper')) {
    result[[n]]<-NA_real_
    result[[n]][ix]<-pred$rows[[n]]
  }
  required_precision<-c('posterior_precision_ok','numerical_ok','posterior_mcse_z')
  if(!all(required_precision%in%names(pred$rows))||anyNA(pred$rows[required_precision])||
     any(!pred$rows$posterior_precision_ok|!pred$rows$numerical_ok|!is.finite(pred$rows$z_mid)))
    stop('Primary state precision unresolved or missing; regenerate predictions')
  result$precision_status<-ifelse(x$eligible,'passed','not_scored')
  for(n in required_precision) {
    result[[n]]<-NA
    result[[n]][ix]<-pred$rows[[n]]
  }
  edge<-!is.finite(pred$rows$z_mid)
  result$normative_z[ix[edge]]<-NA_real_;result$score_status[ix[edge]]<-'numeric_tail_limit'
  if(nrow(result)!=nrow(original))stop('Raw row preservation failed')
  for(n in names(original))if(!identical(result[[n]],original[[n]]))stop('Raw column changed: ',n)
  swm_write_score_table(result,ctx)
  swm_write_csv(as.data.frame(table(result$Diagnosis,result$score_status,useNA='ifany')),ctx$out,'export_status_counts.csv')
  result
}
