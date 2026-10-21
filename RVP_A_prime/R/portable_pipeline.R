# R-only orchestration and fail-closed assembly for the normative analysis.
aprime_hash <- function(path)digest::digest(file=path,algo='sha256')
aprime_json <- function(x,path)jsonlite::write_json(x,path,pretty=TRUE,auto_unbox=TRUE,digits=NA,null='null')
aprime_identity_token <- function(x)digest::digest(jsonlite::toJSON(x,auto_unbox=TRUE,digits=NA,null='null'),algo='sha256',serialize=FALSE)
aprime_inventory <- function(root) {
  paths<-c(file.path(root,'run.R'),file.path(root,'config.json'),file.path(root,'cantab_design_audit.R'),
    list.files(file.path(root,'R'),pattern='[.]R$',full.names=TRUE),file.path(dirname(root),c('R/input_data.R','R/results.R')))
  paths<-sort(normalizePath(paths,mustWork=TRUE,winslash='/'))
  as.list(setNames(vapply(paths,aprime_hash,character(1)),paths))
}
aprime_reference_folds <- function(reference_dir) {
  if(is.null(reference_dir))return(character())
  paths<-c(file.path(reference_dir,'primary_subject_folds.csv'),
    file.path(reference_dir,'reference','primary_subject_folds.csv'))
  paths<-paths[file.exists(paths)]
  if(length(paths)>1L)stop('Ambiguous reference: more than one primary_subject_folds.csv.')
  normalizePath(paths,mustWork=TRUE,winslash='/')
}
aprime_reference_inventory <- function(reference_dir) {
  if(is.null(reference_dir))return(list())
  paths<-c(file.path(reference_dir,c('diagnosis_reference.csv','diagnosis_reference.json')),
    aprime_reference_folds(reference_dir))
  paths<-paths[file.exists(paths)]
  paths<-sort(normalizePath(paths,mustWork=TRUE,winslash='/'))
  as.list(setNames(vapply(paths,aprime_hash,character(1)),paths))
}
aprime_status <- function(run_dir,status,...) {
  aprime_json(c(list(status=status,utc=format(Sys.time(),tz='UTC',usetz=TRUE)),list(...)),
    file.path(run_dir,'run_status.json'))
}
aprime_assert_identity <- function(manifest,root,input_file,reference_dir) {
  if(!identical(manifest$raw_path,input_file)||!identical(manifest$raw_sha256,aprime_hash(input_file)))
    stop('Raw input changed; use a fresh output directory')
  if(!identical(manifest$code,aprime_inventory(root)))stop('Analysis code changed; use a fresh output directory')
  if(!identical(manifest$references,aprime_reference_inventory(reference_dir)))stop('References changed; use a fresh output directory')
}
aprime_seal <- function(run_dir,stage,files) {
  if(!length(files)||any(!file.exists(files)))stop('Missing artifacts for stage ',stage)
  files<-sort(normalizePath(files,winslash='/',mustWork=TRUE))
  relative<-substring(files,nchar(normalizePath(run_dir,winslash='/',mustWork=TRUE))+2L)
  if(any(grepl('^[.][.]',relative)))stop('Artifacts must be within output directory')
  path<-file.path(run_dir,'artifact_seals.rds')
  seals<-if(file.exists(path))readRDS(path) else list()
  seals[[stage]]<-list(identity=trimws(readLines(file.path(run_dir,'identity.sha256'))),
    files=as.list(setNames(vapply(files,aprime_hash,character(1)),relative)))
  saveRDS(seals,path,version=2)
}
aprime_verify_seals <- function(run_dir,stages=NULL) {
  path<-file.path(run_dir,'artifact_seals.rds')
  if(!file.exists(path)) {if(length(stages))stop('Missing stage seals');return(invisible(TRUE))}
  seals<-readRDS(path)
  if(is.null(stages))stages<-names(seals)
  for(stage in stages) {
    seal<-seals[[stage]]
    if(is.null(seal)||!length(seal$files)||!identical(seal$identity,trimws(readLines(file.path(run_dir,'identity.sha256')))))
      stop('Missing or incompatible stage seal: ',stage)
    for(name in names(seal$files))if(!file.exists(file.path(run_dir,name))||
      !identical(seal$files[[name]],aprime_hash(file.path(run_dir,name))))stop('Stage artifact changed: ',name)
  }
  invisible(TRUE)
}
aprime_run_stage <- function(stage,root,run_dir) {
  executable<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
  log<-file.path(run_dir,paste0(stage,'.log'))
  message('A-prime stage: ',stage,'; log: ',log)
  status<-system2(executable,c('--vanilla',shQuote(file.path(root,'R/driver.R')),stage,shQuote(run_dir)),stdout=log,stderr=log)
  if(!identical(status,0L)) {
    lines<-readLines(log,warn=FALSE)
    stop('R stage ',stage,' failed. See ',log,'\n',paste(tail(lines,25L),collapse='\n'),call.=FALSE)
  }
  invisible(log)
}
aprime_assembled_rows <- function(run_dir) {
  raw<-aprime_read_cells(file.path(run_dir,'input/aligned_raw.csv'))
  dat<-aprime_read_cells(file.path(run_dir,'input/prepared_data.csv'))
  scores<-aprime_read_cells(file.path(run_dir,'results/normative_score_precision.csv'))
  truth<-function(x)tolower(x)%in%c('true','1')
  number<-function(x)suppressWarnings(as.numeric(x))
  if(nrow(raw)!=nrow(dat)||anyDuplicated(dat$row_id_global)||
     !identical(number(dat$row_id_global),as.numeric(seq_len(nrow(dat)))))stop('Prepared/raw row identity mismatch')
  for(k in c('Subject_ID','Platform','Diagnosis'))if(!identical(raw[[k]],dat[[k]]))stop('Prepared row identity changed: ',k)
  if(!identical(number(raw$Timepoint),number(dat$Timepoint)))stop('Prepared Timepoint changed')
  ids<-dat$row_id_global[truth(dat$model_eligible)]
  if(anyDuplicated(scores$row_id_global)||!setequal(scores$row_id_global,ids)||
     any(!grepl('^normative_',scores$role))||any(scores$precision_status!='ok')||
     any(!is.finite(number(scores$z))))stop('Score key coverage, role or precision mismatch')
  index<-match(scores$row_id_global,dat$row_id_global)
  if(any(grepl('OOF',scores$role)!=(dat$Diagnosis[index]=='HC')))stop('Patient/HC model source mismatch')
  out<-raw
  for(name in setdiff(names(dat),names(raw)))out[[name]]<-dat[[name]]
  out$row_id_global<-dat$row_id_global
  index<-match(dat$row_id_global,scores$row_id_global)
  take<-function(column){v<-scores[[column]][index];v[is.na(v)]<-'';v}
  out$normative_z<-take('z')
  out$precision_status<-ifelse(is.na(index),'not_scored','ok')
  out$quadrature_error<-take('quadrature_z_difference')
  out$posterior_mcse<-take('posterior_mcse_z')
  outside<-rep(FALSE,nrow(dat))
  for(platform in unique(dat$Platform)) {
    ix<-which(dat$Platform==platform);hc<-which(dat$Platform==platform&dat$Diagnosis=='HC'&truth(dat$model_eligible))
    if(!length(hc)) {
      if(any(truth(dat$model_eligible[ix])))stop('No eligible HC covariate support for platform ',platform)
      outside[ix]<-TRUE
      next
    }
    for(nm in c('Age_baseline_years','Time_since_baseline_years','FSIQ')) {
      r<-range(number(dat[[nm]][hc]));x<-number(dat[[nm]][ix])
      outside[ix]<-outside[ix]|!is.finite(x)|x<r[1]|x>r[2]
    }
  }
  out$outside_hc_support<-ifelse(outside,'TRUE','FALSE')
  out$score_source<-ifelse(dat$Diagnosis=='HC','HC_SUBJECT_OOF','NONHC_FULL_HC_MODEL')
  out$score_status<-ifelse(is.na(index),'not_scored','ok')
  out
}
aprime_assemble <- function(run_dir) {
  rows<-aprime_assembled_rows(run_dir);out<-file.path(run_dir,'results')
  write.csv(rows,file.path(out,'RVP_A_prime_scores.csv'),row.names=FALSE,na='')
  counts<-do.call(rbind,lapply(c('HC','BD','MDD'),function(dx){
    ix<-rows$Diagnosis==dx&is.finite(suppressWarnings(as.numeric(rows$normative_z)))
    data.frame(column='normative_z',Diagnosis=dx,n_rows=sum(ix),n_subjects=length(unique(rows$Subject_ID[ix])))
  }))
  write.csv(counts,file.path(out,'score_availability.csv'),row.names=FALSE)
  dictionary<-c(normative_z='Higher is better relative to the healthy-control predictive distribution given the modeled covariates',
    outside_hc_support='Outside at least one platform-specific HC covariate range, or no eligible HC support for that platform; not a convex hull')
  dictionary<-dictionary[names(dictionary)%in%names(rows)]
  write.csv(data.frame(column=names(dictionary),meaning=unname(dictionary)),file.path(out,'score_dictionary.csv'),row.names=FALSE)
  invisible(rows)
}
aprime_verify_artifacts <- function(run_dir,manifest,root,input_file,reference_dir) {
  aprime_assert_identity(manifest,root,input_file,reference_dir)
  identity<-aprime_identity_token(manifest)
  if(!identical(trimws(readLines(file.path(run_dir,'identity.sha256'))),identity))stop('Run identity corrupted')
  cfg<-jsonlite::fromJSON(file.path(run_dir,'config.json'))
  if(!identical(cfg,manifest$configuration))stop('Saved configuration changed')
  aprime_verify_seals(run_dir,c('prepare','score','prior','assemble','calibrate'))
  # Reconstruct cohort selection and raw alignment from the supplied source;
  # saved prepared data alone is not evidence of source preservation.
  temp<-tempfile('aprime_source_verify_');dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  raw_input<-cantab_prepare_input(input_file,temp,reference_dir=reference_dir,outcome='RVP_A_prime')
  aprime_prepare_raw(raw_input,file.path(temp,'input'))
  if(!identical(aprime_read_cells(file.path(temp,'input/aligned_raw.csv')),
    aprime_read_cells(file.path(run_dir,'input/aligned_raw.csv'))))stop('Aligned source rows changed')
  if(!identical(aprime_read_cells(file.path(run_dir,'results/RVP_A_prime_scores.csv')),aprime_assembled_rows(run_dir)))
    stop('Assembled CSV differs from original raw cells or accepted scores')
  model_rows<-read.csv(file.path(run_dir,'models/model_manifest.csv'),stringsAsFactors=FALSE)
  if(nrow(model_rows)!=6L||anyDuplicated(model_rows$path))stop('Six full/OOF fits are required')
  hashes<-as.list(setNames(vapply(basename(model_rows$path),function(p)aprime_hash(file.path(run_dir,'models',p)),character(1)),basename(model_rows$path)))
  if(!identical(unname(unlist(hashes)),model_rows$file_sha256))stop('Model manifest hashes stale')
  s<-read.csv(file.path(run_dir,'results/normative_score_precision.csv'),stringsAsFactors=FALSE)
  if(anyNA(match(s$model_file,names(hashes)))||any(s$model_sha256!=unlist(hashes)[s$model_file]))stop('Score belongs to another model')
  if(any(s$precision_status!='ok')||any(!is.finite(s$z))||any(!is.finite(s$quadrature_z_difference))||
     any(!is.finite(s$posterior_mcse_z))||any(s$quadrature_z_difference>cfg$scoring$tolerance)||
     any(s$posterior_mcse_z>cfg$scoring$mcse_limit))stop('Normative score precision exceeded')
  stamp<-jsonlite::fromJSON(file.path(run_dir,'models/verification_stamp.json'))
  if(!identical(stamp$identity,identity)||!identical(stamp$model_hashes,hashes))stop('Actual fit verification stale')
  d<-read.csv(file.path(run_dir,'input/prepared_data.csv'),stringsAsFactors=FALSE)
  n_hc<-sum(d$Diagnosis=='HC'&d$model_eligible)
  prior<-read.csv(file.path(run_dir,'results/prior_direct/numerical_diagnostics.csv'),stringsAsFactors=FALSE)
  n<-setNames(prior$count,prior$metric)
  if(anyNA(n[c('independent_parameter_draws','observations','approximated_beta_draws')])||
     n['independent_parameter_draws']!=cfg$prior$draws||n['observations']!=cfg$prior$draws*n_hc||n['approximated_beta_draws']!=0)
    stop('Direct prior simulation incomplete or approximate')
  cal<-read.csv(file.path(run_dir,'results/OOF_calibration.csv'),stringsAsFactors=FALSE)
  overall<-cal[cal$stratum=='overall',,drop=FALSE]
  if(nrow(overall)!=1L||overall$metric!='normative'||overall$n_rows!=n_hc||!is.finite(overall$sd))stop('OOF calibration invalid')
  data.frame(check=c('identity_and_source_preservation','all_score_precision_and_model_identity',
    'actual_fits_and_subject_holdout','prior_and_oof_calibration'),status='PASS',
    detail=c('Input, settings, source and sealed artifacts match',
      'All eligible normative scores accepted against unchanged precision thresholds',
      'Six actual fits pass diagnostics and exact raw training/holdout checks',
      'Direct prior and HC OOF calibration complete'))
}
