resolve_project_path <- function(path,root) {
 if(length(path)!=1L||is.na(path)||!nzchar(path))stop('Path must be a nonempty string')
 path<-path.expand(path)
 absolute<-grepl('^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)',path)
 if(!absolute)path<-file.path(root,path)
 path
}
atomic_text <- function(lines,path) {
 tmp<-paste0(path,'.tmp-',Sys.getpid());writeLines(lines,tmp)
 if(!file.rename(tmp,path))stop('Cannot finalize marker: ',path)
}
acquire_run_lock <- function(ctx) {
 path<-file.path(ctx$run,'.run.lock')
 if(!dir.create(path,showWarnings=FALSE))
  stop('This run is locked: ',path,'. If a previous process crashed, confirm it has stopped before removing that lock directory.')
 writeLines(c(paste('pid',Sys.getpid()),paste('host',Sys.info()[['nodename']]),
  paste('started',Sys.time())),file.path(path,'owner.txt'))
 function() unlink(path,recursive=TRUE)
}
invalidate_stage <- function(ctx,stage) {
 stopifnot(stage%in%c('fit','score'))
 files<-file.path(ctx$run,c('diagnostic_summaries.zip','aggregate/NORMATIVE_COMPLETE.txt'))
 for(f in files)if(file.exists(f)&&!file.remove(f))stop('Cannot invalidate previous completion marker: ',f)
 invisible(NULL)
}
resolve_context <- function(root,config=CONFIG) {
 root<-normalizePath(root,mustWork=TRUE)
 input<-Sys.getenv('SWM_INPUT','')
 if(!nzchar(input))stop('SWM_INPUT is required; use run_outcome() to supply the input CSV.')
 input<-resolve_project_path(input,root)
 input<-normalizePath(input,mustWork=TRUE)
 code<-c(list.files(file.path(root,'R'),pattern='[.]R$',full.names=TRUE),list.files(file.path(root,'vendor'),pattern='[.]R$',full.names=TRUE),file.path(root,c('run_model.R','config.R','cantab_design_audit.R','run.R')))
 hashes<-vapply(code,function(f)digest::digest(file=f,algo='sha256'),character(1))
 names(hashes)<-substring(code,nchar(root)+2)
 pkgs<-c('brms','posterior','dplyr','mgcv','cmdstanr','digest','loo')
 versions<-vapply(pkgs,function(p)as.character(packageVersion(p)),character(1))
 csvhash<-digest::digest(file=input,algo='sha256')
 signature<-digest::digest(list(config=config,code=hashes,input=csvhash,versions=versions,R=as.character(getRversion()),cmdstan=as.character(cmdstanr::cmdstan_version())),algo='sha256')
 run<-Sys.getenv('SWM_OUTPUT_RUN','')
 if(!nzchar(run))stop('SWM_OUTPUT_RUN is required; use run_outcome() to supply the output directory.')
 run<-cantab_output_directory(resolve_project_path(run,root),dirname(root))
 manifest<-file.path(run,'run_manifest.rds')
 if(file.exists(manifest) && !identical(readRDS(manifest)$signature,signature))stop('Run signature differs; start a fresh run. Existing results are preserved.')
 for(p in c('','models','private','logs','aggregate','stan_cache'))dir.create(file.path(run,p),recursive=TRUE,showWarnings=FALSE)
 if(!file.exists(manifest))atomic_rds(list(signature=signature,input_sha256=csvhash,code=hashes,config=config,created=Sys.time()),manifest)
 options(cmdstanr_write_stan_file_dir=file.path(run,'stan_cache'))
 list(root=root,input=input,run=run,config=config,signature=signature,code=hashes,versions=versions,input_sha256=csvhash)
}
preflight <- function(ctx) {
 message('Preflight: raw CSV, units, folds, ',length(ctx$config$models)*(ctx$config$folds+1L),' model designs; no sampling.')
 if(!file.exists(file.path(cmdstanr::cmdstan_path(),'bin',if(.Platform$OS.type=='windows')'stanc.exe'else'stanc')))stop('CmdStan not installed; complete the shared environment setup')
 raw<-read.csv(ctx$input,check.names=FALSE,stringsAsFactors=FALSE,
   colClasses=c(Subject_ID='character'),fileEncoding='UTF-8-BOM',na.strings=c('','NA','NaN'))
 p<-prepare_raw_data(raw,ctx$config)
 fo<-make_grouped_folds(p$hc,ctx$config$folds,ctx$config$fold_seed,TRUE,ctx$config$fold_candidates)
 p$hc$CV_fold<-fo$folds
 fmap<-unique(data.frame(Subject_ID=as.character(p$hc$Subject_ID),CV_fold=fo$folds))
 stopifnot(!anyDuplicated(fmap$Subject_ID))
 p$all$CV_fold<-fmap$CV_fold[match(as.character(p$all$Subject_ID),fmap$Subject_ID)]
 p$valid$CV_fold<-fmap$CV_fold[match(as.character(p$valid$Subject_ID),fmap$Subject_ID)]
 ranks<-list();spec<-character()
 for(id in ctx$config$models){
  d<-p$hc
  for(k in 0:ctx$config$folds){
   tr<-if(k==0)d else d[d$CV_fold!=k,,drop=FALSE]
   va<-if(k==0)d[FALSE,,drop=FALSE]else d[d$CV_fold==k,,drop=FALSE]
   if(length(intersect(as.character(tr$Subject_ID),as.character(va$Subject_ID))))stop('Subject leakage')
   if(k>0&&!nrow(va))stop('Empty validation fold')
   def<-fixed_definition(id,tr);r<-check_design(def,tr);r$model<-id;r$fold<-k
   ranks[[length(ranks)+1]]<-r
   if(k==0){
    spec<-c(spec,paste('Model:',id),capture.output(print(def$formula)),capture.output(print(def$priors)))
    path<-file.path(ctx$run,'stan_cache',paste0(id,'.stan'))
    writeLines(brms::make_stancode(def$formula,data=tr,family=def$family,prior=def$priors,stanvars=def$stanvars),path)
    cmdstanr::cmdstan_model(path,compile=FALSE)$check_syntax()
   }
  }
 }
 write.csv(dplyr::bind_rows(ranks),file.path(ctx$run,'aggregate/design_ranks.csv'),row.names=FALSE)
 write.csv(fo$audit,file.path(ctx$run,'aggregate/fold_balance.csv'),row.names=FALSE)
 write.csv(fmap,file.path(ctx$run,'private/HC_fold_assignments.csv'),row.names=FALSE)
 write.csv(data.frame(file=names(ctx$code),sha256=unname(ctx$code)),file.path(ctx$run,'aggregate/code_manifest.csv'),row.names=FALSE)
 write.csv(data.frame(package=names(ctx$versions),version=ctx$versions),file.path(ctx$run,'aggregate/package_versions.csv'),row.names=FALSE)
 writeLines(spec,file.path(ctx$run,'aggregate/model_specifications.txt'))
 writeLines(c(paste('Input SHA256:',ctx$input_sha256),paste('Run signature:',ctx$signature),paste('Units:',p$units),capture.output(sessionInfo())),file.path(ctx$run,'aggregate/provenance.txt'))
 write.csv(data.frame(group=names(table(p$all$Diagnosis)),rows=as.integer(table(p$all$Diagnosis))),file.path(ctx$run,'aggregate/input_group_counts.csv'),row.names=FALSE)
 cache_save(p,file.path(ctx$run,'private/prepared.rds'),ctx$signature)
 writeLines('PASS: raw preprocessing, leakage checks, complete population design ranks for full HC and all five folds, and Stan syntax; no MCMC.',file.path(ctx$run,'aggregate/PREFLIGHT_PASS.txt'))
 p
}
normative_job <- function(id,ctx,p=preflight(ctx)) {
 if(!id%in%ctx$config$models)stop('Unknown model: ',id)
 file<-file.path(ctx$run,'private',paste0('job_',id,'.rds'))
 done<-cache_load(file,ctx$signature);if(!is.null(done))return(done)
 hc<-p$hc
 all<-p$valid
 nonhc<-all[all$Diagnosis!=ctx$config$control_label,,drop=FALSE]
 parts<-list();diags<-list();full<-data.frame()
 for(k in 0:ctx$config$folds){
  scope<-paste0(id,if(k==0)'_full'else paste0('_fold',k))
  sf<-file.path(ctx$run,'private',paste0('prediction_',scope,'.rds'))
  z<-cache_load(sf,ctx$signature)
  if(is.null(z)){
   tr<-if(k==0)hc else hc[hc$CV_fold!=k,,drop=FALSE]
   va<-if(k==0)nonhc else hc[hc$CV_fold==k,,drop=FALSE]
   stopifnot(!length(intersect(as.character(tr$Subject_ID),as.character(va$Subject_ID))))
   seed<-ctx$config$seed+1000+k*10
   obj<-fit_fixed(id,tr,scope,seed,ctx)
   pred<-predict_state(obj$fit,va,id,seed+1,ctx$config$prediction_draws,coverage=TRUE,
    diagnostics_file=file.path(ctx$run,'private',paste0(scope,'_state_precision.csv')))
   if(nrow(pred))pred$outer_fold<-if(k==0)NA_integer_ else k
   diag<-obj$diagnostics;diag$model<-id;diag$scope<-scope;diag$fit_signature<-obj$fit_signature
   z<-list(pred=pred,diagnostics=diag)
   cache_save(z,sf,ctx$signature);rm(obj);gc(FALSE)
  }
  if(!diagnostic_pass(z$diagnostics))stop('Cached diagnostics failed')
  diags[[scope]]<-z$diagnostics
  if(k==0){full<-z$pred}else{parts[[scope]]<-z$pred}
 }
 prior<-data.frame()
 if(ctx$config$prior_predictive){
  pf<-file.path(ctx$run,'private',paste0('prior_summary_',id,'.rds'));prior<-cache_load(pf,ctx$signature)
  if(is.null(prior)){
   obj<-fit_fixed(id,hc,paste0(id,'_prior'),ctx$config$seed+1000+500,ctx,TRUE)
   prior<-summarize_prior_predictions(obj$fit,hc,ctx$config$seed+1000+501);prior$model<-id
   cache_save(prior,pf,ctx$signature);rm(obj);gc(FALSE)
  }
 }
 done<-list(oof=dplyr::bind_rows(parts),full=full,diagnostics=dplyr::bind_rows(diags),prior=prior)
 stopifnot(setequal(done$oof$row_id_global,hc$row_id_global),setequal(done$full$row_id_global,nonhc$row_id_global))
 cache_save(done,file,ctx$signature);done
}
assemble_scores <- function(p,jobs,config) {
 if(!identical(config$models,'normative') || !identical(names(jobs),'normative'))
  stop('Expected one normative model.')
 out<-p$raw
 derived<-setdiff(names(p$all),names(p$raw));out[derived]<-p$all[derived]
 out$normative_score_source<-ifelse(p$all$Diagnosis==config$control_label,'HC_subject_held_out','nonHC_full_HC')
 out$HC_outer_fold<-p$all$CV_fold
 j<-jobs[['normative']]
 out$normative_z<-choose_analysis_score(p$all$Diagnosis,align_scores(p$all,j$oof),
   align_scores(p$all,j$full),config$control_label)
 combined<-dplyr::bind_rows(j$oof,j$full)
 for(n in c('posterior_precision_ok','numerical_ok','posterior_mcse_z')) {
  if(!n%in%names(combined))stop('Predictions lack precision diagnostics; regenerate predictions')
  out[[n]]<-combined[[n]][match(p$all$row_id_global,combined$row_id_global)]
 }
 scored<-is.finite(out$normative_z)
 if(anyNA(out$posterior_precision_ok[scored]) || anyNA(out$numerical_ok[scored]) ||
    any(!out$posterior_precision_ok[scored] | !out$numerical_ok[scored]))
  stop('Unresolved normative score precision')
 stopifnot(nrow(out)==nrow(p$raw))
 out
}
load_jobs <- function(ctx){
 jobs<-setNames(lapply(ctx$config$models,function(id)cache_load(file.path(ctx$run,'private',paste0('job_',id,'.rds')),ctx$signature)),ctx$config$models)
 if(any(vapply(jobs,is.null,logical(1))))stop('Normative estimation is incomplete; run fit/all first.')
 jobs
}
run_fits <- function(ctx,workers=1) {
 invalidate_stage(ctx,'fit')
 preflight(ctx)
 if(length(workers)!=1L||!is.finite(workers)||workers!=floor(workers)||workers<1||workers>4)stop('workers must be an integer in 1..4')
 # Each child runs config$cores chains concurrently; model workers share this job budget.
 workers<-min(workers,max(1L,cantab_job_cores() %/% ctx$config$cores))
 launch<-function(id){
  log<-file.path(ctx$run,'logs',paste0(id,'.log'))
  status<-system2(file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe'else'Rscript'),c('--vanilla',shQuote(file.path(ctx$root,'run_model.R')),'worker',id),stdout=log,stderr=log)
  list(id=id,status=status,log=log)
 }
 ans<-if(workers>1&&.Platform$OS.type=='unix')parallel::mclapply(ctx$config$models,launch,mc.cores=workers,mc.preschedule=FALSE)else lapply(ctx$config$models,launch)
 for(x in ans)if(!is.list(x)||x$status!=0)stop('Normative worker failed; inspect model log. Resume with the same command.')
 invisible(load_jobs(ctx))
}
# Restore raw-score cell text only in the CSV copy; numeric model data stay intact.
write_score_table <- function(scores,ctx) {
 cells<-read.csv(ctx$input,check.names=FALSE,colClasses='character',stringsAsFactors=FALSE,
   fileEncoding='UTF-8-BOM',na.strings=c('','NA','NaN'))
 if(nrow(cells)!=nrow(scores) || !identical(cells$Subject_ID,as.character(scores$Subject_ID)) ||
    !identical(suppressWarnings(as.numeric(cells$SWM_Between_errors)),as.numeric(scores$SWM_Between_errors)))
  stop('Input rows or raw score values changed before export')
 export<-scores;export$SWM_Between_errors<-cells$SWM_Between_errors
 write.csv(export,file.path(ctx$run,'private/SWM_scores.csv'),row.names=FALSE,na='')
 invisible(scores)
}

run_score <- function(ctx){
 invalidate_stage(ctx,'score')
 p<-preflight(ctx);jobs<-load_jobs(ctx);scores<-assemble_scores(p,jobs,ctx$config)
 write_score_table(scores,ctx)
 cache_save(scores,file.path(ctx$run,'private/scores.rds'),ctx$signature)
 report_normative(p,jobs,scores,ctx)
 atomic_text(c('NORMATIVE_COMPLETE',ctx$signature),file.path(ctx$run,'aggregate/NORMATIVE_COMPLETE.txt'))
 scores
}
collect_diagnostics <- function(ctx){
 marker<-file.path(ctx$run,'aggregate/NORMATIVE_COMPLETE.txt')
 if(!file.exists(marker)||!ctx$signature%in%readLines(marker,warn=FALSE))stop('Normative scoring is incomplete; run score/all before collect.')
 marker_state<-identical(readLines(marker,warn=FALSE)[1],'NORMATIVE_COMPLETE')
 if(!marker_state)stop('Normative completion marker is invalid')
 files<-list.files(file.path(ctx$run,'aggregate'),full.names=TRUE)
 if(!length(files))stop('No aggregate outputs available')
 rel<-substring(files,nchar(ctx$run)+2)
 dest<-file.path(ctx$run,'diagnostic_summaries.zip')
 old<-getwd();on.exit(setwd(old),add=TRUE);setwd(ctx$run)
 # New temporary ZIP avoids stale members from an older collection.
 tmp<-paste0(dest,'.tmp.zip');if(file.exists(tmp))unlink(tmp)
 if(!requireNamespace('zip',quietly=TRUE))stop('Install the R package zip before collecting outputs.')
 zip::zipr(tmp,rel,include_directories=FALSE,mode='mirror')
 if(!file.rename(tmp,dest))stop('Cannot finalize diagnostic-summary ZIP')
 message('Diagnostic summaries (aggregate files only): ',dest);invisible(dest)
}
