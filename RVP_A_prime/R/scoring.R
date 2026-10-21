# Numerical scoring from accepted models with verified subject holdouts.
atomic_rds <- function(x,path) {
  tmp<-tempfile(paste0(basename(path),'.'),tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE);saveRDS(x,tmp)
  if(!file.rename(tmp,path))stop('Atomic RDS save failed: ',path)
}
save_score_cache <- function(path,identity,result) {
  atomic_rds(list(identity=identity,result=result,result_sha256=digest::digest(result,algo='sha256')),path)
}
read_score_cache <- function(path,identity) {
  saved<-readRDS(path)
  if(!identical(saved$identity,identity))stop('Score cache identity mismatch: ',path)
  if(!identical(saved$result_sha256,digest::digest(saved$result,algo='sha256')))
    stop('Score cache payload hash mismatch: ',path)
  saved$result
}
atomic_score_csv <- function(x,path) {
  tmp<-tempfile(paste0(basename(path),'.'),tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE)
  write.csv(x,tmp,row.names=FALSE,na='')
  if(!file.rename(tmp,path))stop('Atomic score CSV save failed: ',path)
}
score_all <- function(prepared,config,run_dir,identity) {
  out<-file.path(run_dir,'results');cache<-file.path(run_dir,'score_cache')
  dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(cache,showWarnings=FALSE)
  d<-prepared$data;score_parts<-list()
  fold_by_subject<-setNames(prepared$folds$CV_fold,as.character(prepared$folds$Subject_ID))
  cached_score<-function(name,model_path,dat,fun) {
    token<-digest::digest(list(identity,name,digest::digest(file=model_path,algo='sha256'),dat),algo='sha256')
    path<-file.path(cache,paste0(name,'.rds'))
    if(file.exists(path)) {
      ans<-read_score_cache(path,token);message('Reusing accepted score cache: ',name)
    } else {
      ans<-fun()
      if(any(ans$precision_status!='ok')) {
        write.csv(ans,file.path(cache,paste0(name,'_UNRESOLVED.csv')),row.names=FALSE)
        stop('Numerical precision unresolved: ',name,'; inspect UNRESOLVED.csv. No scores exported.')
      }
      ans$model_file<-basename(model_path);ans$model_sha256<-digest::digest(file=model_path,algo='sha256')
      save_score_cache(path,token,ans)
    }
    stopifnot(nrow(ans)>0,all(ans$precision_status=='ok'),all(is.finite(ans$z)))
    ans
  }
  for(k in 0:5) {
    population<-d[d$model_eligible,,drop=FALSE]
    nd<-if(k==0)population[population$Diagnosis!='HC',,drop=FALSE] else
      population[population$Diagnosis=='HC' & !is.na(fold_by_subject[as.character(population$Subject_ID)]) &
        fold_by_subject[as.character(population$Subject_ID)]==k,,drop=FALSE]
    if(!nrow(nd))next
    path<-file.path(run_dir,'models',paste0('primary',if(k==0)'_full'else paste0('_fold',k),'.rds'))
    fit<-read_accepted_fit(path,identity)
    if(k>0 && any(as.character(nd$Subject_ID)%in%as.character(fit$data$Subject_ID)))stop('Held-out subject leakage')
    prefix<-paste0('normative',if(k==0)'_full'else paste0('_OOF_fold',k))
    score_parts[[prefix]]<-cached_score(prefix,path,nd,function()
      do.call(score_cdf_quadrature,c(list(fit=fit,data=nd,role=prefix),config$scoring)))
  }
  atomic_score_csv(do.call(rbind,score_parts),file.path(out,'normative_score_precision.csv'))
  invisible(TRUE)
}
