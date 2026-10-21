# Grouped-validation checkpoint and retry policy.
rti_sampling_plan <- function() {
  # Same 4,000 retained draws per chain in both attempts: brms::kfold_predict
  # allocates its prediction matrix from the first fold's retained draw count.
  data.frame(iter=c(6000L,10000L),warmup=c(2000L,2000L),thin=c(1L,2L))
}

rti_atomic_save <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile('writing_',tmpdir=dirname(path))
  on.exit(if(file.exists(tmp))unlink(tmp),add=TRUE)
  saveRDS(x,tmp)
  if(!file.rename(tmp,path))stop('Could not publish checkpoint: ',path)
  invisible(path)
}

rti_attempts <- function(path,sample,diagnose) {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  for(i in 1:2) {
    file<-file.path(path,paste0('attempt_',i,'.rds'))
    if(file.exists(file)) {
      message('Reusing saved CV attempt: ',file)
      value<-readRDS(file)
    } else {
      value<-sample(i)
      # Save the complete fits before any diagnostic can throw an error.
      rti_atomic_save(value,file)
    }
    error<-tryCatch({diagnose(value,i);NULL},error=identity)
    if(is.null(error))return(value)
    reason<-conditionMessage(error)
    writeLines(reason,file.path(path,paste0('attempt_',i,'_failure.txt')))
    if(!startsWith(reason,'CV convergence gate failed'))stop(error)
    message('CV attempt ',i,' failed the unchanged convergence criteria: ',reason)
    rm(value);gc()
  }
  stop('Both CV attempts failed convergence. Saved fits and diagnostics: ',path)
}

rti_run_grouped_kfold <- function(fit,folds,save_fits,seed) {
  plan<-rti_sampling_plan()
  stopifnot(cv_iter==plan$iter[1],cv_warmup==plan$warmup[1])
  # Both model data and actual parent posterior draws identify a CV calculation.
  key<-digest::digest(list(data=fit$data,folds=folds,seed=seed,plan=plan,
    joint=cv_joint_likelihood,chains=cv_chains,adapt_delta=adapt_delta_fit,max_treedepth=max_treedepth_fit,
    parent_draws=digest::digest(posterior::as_draws_array(fit),algo='sha256')),
    algo='sha256')
  path<-file.path(cv_checkpoint_dir,'cv_attempts',key)
  sample<-function(i) {
    settings<-plan[i,]
    message('CV seed=',seed,' attempt=',i,' iter=',settings$iter,
            ' warmup=',settings$warmup,' thin=',settings$thin)
    # All folds use the same retained draw count, including on a longer retry.
    # Use the same seed: retries increase effort, never search across seeds.
    set.seed(seed)
    cantab_kfold(fit,K=max(folds),folds=folds,group='Subject_ID',
      joint=cv_joint_likelihood,save_fits=TRUE,compare=FALSE,
      chains=cv_chains,iter=settings$iter,warmup=settings$warmup,thin=settings$thin,
      cores=cv_cores,seed=seed,refresh=100,file=NULL,
      control=list(adapt_delta=adapt_delta_fit,max_treedepth=max_treedepth_fit))
  }
  diagnose<-function(kf,i) {
    rti_validate_kfold(kf,fit$data,paste0(movement_outcome,'_rti_',
      fit$family$family,'_attempt',i),seed,diagnose_brms_fit,max_treedepth_fit)
  }
  kf<-rti_attempts(path,sample,diagnose)
  # Disk checkpoint always keeps fits, even when a comparison only needs ELPD.
  if(!save_fits)kf$fits<-NULL
  kf
}
