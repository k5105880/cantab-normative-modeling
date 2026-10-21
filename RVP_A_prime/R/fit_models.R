# All posterior fits are signed envelopes. A cache is never trusted on its name.
model_hash <- function(x) digest::digest(x,algo="sha256",serialize=TRUE)
atomic_model_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(paste0(basename(path),"."),tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(x,tmp,version=3,compress="gzip")
  if(!file.rename(tmp,path)) stop("Cannot atomically persist fit: ",path)
  invisible(path)
}

write_model_envelope <- function(fit,metadata,identity,diagnostics,accepted,path) {
  # Hash the immutable serialized payload that will actually be persisted.
  # Re-serializing a reconstructed brmsfit is not a file-integrity check: R can
  # restore different reference/promise representations without changing draws.
  payload <- serialize(fit,NULL,version=3,xdr=TRUE)
  e <- list(cache_format='serialized-fit-payload-v2',identity=identity,
    fit_payload=payload,fit_sha256=digest::digest(payload,algo='sha256',serialize=FALSE),
    metadata=metadata,metadata_sha256=model_hash(metadata),diagnostics=diagnostics,accepted=accepted)
  atomic_model_rds(e,path)
}

accept_sampler_diagnostics <- function(d) {
  needed <- c("n_divergent","n_max_treedepth","max_rhat","min_ess_bulk","min_ess_tail","min_EBFMI")
  if(nrow(d)!=1L || !all(needed%in%names(d)) ||
     any(!is.finite(unlist(d[needed]))) || d$n_divergent!=0 || d$n_max_treedepth!=0 ||
     d$max_rhat>=1.01 || d$min_ess_bulk<400 || d$min_ess_tail<400 || d$min_EBFMI<=.3)
    stop("Posterior fit failed sampler acceptance")
  invisible(TRUE)
}

accept_draw_shape <- function(iterations,chains,sampling) {
  if(iterations!=sampling$iter-sampling$warmup || chains!=sampling$chains)
    stop("Posterior fit has incomplete chain draws")
  invisible(TRUE)
}

validate_parameter_diagnostics <- function(summary) {
  fields<-c('rhat','ess_bulk','ess_tail')
  if(!nrow(summary) || !all(fields%in%names(summary)) || any(!is.finite(as.matrix(summary[fields]))))
    stop('Nonfinite parameter diagnostics; a stuck or undefined parameter cannot be ignored')
  invisible(TRUE)
}
summarise_monitored_draws <- function(draws) {
  # For a 2x2 correlation Cholesky factor the first row is deterministically (1,0).
  # Verify exact constants before excluding only these two structural entries.
  variables<-posterior::variables(draws)
  structural<-c('L_1[1,1]','L_1[1,2]')
  if(any(structural%in%variables)) {
    stopifnot(all(c(structural,'L_1[2,1]','L_1[2,2]')%in%variables),
      all(draws[,,'L_1[1,1]']==1),all(draws[,,'L_1[1,2]']==0))
    draws<-posterior::subset_draws(draws,variable=setdiff(variables,structural))
  }
  result<-posterior::summarise_draws(draws,'rhat','ess_bulk','ess_tail')
  validate_parameter_diagnostics(result)
  result
}
validate_nuts_records <- function(np,iterations,chains) {
  if(!all(c('Parameter','Chain','Iteration','Value')%in%names(np)))stop('Missing sampler columns')
  expected<-as.vector(outer(seq_len(chains),seq_len(iterations),paste,sep=':'))
  for(parameter in c('divergent__','treedepth__','energy__')) {
    rows<-np[np$Parameter==parameter,,drop=FALSE]
    keys<-paste(rows$Chain,rows$Iteration,sep=':')
    if(nrow(rows)!=iterations*chains || anyDuplicated(keys) || !setequal(keys,expected) || any(!is.finite(rows$Value)))
      stop('Incomplete or invalid sampler records: ',parameter)
  }
  invisible(TRUE)
}
actual_fit_diagnostics <- function(fit,sampling,enforce=TRUE) {
  if(!inherits(fit,"brmsfit")) stop("Cache does not contain a brmsfit")
  dr <- posterior::as_draws_array(fit)
  accept_draw_shape(posterior::niterations(dr),posterior::nchains(dr),sampling)
  s <- summarise_monitored_draws(dr)
  np <- brms::nuts_params(fit)
  validate_nuts_records(np,posterior::niterations(dr),posterior::nchains(dr))
  energy <- split(np$Value[np$Parameter=="energy__"],np$Chain[np$Parameter=="energy__"])
  if(length(energy)!=sampling$chains) stop("Missing chain energy diagnostics")
  ebfmi <- vapply(energy,function(x) mean(diff(x)^2)/var(x),numeric(1))
  d <- data.frame(n_posterior_draws=posterior::ndraws(dr),n_chains=posterior::nchains(dr),
    n_iterations=posterior::niterations(dr),
    n_divergent=sum(np$Parameter=="divergent__" & np$Value==1),
    n_max_treedepth=sum(np$Parameter=="treedepth__" & np$Value>=sampling$max_treedepth),
    max_rhat=max(s$rhat),min_ess_bulk=min(s$ess_bulk),
    min_ess_tail=min(s$ess_tail),min_EBFMI=min(ebfmi))
  if(enforce) accept_sampler_diagnostics(d)
  d
}

read_model_envelope <- function(path,expected_identity=NULL,expected_accepted=TRUE) {
  stopifnot(identical(expected_accepted,TRUE)||identical(expected_accepted,FALSE))
  e <- readRDS(path)
  required <- c("identity","fit_sha256","metadata_sha256","metadata","diagnostics","accepted")
  if(!is.list(e) || !all(required%in%names(e)) || !identical(e$accepted,expected_accepted)) stop("Unsigned or rejected model cache: ",path)
  if(!is.null(expected_identity) && !identical(e$identity,expected_identity)) stop("Model cache identity mismatch: ",path)
  if(!is.null(e$cache_format)) {
    if(!identical(e$cache_format,'serialized-fit-payload-v2'))stop('Unknown model cache format: ',path)
    if(!is.raw(e$fit_payload)||'fit'%in%names(e))stop('Invalid model payload: ',path)
    fit_hash <- digest::digest(e$fit_payload,algo='sha256',serialize=FALSE)
  } else {
    if(!'fit'%in%names(e))stop('Unsigned model cache: ',path)
    fit_hash <- model_hash(e$fit)
  }
  if(!identical(e$fit_sha256,fit_hash) || !identical(e$metadata_sha256,model_hash(e$metadata)))
    stop("Model cache hash mismatch: ",path)
  if(!is.null(e$cache_format)) {
    e$fit <- unserialize(e$fit_payload)
    e$fit_payload <- NULL
  }
  if(expected_accepted) {
    e$actual_diagnostics <- actual_fit_diagnostics(e$fit,e$metadata$sampling)
  } else {
    d<-tryCatch(actual_fit_diagnostics(e$fit,e$metadata$sampling,enforce=FALSE),error=function(error)error)
    if(inherits(d,'error'))e$diagnostic_error<-conditionMessage(d) else e$actual_diagnostics<-d
  }
  e
}

read_accepted_fit <- function(path,expected_identity=NULL) {
  e <- read_model_envelope(path,expected_identity)
  fit <- e$fit
  attr(fit,"model_metadata") <- e$metadata
  attr(fit,"actual_diagnostics") <- e$actual_diagnostics
  fit
}

fit_tasks <- function(prepared,config,run_dir) {
  sampling <- config$sampling
  needed <- c("chains","iter","warmup","cores","adapt_delta","max_treedepth")
  if(!all(needed%in%names(sampling))) stop("Missing sampling configuration")
  if(sampling$chains<2L || sampling$iter<=sampling$warmup || sampling$cores<1L)
    stop("Invalid posterior sampling configuration")
  folds <- prepared$folds
  hc <- prepared$hc
  if(is.list(folds) && !is.data.frame(folds) && !is.null(folds$assignment)) folds <- folds$assignment
  if(is.list(folds) && !is.data.frame(folds) && !is.null(folds$folds)) folds <- folds$folds
  if(!is.data.frame(folds)) folds <- data.frame(Subject_ID=hc$Subject_ID,CV_fold=as.integer(folds))
  if(!all(c("Subject_ID","CV_fold")%in%names(folds))) stop("Fold map requires Subject_ID and CV_fold")
  map <- unique(data.frame(Subject_ID=as.character(folds$Subject_ID),CV_fold=as.integer(folds$CV_fold)))
  if(anyDuplicated(map$Subject_ID) || anyNA(map$CV_fold) || !setequal(unique(map$CV_fold),1:5))
    stop("Invalid grouped five-fold map")
  if(!setequal(map$Subject_ID,as.character(hc$Subject_ID))) stop("Fold subjects do not match complete-case HC")
  tasks <- list()
  {
    kind <- "primary"
    d <- hc
    if(!nrow(d) || any(as.character(d$Diagnosis)!="HC") || anyDuplicated(d$row_id_global))
      stop("Invalid HC training data for ",kind)
    f <- map$CV_fold[match(as.character(d$Subject_ID),map$Subject_ID)]
    if(anyNA(f) || !setequal(unique(f),1:5)) stop("Training data must cover all five folds")
    for(k in 0:5) {
      ix <- if(k==0L) seq_len(nrow(d)) else which(f!=k)
      held <- if(k==0L) character() else sort(unique(as.character(d$Subject_ID[f==k])))
      train <- droplevels(d[ix,,drop=FALSE])
      if(length(intersect(as.character(train$Subject_ID),held))) stop("OOF subject leakage")
      task <- list(kind=kind,fold=k,train=train,
        path=file.path(run_dir,"models",paste0(kind,if(k==0L) "_full" else paste0("_fold",k),".rds")),
        metadata=list(kind=kind,fold=k,train_ids=as.integer(train$row_id_global),
          held_subjects=held,train_sha256=model_hash(train),sampling=sampling,
          seed=as.integer(config$seed+10000L+k*101L)))
      if(!is.null(config$sampling_recovery)) {
        policy<-config$sampling_recovery
        stopifnot(identical(policy$trigger,'divergence_only'),policy$max_retries==1L,
          is.finite(policy$warmup),policy$warmup==floor(policy$warmup),policy$warmup>sampling$warmup,
          is.finite(policy$adapt_delta),policy$adapt_delta>sampling$adapt_delta,policy$adapt_delta<1)
        task$recovery_sampling<-sampling
        task$recovery_sampling$warmup<-as.integer(policy$warmup)
        task$recovery_sampling$iter<-as.integer(policy$warmup+sampling$iter-sampling$warmup)
        task$recovery_sampling$adapt_delta<-policy$adapt_delta
      }
      tasks[[length(tasks)+1L]]<-task
    }
  }
  tasks
}

validate_task_fit <- function(envelope,task) {
  matches<-identical(envelope$metadata,task$metadata)||
    (!is.null(task$recovery_sampling)&&identical(envelope$metadata,recovery_task(task)$metadata))
  if(!matches) stop("Fit cache does not match expected raw training rows, held subjects, or settings")
  fit <- envelope$fit
  if(!identical(fit$family$family,"xbeta")) stop("Unexpected posterior family")
  expected <- build_model(task$train)
  if(!identical(deparse(fit$formula$formula),deparse(expected$formula$formula)) ||
     !identical(deparse(fit$formula$pforms$phi),deparse(expected$formula$pforms$phi)))
    stop("Unexpected fitted predictors")
  for(nm in names(fit$data)) {
    if(!nm%in%names(task$train) || !identical(as.character(fit$data[[nm]]),as.character(task$train[[nm]])))
      stop("Fit training data mismatch in ",nm)
  }
  if(length(intersect(as.character(fit$data$Subject_ID),task$metadata$held_subjects))) stop("OOF subject leakage in fitted object")
  invisible(TRUE)
}

write_model_manifest <- function(tasks,envelopes,run_dir) {
  manifest <- do.call(rbind,Map(function(t,e) data.frame(kind=t$kind,fold=t$fold,path=normalizePath(t$path),
    file_sha256=digest::digest(file=t$path,algo="sha256"),fit_sha256=e$fit_sha256,
    n_train_rows=length(e$metadata$train_ids),n_train_subjects=length(unique(t$train$Subject_ID)),
    n_held_subjects=length(e$metadata$held_subjects),seed=e$metadata$seed,
    sampling_attempt=if(identical(e$metadata, t$metadata))'initial' else 'retry1',
    chains=e$metadata$sampling$chains,iter=e$metadata$sampling$iter,warmup=e$metadata$sampling$warmup,
    adapt_delta=e$metadata$sampling$adapt_delta,max_treedepth=e$metadata$sampling$max_treedepth,
    initial_rejection_sha256=if(file.exists(rejection_path(t)))digest::digest(file=rejection_path(t),algo='sha256') else '',
    stringsAsFactors=FALSE),tasks,envelopes))
  diagnostics <- do.call(rbind,Map(function(t,e)cbind(kind=t$kind,fold=t$fold,seed=e$metadata$seed,e$actual_diagnostics),tasks,envelopes))
  dir.create(file.path(run_dir,"models"),recursive=TRUE,showWarnings=FALSE)
  for(item in c("manifest","diagnostics")) {
    path <- file.path(run_dir,"models",paste0("model_",item,".csv"))
    tmp <- tempfile(tmpdir=dirname(path))
    write.csv(get(item),tmp,row.names=FALSE,na="")
    if(!file.rename(tmp,path)) stop("Cannot persist model audit")
  }
  manifest
}

# The same arguments are used to generate code serially and to sample each task.
fit_task_arguments <- function(task) {
  model <- build_model(task$train)
  model_design_audit(model,task$train)
  s <- task$metadata$sampling
  list(formula=model$formula,data=task$train,family=model$family,prior=model$prior,
    chains=s$chains,iter=s$iter,warmup=s$warmup,cores=s$cores,seed=task$metadata$seed,
    init=0,backend="cmdstanr",save_pars=brms::save_pars(all=TRUE),file=NULL,
    refresh=if(is.null(s$refresh))200L else s$refresh,
    control=list(adapt_delta=s$adapt_delta,max_treedepth=s$max_treedepth))
}

recovery_task <- function(t) {
  if(is.null(t$recovery_sampling))stop('No sampling recovery is configured')
  out<-t;out$metadata$sampling<-t$recovery_sampling;out$recovery_sampling<-NULL
  out
}
rejection_path <- function(t,retry=FALSE)sub('[.]rds$',if(retry)'.retry1.rejected.rds' else '.rejected.rds',t$path)
stop_sampler_rejection <- function(message) {
  stop(structure(list(message=message,call=NULL),class=c('aprime_sampler_rejection','error','condition')))
}
divergence_only_failure <- function(d) {
  if(!is.data.frame(d)||nrow(d)!=1L||!is.finite(d$n_divergent)||d$n_divergent<=0)return(FALSE)
  other<-d;other$n_divergent<-0L
  isTRUE(tryCatch(accept_sampler_diagnostics(other),error=function(e)FALSE))
}
read_initial_rejection <- function(t,identity) {
  e<-read_model_envelope(rejection_path(t),identity,expected_accepted=FALSE)
  base<-t;base$recovery_sampling<-NULL
  validate_task_fit(e,base)
  if(!is.null(e$diagnostic_error))stop_sampler_rejection(paste('Stored sampler diagnostics require investigation:',e$diagnostic_error))
  if(!identical(e$diagnostics,e$actual_diagnostics))stop('Stored rejection diagnostics differ from actual draws')
  e
}
sample_fit_task <- function(t) {
    args <- fit_task_arguments(t)
    message('Fitting ',t$kind,if(t$fold==0L)' full HC' else paste0(' held-out fold ',t$fold))
    do.call(brms::brm,args)
}
execute_fit_task <- function(t,identity,verify_only=FALSE) {
  if(file.exists(rejection_path(t,TRUE))&&(!file.exists(rejection_path(t))||file.exists(t$path)))
    stop('Inconsistent retry cache state; preserve evidence and investigate: ',t$path)
  if(!file.exists(t$path)) {
    if(verify_only)stop('Missing posterior fit: ',t$path)
    attempt<-t;retry<-FALSE
    if(file.exists(rejection_path(t))) {
      initial<-read_initial_rejection(t,identity)
      if(is.null(t$recovery_sampling)||!divergence_only_failure(initial$actual_diagnostics))
        stop_sampler_rejection(paste('Stored diagnostics require investigation:',rejection_path(t)))
      if(file.exists(rejection_path(t,TRUE))) {
        failed<-read_model_envelope(rejection_path(t,TRUE),identity,expected_accepted=FALSE)
        validate_task_fit(failed,recovery_task(t))
        stop_sampler_rejection(paste('The single sampling retry already failed; no automatic repetition:',rejection_path(t,TRUE)))
      }
      attempt<-recovery_task(t);retry<-TRUE
    }
    repeat {
    if(retry)message('One bounded divergence retry: ',basename(t$path),'; warmup=',attempt$metadata$sampling$warmup,
      ', adapt_delta=',attempt$metadata$sampling$adapt_delta)
    fit <- sample_fit_task(attempt)
    diagnostics <- tryCatch(actual_fit_diagnostics(fit,attempt$metadata$sampling,enforce=FALSE),error=function(e)e)
    accepted <- !inherits(diagnostics,'error') &&
      isTRUE(tryCatch(accept_sampler_diagnostics(diagnostics),error=function(e)FALSE))
    if(!accepted) {
      rejected<-rejection_path(t,retry)
      write_model_envelope(fit,attempt$metadata,identity,
        if(inherits(diagnostics,'error'))data.frame(error=conditionMessage(diagnostics)) else diagnostics,
        FALSE,rejected)
      if(!retry&&!is.null(t$recovery_sampling)&&!inherits(diagnostics,'error')&&divergence_only_failure(diagnostics)) {
        initial<-read_initial_rejection(t,identity)
        attempt<-recovery_task(t);retry<-TRUE;next
      }
      stop_sampler_rejection(paste('Fit diagnostics failed; retained rejected fit for inspection:',rejected))
    }
    validate_task_fit(list(fit=fit,metadata=attempt$metadata),attempt)
    write_model_envelope(fit,attempt$metadata,identity,diagnostics,TRUE,t$path)
    break
    }
  }
  e <- read_model_envelope(t$path,identity)
  validate_task_fit(e,t)
  if(!identical(e$metadata,t$metadata)) {
    initial<-read_initial_rejection(t,identity)
    if(!divergence_only_failure(initial$actual_diagnostics))stop('Retry lacks a verified divergence-only trigger')
  }
  # Workers return small audit records; posterior draws remain in their own RDS.
  e$fit <- NULL
  e
}

fit_all_models <- function(prepared,config,run_dir,identity,verify_only=FALSE) {
  tasks <- fit_tasks(prepared,config,run_dir)
  # Preserve each accepted model before continuing; aggregate diagnostics only
  # after every independent full/held-out fit has been checked.
  work <- function(t)tryCatch(execute_fit_task(t,identity,verify_only),
    aprime_sampler_rejection=function(e)list(sampler_rejected=TRUE,error=conditionMessage(e)))
  envelopes <- lapply(tasks,work)
  rejected<-vapply(envelopes,function(e)isTRUE(e$sampler_rejected),logical(1))
  if(any(rejected))stop('Sampler diagnostics remain rejected after independent fits finished: ',
    paste(vapply(envelopes[rejected],function(e)e$error,character(1)),collapse='; '),call.=FALSE)
  if(verify_only)return(invisible(do.call(rbind,Map(function(t,e)cbind(kind=t$kind,fold=t$fold,e$actual_diagnostics),tasks,envelopes))))
  # Only the parent writes these aggregates, in the original kind/fold order.
  write_model_manifest(tasks,envelopes,run_dir)
}

validate_all_models <- function(prepared,config,run_dir,identity) {
  fit_all_models(prepared,config,run_dir,identity,verify_only=TRUE)
}
