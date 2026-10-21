# Runtime allocation only: chains, draws, seeds and model definitions are unchanged.
cantab_job_cores <- function() {
  raw <- Sys.getenv('CANTAB_JOB_CORES', unset='1')
  value <- suppressWarnings(as.numeric(raw))
  if (!grepl('^[1-9][0-9]*$', raw) || !is.finite(value) || value > .Machine$integer.max)
    stop('CANTAB_JOB_CORES must be a positive integer')
  as.integer(value)
}
cantab_worker_cores <- function(requested=1L) {
  if (length(requested)!=1L || !is.finite(requested) || requested<1 || requested!=floor(requested))
    stop('Requested worker count must be a positive integer')
  as.integer(min(requested, cantab_job_cores()))
}
# Identifiable common time smooth with optional platform-specific linear slopes.
cantab_time_term <- function(smooth, platform_slopes) {
  if (!smooth) return('Time_since_baseline_1y')
  if (platform_slopes) return('s(Time_since_baseline_1y, k = 4, bs = "tp", m = c(2, 0))')
  # Without platform slopes, the ordinary smooth supplies the common linear effect.
  's(Time_since_baseline_1y, k = 4)'
}

cantab_brm <- function(formula, data, family, prior, ...) {
  dots <- list(...)
  sargs <- dots[intersect(names(dots), c('sample_prior','data2','knots','drop_unused_levels','stanvars'))]
  sd <- do.call(brms::make_standata, c(list(formula=formula,data=data,family=family,prior=prior),sargs))
  label <- if (!is.null(dots$file)) as.character(dots$file) else 'full_fit'
  cantab_audit_standata(sd, label, file.path(getwd(),'design_rank_audit.csv'))
  dots$cores <- cantab_worker_cores(if (is.null(dots$cores)) 1L else dots$cores)
  dots$threads <- brms::threading(1)
  do.call(brms::brm,c(list(formula=formula,data=data,family=family,prior=prior),dots))
}

cantab_kfold <- function(x, K, folds, ...) {
  # brms::kfold updates the parent formula with these exact training rows. Rebuild
  # the same training bases (no full-data basis reuse, no held-out outcome input).
  folds <- as.integer(factor(folds))
  if (length(folds)!=nrow(x$data) || anyNA(folds)) stop('Invalid fold vector')
  drop_unused <- attr(x$data,'drop_unused_levels',exact=TRUE)
  if (is.null(drop_unused)) drop_unused <- TRUE
  for (k in sort(unique(folds))) {
    train <- x$data[folds!=k,,drop=FALSE]
    test <- x$data[folds==k,,drop=FALSE]
    if (any(as.character(train$Subject_ID) %in% as.character(test$Subject_ID)))
      stop('Subject leakage in training design audit: fold ',k)
    sd <- brms::make_standata(x$formula,data=train,family=x$family,prior=x$prior,
      knots=attr(x$data,'knots',exact=TRUE),drop_unused_levels=drop_unused,
      stanvars=x$stanvars)
    label <- paste('CV',x$family$family,'K',max(folds),'fold',k,
      digest::digest(list(x$formula,folds,train$row_id_global),algo='xxhash64'))
    cantab_audit_standata(sd,label,file.path(getwd(),'design_rank_audit.csv'))
  }
  # Folds are sequential; all four chains run within the job budget, one Stan thread each.
  future::plan(future::sequential)
  dots <- list(...);dots$cores <- cantab_worker_cores(if (is.null(dots$cores)) 1L else dots$cores)
  dots$threads <- brms::threading(1)
  do.call(brms::kfold,c(list(x=x,K=K,folds=folds),dots))
}
