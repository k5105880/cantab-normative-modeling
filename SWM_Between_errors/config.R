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
# Normative model settings. Changes require a fresh run.
# See ../docs/SETTINGS.md for every setting, units and fixed constraints.
CONFIG <- list(
 # Diagnosis labels and platform exposure mapping are part of the data/model contract.
 control_label='HC', groups=c('HC','MDD','BD'),
 trials=c(RS6=4,Connect=1),
 # Partition participants, keeping all their visits together.
 folds=5L, fold_seed=1334L, fold_candidates=500L,
 # Only this fixed model ID is supported; not a list of alternative models.
 models='normative',
 # Counts are per chain; iter includes warmup. cores caps concurrent chains.
 chains=4L, cores=1L, iter=3000L, warmup=1000L,
 # One retry after posterior sampling diagnostics fail.
 retry_iter=5000L, retry_warmup=1500L,
 adapt_delta=.995, retry_adapt_delta=.999,
 max_treedepth=18L, retry_max_treedepth=20L,
 # Sampled interval-coverage diagnostics only; normative z uses all retained draws.
 prediction_draws=8000L,
 # Subject bootstrap repetitions; optional prior fit uses 1200/600 total/warmup.
 bootstrap_reps=2000L, prior_predictive=TRUE,
 # Fixed role-specific offsets are added for fitting/prediction/bootstrap.
 seed=202609210L
)

CONFIG$cores <- cantab_worker_cores(CONFIG$cores)
