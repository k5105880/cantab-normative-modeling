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
# Thread limits are applied only during execution, not while reading settings.
invisible(cantab_job_cores())

# Default sampling settings for the normative model.
# See ../docs/SETTINGS.md for units, constraints and the full settings dictionary.
# iter includes warmup per chain. Four chains are required; cores only caps concurrency.
# Divergence-only retry: warmup >=2500 and adapt_delta >=.9999, same retained draws.
# init_radius is the Stan initialization range; refresh is progress-message frequency.
# ndraws controls sampled probability-mass diagnostics; normative z uses all draws.
# fold_seed/fold_candidates build the fixed five participant-held-out folds.
# prediction_chunk_subjects batches non-control participants only (not control folds).
NORMATIVE_CONFIG <- list(chains=4L,iter=4000L,warmup=1500L,cores=1L,adapt_delta=.999,
 max_treedepth=12L,init_radius=.1,refresh=25L,ndraws=8000L,seed=20260921L,
 fold_seed=1334L,fold_candidates=500L,prediction_chunk_subjects=25L)
# Use run_cantab() or run_outcome() for input/output paths. A changed input,
# configuration, implementation or environment requires a new output directory.

NORMATIVE_CONFIG$cores <- cantab_worker_cores(NORMATIVE_CONFIG$cores)
