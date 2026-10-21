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

# Fixed normative model and reproducible sampling settings.
# See ../docs/SETTINGS.md for units, constraints and the full settings dictionary.
rti_config <- list(
  # version is a provenance label. folds are IDs; the workflow requires all five.
  version = 'reaction_normative_v1', seed = 20260921L, folds = 1:5,
  # cores limits concurrent chains; it does not reduce the four fitted chains.
  chains = 4L, cores = 1L,
  # Total iterations per chain INCLUDE warmup; warmup draws are discarded.
  full_iter = 5000L, full_warmup = 1500L,
  cv_iter = 3000L, cv_warmup = 1000L,
  # One diagnostic retry; the multiplier increases both iterations and warmup.
  adapt_delta = .999, retry_adapt_delta = .99999, retry_iteration_multiplier = 2L,
  max_treedepth = 20L,
  # Acceptance gates: R-hat strictly below max_rhat; bulk/tail ESS and BFMI >= minima.
  max_rhat = 1.01, min_ess = 400, min_bfmi = .3,
  # All ordered post-warmup draws are used. MCSE is simulation error of z.
  posterior_mcse_z_limit = .03, quad_nodes = 40L, check_nodes = 80L,
  # Deterministic integration refinement limits, separate from posterior MCSE.
  quad_z_tolerance = .005, quad_log_score_tolerance = .0001, quad_max_nodes = 320L,
  # Rows scored together, not posterior draws or participants.
  cdf_chunk_size = 40L,
  # prior_draws controls predictive summaries. Prior fit: 2 chains, 1500/750 total/warmup.
  run_prior_predictive = TRUE, prior_draws = 1000L
)
rti_config$cores <- cantab_worker_cores(rti_config$cores)
