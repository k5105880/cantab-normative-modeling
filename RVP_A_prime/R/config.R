# Runtime allocation only: chains, draws, seeds and model definitions are unchanged.
# User-readable defaults live in ../config.json; descriptions are in ../../docs/SETTINGS.md.
# The JSON is intentionally strict: four chains, fixed quadrature nodes, at most
# one divergence-only retry, >=40000 direct prior draws and >=2000 bootstraps.
# analysis is a provenance label, not a model selector. iter includes warmup.
# scoring uses all retained draws; prior.draws instead counts independent direct draws.
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
validate_config <- function(config) {
  invalid <- function(field, message) {
    stop("Invalid configuration field '", field, "': ", message, call. = FALSE)
  }
  check_fields <- function(value, required, field, optional = character()) {
    if (!is.list(value) || is.object(value) || !is.null(dim(value)))
      invalid(field, "must be a named JSON object")
    keys <- names(value)
    if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)))
      invalid(field, "must have nonempty field names")
    if (anyDuplicated(keys))
      invalid(field, paste("duplicate field(s):", paste(unique(keys[duplicated(keys)]), collapse = ", ")))
    unknown <- setdiff(keys, c(required, optional))
    if (length(unknown))
      invalid(field, paste("unknown field(s):", paste(unknown, collapse = ", ")))
    missing <- setdiff(required, keys)
    if (length(missing))
      invalid(field, paste("missing required field(s):", paste(missing, collapse = ", ")))
  }
  plain_numeric <- function(value) {
    typeof(value) %in% c("integer", "double") && !is.object(value) && is.null(dim(value))
  }
  check_number <- function(value, field) {
    if (!plain_numeric(value) || length(value) != 1L || !is.finite(value))
      invalid(field, "must be one finite numeric value")
  }
  check_integer <- function(value, field, minimum = 1) {
    check_number(value, field)
    if (value != floor(value) || value < minimum || value > .Machine$integer.max)
      invalid(field, paste0("must be an integer from ", minimum, " to ", .Machine$integer.max))
  }
  check_string <- function(value, field) {
    if (!is.character(value) || is.object(value) || !is.null(dim(value)) ||
        length(value) != 1L || is.na(value) || !nzchar(value))
      invalid(field, "must be one nonempty character string")
  }

  check_fields(config, c("schema_version", "analysis", "seed", "sampling", "scoring",
                        "prior", "calibration_bootstraps"), "config", "sampling_recovery")
  check_fields(config$sampling, c("chains", "iter", "warmup", "cores", "adapt_delta",
                                 "max_treedepth", "refresh"), "sampling")
  check_fields(config$scoring, c("nodes", "tolerance", "mcse_limit", "workers", "chunk"), "scoring")
  check_fields(config$prior, c("draws", "batch_size"), "prior")
  check_integer(config$schema_version, "schema_version")
  if (config$schema_version != 1)
    invalid("schema_version", "must equal 1")
  check_string(config$analysis, "analysis")
  check_integer(config$seed, "seed")
  if (config$seed >= 2e9)
    invalid("seed", "must be less than 2000000000")
  check_integer(config$calibration_bootstraps, "calibration_bootstraps", 2000)

  s <- config$sampling
  for (name in c("chains", "iter", "warmup", "cores", "max_treedepth", "refresh"))
    check_integer(s[[name]], paste0("sampling.", name))
  check_number(s$adapt_delta, "sampling.adapt_delta")
  if (s$chains != 4)
    invalid("sampling.chains", "must equal 4")
  if (s$iter - s$warmup < 1000)
    invalid("sampling.iter", "must exceed sampling.warmup by at least 1000")
  if (s$cores > s$chains)
    invalid("sampling.cores", "must not exceed sampling.chains")
  if (s$adapt_delta < .99 || s$adapt_delta >= 1)
    invalid("sampling.adapt_delta", "must be at least 0.99 and less than 1")
  if (s$max_treedepth < 15)
    invalid("sampling.max_treedepth", "must be at least 15")

  scoring <- config$scoring
  if (!plain_numeric(scoring$nodes) ||
      !identical(as.numeric(scoring$nodes), c(32, 64, 128, 256)))
    invalid("scoring.nodes", "must be the numeric vector [32, 64, 128, 256]")
  for (name in c("workers", "chunk"))
    check_integer(scoring[[name]], paste0("scoring.", name))
  check_number(scoring$tolerance, "scoring.tolerance")
  if (scoring$tolerance <= 0 || scoring$tolerance > .002)
    invalid("scoring.tolerance", "must be greater than 0 and at most 0.002")
  check_number(scoring$mcse_limit, "scoring.mcse_limit")
  if (scoring$mcse_limit <= 0 || scoring$mcse_limit > .03)
    invalid("scoring.mcse_limit", "must be greater than 0 and at most 0.03")

  check_integer(config$prior$draws, "prior.draws", 40000)
  check_integer(config$prior$batch_size, "prior.batch_size")
  if (config$prior$draws %% config$prior$batch_size != 0)
    invalid("prior.draws", "must be divisible by prior.batch_size")

  # Both omission and JSON null preserve the existing no-retry behavior.
  if (!is.null(config$sampling_recovery)) {
    p <- config$sampling_recovery
    check_fields(p, c("trigger", "max_retries", "warmup", "adapt_delta"), "sampling_recovery")
    check_string(p$trigger, "sampling_recovery.trigger")
    if (p$trigger != "divergence_only")
      invalid("sampling_recovery.trigger", "must equal 'divergence_only'")
    check_integer(p$max_retries, "sampling_recovery.max_retries")
    if (p$max_retries != 1)
      invalid("sampling_recovery.max_retries", "must equal 1")
    check_integer(p$warmup, "sampling_recovery.warmup")
    if (p$warmup <= s$warmup)
      invalid("sampling_recovery.warmup", "must exceed sampling.warmup")
    if (as.double(p$warmup) + s$iter - s$warmup > .Machine$integer.max)
      invalid("sampling_recovery.warmup", "retry total iterations must not exceed the integer limit")
    check_number(p$adapt_delta, "sampling_recovery.adapt_delta")
    if (p$adapt_delta <= s$adapt_delta || p$adapt_delta >= 1)
      invalid("sampling_recovery.adapt_delta", "must exceed sampling.adapt_delta and be less than 1")
  }
  invisible(TRUE)
}

aprime_resource_plan <- function(config) {
  # Sampling and scoring caps are separate; Windows scoring runs serially.
  budget <- cantab_job_cores()
  list(job_cores=budget,
    chain_cores=as.integer(min(config$sampling$cores,config$sampling$chains,budget)),
    scoring_workers=cantab_worker_cores(config$scoring$workers))
}

apply_runtime_budget <- function(config) {
  plan <- aprime_resource_plan(config)
  config$sampling$cores <- plan$chain_cores
  config$scoring$workers <- plan$scoring_workers
  config
}
