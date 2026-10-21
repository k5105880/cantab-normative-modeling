# The normal interface needs no CPU setting; four chains still run serially.
local({
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  before <- Sys.getenv("CANTAB_JOB_CORES", unset = NA_character_)
  on.exit(if (is.na(before)) Sys.unsetenv("CANTAB_JOB_CORES") else
    Sys.setenv(CANTAB_JOB_CORES = before), add = TRUE)
  Sys.unsetenv("CANTAB_JOB_CORES")
  source(file.path(root, "run_analysis.R"), local = TRUE)
  source(file.path(root, "environment/setup.R"), local = TRUE)
  stopifnot(!"cores" %in% names(formals(run_cantab)),
            !"cores" %in% names(formals(install_cantab_packages)))
  outcomes <- c("RTI_reaction", "RTI_movement", "RVP_A_prime", "RVP_latency",
                "SWM_Between_errors", "SWM_Strategy")
  for (outcome in outcomes) {
    e <- new.env(parent = globalenv())
    source(file.path(root, outcome, "run.R"), local = e)
    stopifnot(!"cores" %in% names(formals(e$run_outcome)))
    runtime <- switch(outcome,
      RTI_movement = "identifiability_runtime.R", RVP_latency = "identifiability_runtime.R",
      RVP_A_prime = "R/config.R", "config.R")
    source(file.path(root, outcome, runtime), local = e)
    stopifnot(e$cantab_job_cores() == 1L, e$cantab_worker_cores(4L) == 1L)
    if (outcome %in% c("RTI_movement", "RVP_latency")) {
      source(file.path(root, outcome, "model_settings.R"), local = e)
      stopifnot(e$chains_fit == 4L, e$cv_chains == 4L,
                e$cantab_worker_cores(e$chains_fit) == 1L,
                e$cantab_worker_cores(e$cv_chains) == 1L)
    } else if (outcome == "RVP_A_prime") {
      cfg <- jsonlite::fromJSON(file.path(root, outcome, "config.json"))
      e$validate_config(cfg)
      effective <- e$apply_runtime_budget(cfg)
      stopifnot(cfg$sampling$cores == 1L, cfg$scoring$workers == 1L,
                effective$sampling$cores == 1L, effective$scoring$workers == 1L,
                effective$sampling$chains == 4L)
      source(file.path(root, outcome, "R/precision_helpers.R"), local = e)
      stopifnot(formals(e$score_cdf_quadrature)$workers == 1L)
    } else {
      cfg <- switch(outcome, RTI_reaction = e$rti_config,
        SWM_Between_errors = e$CONFIG, SWM_Strategy = e$NORMATIVE_CONFIG)
      stopifnot(cfg$cores == 1L, cfg$chains == 4L)
    }
  }
  cat("PASS: no public core setting; six runtime defaults and scoring use one worker while preserving four sampling chains.\n")
})
