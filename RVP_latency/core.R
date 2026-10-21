# Sourceable context and engine for the RVP_latency normative model.
normative_context <- local({
  paths <- unlist(lapply(sys.frames(),function(f)f$ofile),use.names=FALSE)
  root <- dirname(normalizePath(tail(paths,1L),mustWork=TRUE))
  function(input_file,output_dir,prepare_run=FALSE) {
    e <- new.env(parent=globalenv())
    source(file.path(dirname(root),'R/input_data.R'),local=e)
    output_dir <- e$cantab_output_directory(output_dir,dirname(root))
    input_file <- normalizePath(input_file,mustWork=TRUE)
    if(prepare_run)dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    for(package in c('dplyr','brms','posterior','loo'))
      suppressPackageStartupMessages(library(package,character.only=TRUE))
    for(file in c('engine_helpers.R','rvp_validation_helpers.R',
      'cantab_design_audit.R','identifiability_runtime.R','state_cdf.R','model_settings.R',
      'prepare_data.R','model_functions.R','pipeline.R','audit_run.R'))source(file.path(root,file),local=e)
    e$engine_check_packages()
    e$cores_fit <- e$cantab_worker_cores(e$chains_fit);e$cv_cores <- e$cantab_worker_cores(e$cv_chains)
    e$script_dir <- root;e$infile <- input_file;e$run_dir <- output_dir
    e$analysis_version <- 'RVP_latency_normative';e$model_cache_version <- 'normative'
    e$output_prefix <- file.path(output_dir,'RVP_latency')
    e$checkpoint_dir <- file.path(output_dir,'model_checkpoints')
    e$cv_checkpoint_dir <- file.path(output_dir,'cv_checkpoints')
    e$verbose_checks <- TRUE
    e$movement_outcome <- 'RVP_Median_latency'
    e$outcome_specs <- list(list(raw_col='RVP_Median_latency',outcome_prefix='RVP_Median_latency',
      primary_family='gamma_raw',higher_is_better=FALSE,analysis_role='primary',
      support_class='positive_continuous_raw'))
    prepared <- e$prepare_normative_data(input_file)
    e$data <- prepared$data;e$control_label <- prepared$control_label
    e
  }
})
run_normative_engine <- function(input_file,output_dir) {
  e <- normative_context(input_file,output_dir,prepare_run=TRUE)
  previous_dir <- setwd(e$run_dir);on.exit(setwd(previous_dir),add=TRUE)
  previous_options <- options();on.exit(options(previous_options),add=TRUE)
  e$engine_guard(e$script_dir,e$infile,e$run_dir)
  if(e$engine_complete(e$run_dir))return(invisible(e$run_dir))
  if(grepl('=',tempdir(),fixed=TRUE))stop('R temporary directory must not contain an equals sign')
  cmdstanr::check_cmdstan_toolchain(fix=FALSE,quiet=TRUE)
  cache <- file.path(e$run_dir,'stan_cache');dir.create(cache,showWarnings=FALSE)
  options(brms.backend='cmdstanr',cmdstanr_write_stan_file_dir=cache)
  dir.create(e$checkpoint_dir,showWarnings=FALSE)
  e$run_normative_pipeline()
  e$audit_normative_run(e$run_dir,file.path(e$run_dir,'audit'),e$infile)
  e$engine_mark_complete(e$run_dir)
  invisible(e$run_dir)
}
