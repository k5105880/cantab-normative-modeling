# Source this file, then call run_outcome(). Sourcing never runs an analysis.
run_outcome <- local({
  sources <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
  sources <- sources[basename(sources) == 'run.R']
  if (!length(sources)) stop('Load this entry point with source("RTI_movement/run.R").')
  root <- dirname(normalizePath(tail(sources, 1L), mustWork = TRUE))
  outcome <- 'RTI_movement'

  function(input_file, output_dir, mode = c('check', 'fresh', 'resume'),
           reference_dir = NULL) {
    previous_options <- options()
    on.exit({
      added <- setdiff(names(options()),names(previous_options))
      if(length(added))options(setNames(rep(list(NULL),length(added)),added))
      options(previous_options)
    },add=TRUE)
    mode <- match.arg(mode)
    scalar_path <- function(value, label) {
      if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
        stop(label, ' must be one nonempty path.')
      path.expand(value)
    }
    input_file <- normalizePath(scalar_path(input_file, 'input_file'), mustWork = TRUE)
    if (dir.exists(input_file)) stop('input_file must be a CSV file.')
    output_dir <- scalar_path(output_dir, 'output_dir')
    # One worker; all four sampling chains and iteration counts are preserved.
    cores <- 1L
    if (!is.null(reference_dir)) {
      reference_dir <- normalizePath(scalar_path(reference_dir, 'reference_dir'), mustWork = TRUE)
      if (!dir.exists(reference_dir)) stop('reference_dir must be a directory.')
    }
    packages <- c('digest', 'callr')
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', '), '. Source environment/setup.R, then call install_cantab_packages().')
    helper <- file.path(dirname(root), 'R', 'input_data.R')
    if (!file.exists(helper)) stop('Shared input helper is missing: ', helper)
    helper_env <- new.env(parent = globalenv())
    source(helper, local = helper_env)
    output_dir <- helper_env$cantab_output_directory(output_dir, dirname(root))
    results_helper <- file.path(dirname(root), 'R', 'results.R')
    source(results_helper, local = helper_env)
    if (!is.function(helper_env$cantab_prepare_input)) stop('Shared input helper is incomplete.')
    if (mode == 'resume') {
      if (!dir.exists(output_dir)) stop('Resume requires an existing output directory.')
    } else {
      if (file.exists(output_dir)) stop('Use a new output directory for ', mode, ': ', output_dir)
      if (!dir.create(output_dir, recursive = TRUE)) stop('Cannot create output directory: ', output_dir)
    }
    output_dir <- normalizePath(output_dir, winslash = '/', mustWork = TRUE)
    lock <- file.path(output_dir, '.run.lock')
    if (!dir.create(lock, showWarnings = FALSE))
      stop('This output is locked: ', lock, '. Confirm its previous R process has stopped before removing a stale lock.')
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    writeLines(c(paste('PID', Sys.getpid()), paste('Host', Sys.info()[['nodename']]),
                 as.character(Sys.time())), file.path(lock, 'owner.txt'))
    hash <- function(path) unname(vapply(path, function(p) digest::digest(file = p, algo = 'sha256'), character(1)))
    code <- sort(unique(c(list.files(root, pattern = '[.]R$', full.names = TRUE),
      unlist(lapply(file.path(root, c('R', 'vendor', 'stan')), function(p)
        list.files(p, pattern = '[.](R|stan)$', recursive = TRUE, full.names = TRUE)), use.names = FALSE))))
    references <- if (is.null(reference_dir)) character() else
      sort(list.files(reference_dir, pattern = '^diagnosis_reference([.].*)?$', full.names = TRUE))
    identity <- list(schema = 'cantab-normative-v1', outcome = outcome,
      raw_input = input_file, raw_sha256 = hash(input_file), reference_dir = reference_dir,
      references = setNames(hash(references), basename(references)),
      code = setNames(hash(code), substring(code, nchar(root) + 2L)),
      input_helper_sha256 = hash(helper), results_helper_sha256 = hash(results_helper),
      R = R.version.string,
      platform = R.version$platform, cores = cores)
    manifest_path <- file.path(output_dir, 'run_identity.rds')
    analysis_input <- file.path(output_dir, 'analysis_input.csv')
    if (mode == 'resume') {
      if (!file.exists(manifest_path)) stop('No portable run manifest. Start a fresh run; a run identity is required.')
      previous <- readRDS(manifest_path)
      if (!identical(previous$kind, 'fitting')) stop('A check directory cannot be resumed as a fitting run.')
      if (!identical(previous$identity, identity)) stop('Input, reference, code, R environment or execution settings changed. Use a new output directory.')
      if (!file.exists(analysis_input) || !identical(hash(analysis_input), previous$analysis_sha256))
        stop('Prepared analysis input is missing or changed. Resume refused.')
    }
    prepared <- helper_env$cantab_prepare_input(input_file = input_file,
      output_dir = output_dir, reference_dir = reference_dir, outcome = outcome)
    analysis_input <- normalizePath(prepared, winslash = '/', mustWork = TRUE)
    if (mode == 'resume' && !identical(hash(analysis_input), previous$analysis_sha256))
      stop('Prepared analysis input differs from the saved run. Resume refused.')
    if (mode != 'resume') saveRDS(list(identity = identity,
      kind = if (mode == 'check') 'check' else 'fitting',
      analysis_sha256 = hash(analysis_input), created = Sys.time()), manifest_path)
    env <- c(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),
      CANTAB_JOB_CORES=as.character(cores),NM_INPUT_FILE=analysis_input,NM_RUN_DIR=output_dir,
      OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
      VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',STAN_NUM_THREADS='1',
      MAKEFLAGS='-j1',CMAKE_BUILD_PARALLEL_LEVEL='1')
    log <- file.path(output_dir, paste0('console_', mode, '_', format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid(), '.log'))
    status_file <- file.path(output_dir, 'run_status.txt')
    writeLines(c(paste('RUNNING', toupper(mode)), as.character(Sys.time())), status_file)
    message(outcome, ': ', mode, '; output: ', output_dir)
    args <- paste0('--', mode)
    scores <- NULL
    tryCatch({
      callr::rscript(file.path(root, 'main.R'), cmdargs = args, wd = output_dir,
        libpath = .libPaths(), env = env, stdout = log, stderr = '2>&1',
        show = TRUE, system_profile = FALSE, user_profile = FALSE)
      if (mode != 'check') scores <- helper_env$cantab_export_results(outcome, analysis_input, output_dir)
      writeLines(c(if (mode == 'check') 'CHECK_COMPLETE: no MCMC or normative scores.' else 'COMPLETE',
        paste('Outcome:', outcome), as.character(Sys.time())), status_file)
    }, interrupt = function(e) {
      writeLines(c('INTERRUPTED', 'Inspect the logs before resuming with unchanged input and code.'), status_file)
      stop(e)
    }, error = function(e) {
      writeLines(c('FAILED', conditionMessage(e), paste('Log:', log)), status_file)
      stop('Analysis failed. Inspect ', log, '\n', conditionMessage(e), call. = FALSE)
    })
    invisible(list(outcome = outcome, mode = mode, output_dir = output_dir,
      analysis_input = analysis_input, log = log, scores = scores,
      status = if (mode == 'check') 'checked' else 'completed'))
  }
})
