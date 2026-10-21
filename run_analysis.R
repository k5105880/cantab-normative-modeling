# Source this file in R or RStudio. Working directory may be anywhere.
.cantab_source_files <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
.cantab_source_files <- .cantab_source_files[basename(.cantab_source_files) == "run_analysis.R"]
if (!length(.cantab_source_files)) stop("Use source('path/to/run_analysis.R') from R.")
.cantab_repository <- dirname(normalizePath(tail(.cantab_source_files, 1L), winslash = "/"))
source(file.path(.cantab_repository, "R", "input_data.R"), local = TRUE)

run_cantab <- local({
  repository <- .cantab_repository
  outcomes <- c("RTI_reaction", "RTI_movement", "RVP_A_prime", "RVP_latency",
                "SWM_Between_errors", "SWM_Strategy")
  function(outcome, input_file, output_dir, mode = c("check", "fresh", "resume"),
           reference_dir = NULL) {
    if (length(outcome) != 1L || !outcome %in% outcomes)
      stop("outcome must be one of: ", paste(outcomes, collapse = ", "))
    mode <- match.arg(mode)
    input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
    if (!is.null(reference_dir))
      reference_dir <- normalizePath(reference_dir, winslash = "/", mustWork = TRUE)
    destination <- cantab_output_directory(output_dir, repository)
    env <- new.env(parent = globalenv())
    source(file.path(repository, outcome, "run.R"), local = env)
    env$run_outcome(input_file = input_file, output_dir = destination, mode = mode,
      reference_dir = reference_dir)
  }
})
rm(.cantab_source_files, .cantab_repository)
