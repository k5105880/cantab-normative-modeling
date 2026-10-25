# Output protection applies to the common launcher and all six direct entry points.
local({
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  source(file.path(root, "R/input_data.R"), local = TRUE)
  source(file.path(root, "run_analysis.R"), local = TRUE)
  work <- tempfile("output_path_checks_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  repo <- file.path(work, "working_directory")
  dir.create(repo)
  rejected <- function(expr, pattern) {
    result <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  for (path in list(NA_character_, "", " ", character(), c("one", "two"), 3))
    rejected(cantab_output_directory(path, repo), "Supply one output")
  # POSIX treats drive-prefixed paths as relative; accepting them as absolute
  # could put a run inside the code repository when it is the working directory.
  if (.Platform$OS.type != "windows") {
    for (path in c("C:/new_output", "C:\\new_output", "z:/", "z:\\"))
      rejected(cantab_output_directory(path, repo), "Windows drive paths")
  }
  for (path in c(repo, file.path(repo, "new", "nested"),
                 file.path(work, "output", "..", "working_directory", "new")))
    rejected(cantab_output_directory(path, repo), "outside the public repository")
  allowed <- file.path(work, "output", "run")
  stopifnot(identical(cantab_output_directory(allowed, repo), allowed),
    identical(cantab_output_directory(file.path(repo, "..", "output", "run"), repo), allowed),
    !dir.exists(allowed))
  file <- file.path(work, "a_file")
  writeLines("synthetic", file)
  rejected(cantab_output_directory(file.path(file, "run"), repo), "component is a file")
  alias <- file.path(work, "alias")
  if (isTRUE(suppressWarnings(file.symlink(repo, alias)))) {
    rejected(cantab_output_directory(file.path(alias, "nested", "run"), repo),
             "outside the public repository")
    stopifnot(identical(cantab_output_directory(file.path(alias, "..", "output", "run"), repo),
                        allowed))
  }
  input <- file.path(work, "input.csv")
  writeLines("synthetic header only", input)
  outcomes <- c("RTI_reaction", "RTI_movement", "RVP_A_prime", "RVP_latency",
                "SWM_Between_errors", "SWM_Strategy")
  forbidden <- file.path(root, "must_not_create_test_output")
  stopifnot(!file.exists(forbidden))
  for (outcome in outcomes) {
    e <- new.env(parent = globalenv())
    source(file.path(root, outcome, "run.R"), local = e)
    rejected(e$run_outcome(input, forbidden, mode = "check"),
             "outside the public repository")
    rejected(run_cantab(outcome, input, forbidden, mode = "check"),
             "outside the public repository")
  }
  stopifnot(!file.exists(forbidden))
  cat("PASS: path normalization, nested paths, symlinks, invalid paths and all seven entry points reject output inside the code repository before writing.\n")
})
