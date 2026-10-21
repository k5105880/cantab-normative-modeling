# Source from the repository root, or Rscript scripts/run_tests.R.
.test_files <- as.character(unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE))
.test_self <- .test_files[basename(.test_files) == "run_tests.R"]
if (!length(.test_self)) {
  .test_self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .test_self <- gsub("~+~", " ", .test_self, fixed = TRUE)
}
.test_root <- dirname(dirname(normalizePath(tail(.test_self, 1L), winslash = "/", mustWork = TRUE)))
if (!requireNamespace("callr", quietly = TRUE)) stop("Install the R package callr.")
.test_paths <- c("tests/test_input_data.R", readLines(file.path(.test_root, "scripts", "test_files.txt")))
.test_paths <- .test_paths[nzchar(.test_paths) & !startsWith(.test_paths, "#")]
for (.test_path in .test_paths) {
  message("TEST: ", .test_path)
  callr::rscript(file.path(.test_root, .test_path), wd = .test_root,
                 libpath = .libPaths(), show = TRUE,
                 system_profile = FALSE, user_profile = FALSE)
}
message("PASS: ", length(.test_paths), " public test scripts; synthetic data only; no production MCMC.")
rm(list = ls(pattern = "^\\.test_"))
