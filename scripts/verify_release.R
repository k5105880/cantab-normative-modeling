# Source this file to verify a pristine public release using R.
verify_cantab_release <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  for (pkg in c("digest", "jsonlite"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install the R package ", pkg)
  manifest <- jsonlite::fromJSON(file.path(root, "RELEASE_MANIFEST.json"))
  expected <- manifest$files
  actual <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE, include.dirs = FALSE)
  actual <- actual[!grepl("(^|/)\\.git(/|$)", actual) & actual != "RELEASE_MANIFEST.json"]
  if (!setequal(actual, expected$path))
    stop("Release inventory differs. Unexpected: ", paste(setdiff(actual, expected$path), collapse = ", "),
         "; missing: ", paste(setdiff(expected$path, actual), collapse = ", "))
  for (i in seq_len(nrow(expected))) {
    path <- file.path(root, expected$path[i])
    if (file.info(path)$size != expected$bytes[i] ||
        digest::digest(file = path, algo = "sha256", serialize = FALSE) != expected$sha256[i])
      stop("Changed release file: ", expected$path[i])
  }
  sources <- read.csv(file.path(root, "provenance", "source_manifest.csv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sources)))
    if (digest::digest(file = file.path(root, sources$path[i]), algo = "sha256",
                       serialize = FALSE) != sources$sha256[i])
      stop("Source inventory mismatch: ", sources$path[i])
  files <- list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  for (path in files) parse(file = path)
  message("PASS: ", nrow(expected), " release entries; ", length(files), " R files parsed.")
  invisible(TRUE)
}
.verify_files <- as.character(unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE))
.verify_self <- .verify_files[basename(.verify_files) == "verify_release.R"]
if (!length(.verify_self)) {
  .verify_self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  .verify_self <- gsub("~+~", " ", .verify_self, fixed = TRUE)
}
verify_cantab_release(dirname(dirname(normalizePath(tail(.verify_self, 1L), mustWork = TRUE))))
rm(.verify_files, .verify_self)
