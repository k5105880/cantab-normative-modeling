# Public input preparation. Reading/creating a CSV uses base R only.

# Resolve an output path before creating anything. This also covers direct calls
# to an outcome entry point and paths passing through symbolic links.
cantab_output_directory <- function(output_dir, repository) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(trimws(output_dir)))
    stop("Supply one output directory outside the public repository.")
  repository <- normalizePath(repository, winslash = "/", mustWork = TRUE)
  destination <- path.expand(output_dir)
  if (.Platform$OS.type != "windows" && grepl("^[A-Za-z]:[/\\\\]", destination))
    stop("Windows drive paths are not absolute on this system; supply a local output path outside the public repository.")
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", destination))
    destination <- file.path(getwd(), destination)
  # Resolve existing components in order: symlink/.. must follow the real path.
  destination <- gsub("\\\\", "/", destination)
  anchor <- if (grepl("^[A-Za-z]:/", destination)) substr(destination, 1L, 3L) else
    if (startsWith(destination, "//")) {
      parts <- strsplit(sub("^//", "", destination), "/", fixed = TRUE)[[1]]
      if (length(parts) < 2L) stop("Incomplete network output path.")
      paste0("//", parts[1], "/", parts[2], "/")
    } else "/"
  parts <- strsplit(substring(destination, nchar(anchor) + 1L), "/", fixed = TRUE)[[1]]
  resolved <- anchor
  for (part in parts) {
    if (!nzchar(part) || part == ".") next
    if (part == "..") resolved <- dirname(resolved)
    else resolved <- file.path(resolved, part)
    if (file.exists(resolved)) {
      if (!dir.exists(resolved)) stop("An output path component is a file: ", resolved)
      resolved <- normalizePath(resolved, winslash = "/", mustWork = TRUE)
    } else if (!is.na(Sys.readlink(resolved)) && nzchar(Sys.readlink(resolved))) {
      stop("Output path contains an unresolved symbolic link: ", resolved)
    }
  }
  # Windows filesystem comparisons are case insensitive.
  compare <- function(x) if (.Platform$OS.type == "windows") tolower(x) else x
  if (identical(compare(resolved), compare(repository)) ||
      startsWith(compare(resolved), paste0(compare(repository), "/")))
    stop("Choose output_dir outside the public repository; run outputs contain participant data.")
  resolved
}

cantab_outcome_columns <- c(
  "SWM_Between_errors", "SWM_Strategy", "RVP_A_prime", "RVP_Median_latency",
  "RTI_Median_five-choice_reaction_time", "RTI_Median_five-choice_movement_time")
cantab_required_columns <- c("Subject_ID", "Diagnosis", "Timepoint", "Platform",
  "Gender", "FSIQ", "Age_baseline_years", "Time_since_baseline_years",
  cantab_outcome_columns)

cantab_read_characters <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  # read.csv can return partial rows after an unmatched quote, or invent a
  # row.names column when data records are wider than the header. Reject both.
  withCallingHandlers(
    invisible(scan(path, what = character(), sep = ",", quote = '"',
      comment.char = "", blank.lines.skip = FALSE, na.strings = NULL,
      fileEncoding = "UTF-8-BOM", quiet = TRUE)),
    warning = function(w) stop("Malformed CSV: ", conditionMessage(w), call. = FALSE))
  fields <- count.fields(path, sep = ",", quote = '"', comment.char = "",
                         blank.lines.skip = FALSE)
  fields <- fields[!is.na(fields)] # A quoted multiline field is one record.
  if (!length(fields) || length(unique(fields)) != 1L)
    stop("Malformed CSV: every record must have the same number of fields as the header.")
  x <- read.csv(path, check.names = FALSE, colClasses = "character",
    na.strings = c("", "NA", "NaN"), fileEncoding = "UTF-8-BOM",
    stringsAsFactors = FALSE, fill = FALSE, blank.lines.skip = FALSE,
    row.names = NULL)
  if (!nrow(x)) stop("Input CSV has no data rows.")
  if (anyDuplicated(names(x)) || any(!nzchar(names(x))) ||
      any(names(x) != trimws(names(x))))
    stop("Duplicate, empty or whitespace-padded column names.")
  x
}

read_cantab_input <- function(path) {
  x <- cantab_read_characters(path)
  # Accept documented alternative day columns, but always expose years to callers.
  units <- c(Age_baseline_years = "Age_baseline",
             Time_since_baseline_years = "Time_since_baseline")
  missing <- setdiff(cantab_required_columns,
    c(names(x), names(units)[unname(units) %in% names(x)]))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "),
    ". Preserve hyphens in the RTI names (check.names = FALSE).")
  original <- x
  if (anyNA(x$Subject_ID) || any(x$Subject_ID != trimws(x$Subject_ID)) ||
      any(grepl("[[:cntrl:]]", x$Subject_ID)))
    stop("Subject_ID must be nonempty text without surrounding spaces or control characters.")
  for (nm in intersect(c(cantab_required_columns, "Age_years"), names(x))) {
    if (any(!is.na(x[[nm]]) & x[[nm]] != trimws(x[[nm]])))
      stop(nm, " contains surrounding spaces; correct the source CSV before analysis.")
  }
  x[] <- lapply(x, trimws)
  x[] <- lapply(x, function(z) { z[z %in% c("", "NA", "NaN")] <- NA_character_; z })
  if (anyNA(x$Diagnosis) || any(!x$Diagnosis %in% c("HC", "BD", "MDD")))
    stop("Diagnosis must be HC, MDD or BD; labels are never automatically changed.")
  if (anyNA(x$Platform) || any(!x$Platform %in% c("RS6", "Connect")))
    stop("Platform must be exactly RS6 or Connect.")
  if (any(!is.na(x$Gender) & !x$Gender %in% c("F", "M")))
    stop("Gender must be F, M or missing.")
  numeric_names <- intersect(names(x), c(cantab_outcome_columns, names(units),
    unname(units), "Timepoint", "Age_years", "FSIQ", "row_id_global"))
  for (nm in numeric_names) {
    parsed <- suppressWarnings(as.numeric(x[[nm]]))
    if (any(!is.na(x[[nm]]) & (!is.finite(parsed) | is.na(parsed))))
      stop(nm, " must contain finite numbers or blank/NA missing values.")
    x[[nm]] <- parsed
  }
  for (nm in names(units)) {
    day <- units[[nm]]
    if (!nm %in% names(x)) x[[nm]] <- x[[day]] / 365.25
    else if (day %in% names(x)) {
      conflict <- xor(is.na(x[[nm]]), is.na(x[[day]])) |
        (!is.na(x[[nm]]) & !is.na(x[[day]]) & abs(x[[nm]] - x[[day]] / 365.25) > 1e-6)
      if (any(conflict)) stop("Conflicting day/year columns: ", nm)
    }
  }
  if (anyNA(x$Timepoint) || any(x$Timepoint < 1 | x$Timepoint != floor(x$Timepoint)))
    stop("Timepoint must be a positive integer.")
  if (anyDuplicated(paste(x$Subject_ID, x$Timepoint, sep = "\r")))
    stop("Duplicate Subject_ID/Timepoint. One participant has one row per test occasion.")
  if ("row_id_global" %in% names(x) &&
      (anyNA(x$row_id_global) || any(x$row_id_global < 1 |
       x$row_id_global != floor(x$row_id_global)) || anyDuplicated(x$row_id_global)))
    stop("row_id_global must contain unique positive integers.")
  for (nm in intersect(c("Age_baseline_years", "Age_years"), names(x))) {
    if (any(!is.na(x[[nm]]) & (x[[nm]] < 5 | x[[nm]] > 110)))
      stop(nm, " must be in years (supported input check: 5 to 110).")
  }
  if (any(!is.na(x$Time_since_baseline_years) &
          (x$Time_since_baseline_years < 0 | x$Time_since_baseline_years > 50)))
    stop("Time_since_baseline_years must be in years between 0 and 50.")
  if (any(!is.na(x$FSIQ) & x$FSIQ <= 0)) stop("FSIQ must be positive or missing.")
  for (nm in cantab_outcome_columns[4:6])
    if (any(!is.na(x[[nm]]) & x[[nm]] <= 0))
      stop(nm, " must be positive milliseconds or missing.")
  if (any(!is.na(x$RVP_A_prime) & (x$RVP_A_prime < 0 | x$RVP_A_prime > 1)))
    stop("RVP_A_prime must be within [0, 1], not a percentage.")
  if (any(!is.na(x$SWM_Between_errors) &
          (x$SWM_Between_errors < 0 | x$SWM_Between_errors != floor(x$SWM_Between_errors))))
    stop("SWM_Between_errors must be a nonnegative integer or missing.")
  for (ii in split(seq_len(nrow(x)), x$Subject_ID)) {
    z <- x[ii, , drop = FALSE]
    if (length(unique(z$Diagnosis)) != 1) stop("Diagnosis changes within a participant.")
    ages <- na.omit(z$Age_baseline_years)
    if (length(ages) && diff(range(ages)) > 1e-5)
      stop("Inconsistent Age_baseline_years within a participant.")
    for (nm in c("FSIQ", "Gender"))
      if (length(unique(na.omit(z[[nm]]))) > 1)
        stop("Inconsistent fixed baseline ", nm, " within a participant.")
    elapsed <- z$Time_since_baseline_years
    if (any(!is.na(elapsed) & ((z$Timepoint == 1 & abs(elapsed) > 1e-6) |
                               (z$Timepoint > 1 & elapsed <= 0))))
      stop("Baseline/Timepoint mismatch in Time_since_baseline_years.")
    observed <- na.omit(elapsed[order(z$Timepoint)])
    if (length(observed) > 1 && any(diff(observed) <= 0))
      stop("Elapsed times must increase with Timepoint within a participant.")
    if ("Age_years" %in% names(z)) {
      complete <- complete.cases(z[, c("Age_years", "Age_baseline_years",
                                     "Time_since_baseline_years")])
      if (any(abs(z$Age_years[complete] - z$Age_baseline_years[complete] -
                  elapsed[complete]) > 1e-5))
        stop("Age_years differs from baseline age plus elapsed years.")
    }
  }
  attr(x, "raw_characters") <- original
  x
}

# participants: one row per person; assessments: one row per test occasion.
create_cantab_input <- function(participants, assessments, output_file, overwrite = FALSE) {
  if (file.exists(output_file) && !isTRUE(overwrite)) stop("Output file already exists: ", output_file)
  p <- if (is.character(participants) && length(participants) == 1L)
    cantab_read_characters(participants) else as.data.frame(participants, stringsAsFactors = FALSE)
  a <- if (is.character(assessments) && length(assessments) == 1L)
    cantab_read_characters(assessments) else as.data.frame(assessments, stringsAsFactors = FALSE)
  p_names <- c("Subject_ID", "Diagnosis", "Gender", "FSIQ", "Age_baseline_years", "Baseline_date")
  a_names <- c("Subject_ID", "Assessment_date", "Timepoint", "Platform", cantab_outcome_columns)
  for (item in list(list(x = p, required = p_names, label = "participants"),
                    list(x = a, required = a_names, label = "assessments"))) {
    if (anyDuplicated(names(item$x))) stop("Duplicate columns in ", item$label)
    absent <- setdiff(item$required, names(item$x))
    if (length(absent)) stop("Missing ", item$label, " columns: ", paste(absent, collapse = ", "))
    if (!nrow(item$x)) stop(item$label, " has no data rows.")
  }
  if (anyNA(p$Subject_ID) || anyDuplicated(p$Subject_ID))
    stop("Missing or duplicate Subject_ID in participants.")
  match_id <- match(as.character(a$Subject_ID), as.character(p$Subject_ID))
  if (anyNA(match_id)) stop("Every assessment must match one participant in participants.")
  parse_date <- function(v, label) {
    s <- as.character(v)
    parsed <- suppressWarnings(as.Date(s, format = "%Y-%m-%d"))
    if (anyNA(parsed) || any(format(parsed, "%Y-%m-%d") != s))
      stop(label, " must contain valid ISO dates: YYYY-MM-DD.")
    parsed
  }
  baseline <- parse_date(p$Baseline_date, "Baseline_date")
  visit <- parse_date(a$Assessment_date, "Assessment_date")
  elapsed <- as.numeric(visit - baseline[match_id]) / 365.25
  age <- suppressWarnings(as.numeric(as.character(p$Age_baseline_years[match_id])))
  x <- data.frame(Subject_ID = as.character(a$Subject_ID), Diagnosis = p$Diagnosis[match_id],
    Timepoint = a$Timepoint, Platform = a$Platform, Gender = p$Gender[match_id],
    FSIQ = p$FSIQ[match_id], Age_baseline_years = p$Age_baseline_years[match_id],
    Time_since_baseline_years = elapsed, stringsAsFactors = FALSE, check.names = FALSE)
  for (nm in cantab_outcome_columns) x[[nm]] <- a[[nm]]
  x$Age_years <- age + elapsed
  tmp <- tempfile(fileext = ".csv"); on.exit(unlink(tmp), add = TRUE)
  write.csv(x, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  read_cantab_input(tmp) # validate before writing to the requested destination
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
  if (!file.copy(tmp, output_file, overwrite = overwrite)) stop("Could not write output file.")
  message("Created ", nrow(x), " observations / ", length(unique(x$Subject_ID)),
          " participants: ", output_file)
  invisible(normalizePath(output_file, winslash = "/", mustWork = TRUE))
}

cantab_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install the R package digest.")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

# Optional reference_dir fixes the included cohort; it never changes raw diagnoses.
# Without it, all supplied participants form the requested analysis cohort.
cantab_prepare_input <- function(input_file, output_dir, reference_dir = NULL, outcome) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install the R package jsonlite.")
  before <- cantab_sha256(input_file)
  dat <- read_cantab_input(input_file)
  keep <- rep(TRUE, nrow(dat))
  refs <- list()
  if (!is.null(reference_dir)) {
    reference_dir <- normalizePath(reference_dir, winslash = "/", mustWork = TRUE)
    ref_file <- file.path(reference_dir, "diagnosis_reference.csv")
    if (!file.exists(ref_file))
      stop("reference_dir must contain diagnosis_reference.csv; use NULL for your full input cohort.")
    ref <- cantab_read_characters(ref_file)
    if (!identical(names(ref), c("Subject_ID", "Diagnosis")) ||
        anyNA(ref) || anyDuplicated(ref$Subject_ID) ||
        any(ref$Subject_ID != trimws(ref$Subject_ID)) ||
        any(grepl("[[:cntrl:]]", ref$Subject_ID)) ||
        any(!ref$Diagnosis %in% c("HC", "BD", "MDD")))
      stop("Invalid diagnosis_reference.csv: unique Subject_ID,Diagnosis columns required.")
    refs$diagnosis_reference <- cantab_sha256(ref_file)
    meta_file <- file.path(reference_dir, "diagnosis_reference.json")
    if (file.exists(meta_file)) {
      meta <- jsonlite::fromJSON(meta_file)
      if (!identical(meta$schema, "cantab-diagnosis-reference-v1") ||
          !identical(meta$key, "Subject_ID") || !identical(meta$column, "Diagnosis") ||
          !identical(meta$uncovered_policy, "exclude_not_in_reference") ||
          !identical(meta$reference_sha256, refs$diagnosis_reference) ||
          !identical(as.integer(meta$reference_participants), nrow(ref)) ||
          !is.character(meta$source_sha256) ||
          !grepl("^[a-f0-9]{64}$", meta$source_sha256))
        stop("Invalid or changed diagnosis reference metadata/hash.")
      refs$metadata <- cantab_sha256(meta_file)
    }
    idx <- match(dat$Subject_ID, ref$Subject_ID)
    keep <- !is.na(idx)
    if (any(dat$Diagnosis[keep] != ref$Diagnosis[idx[keep]]))
      stop("Input disagrees with diagnosis_reference.csv; reconcile the source before fitting.")
    if (!any(keep)) stop("No reference participants occur in input.")
  }
  if (!identical(before, cantab_sha256(input_file))) stop("Input changed during validation.")
  identity <- list(schema = "cantab-portable-R-input-v1", input_sha256 = before,
    reference_sha256 = refs, outcome = outcome)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(output_dir, "analysis_input.csv")
  audit_path <- file.path(output_dir, "INPUT_AUDIT_R.rds")
  if (file.exists(audit_path)) {
    old <- readRDS(audit_path)
    if (!identical(old$identity, identity) || !file.exists(target) ||
        !identical(old$analysis_sha256, cantab_sha256(target)))
      stop("Input/reference/output changed; start a fresh run in a new directory.")
    return(normalizePath(target, winslash = "/", mustWork = TRUE))
  }
  if (file.exists(target)) stop("analysis_input.csv exists without its R input audit; use a new directory.")
  if (all(keep)) {
    if (!file.copy(input_file, target, overwrite = FALSE)) stop("Could not copy input.")
  } else {
    # Preserve nonmissing cell strings, extra columns and source row order.
    # Accepted missing-value markers are written as empty CSV cells.
    write.csv(attr(dat, "raw_characters")[keep, , drop = FALSE], target,
              row.names = FALSE, na = "", fileEncoding = "UTF-8")
  }
  audit <- list(identity = identity, analysis_sha256 = cantab_sha256(target),
    source_rows = nrow(dat), source_participants = length(unique(dat$Subject_ID)),
    analysis_rows = sum(keep), analysis_participants = length(unique(dat$Subject_ID[keep])),
    excluded_rows = sum(!keep), cohort = if (is.null(reference_dir)) "all supplied participants" else
      "participants listed in supplied diagnosis reference", diagnosis_labels_changed = 0L,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  saveRDS(audit, audit_path)
  jsonlite::write_json(audit, file.path(output_dir, "INPUT_AUDIT_R.json"),
                       pretty = TRUE, auto_unbox = TRUE)
  row_index <- rep(NA_integer_, nrow(dat)); row_index[keep] <- seq_len(sum(keep))
  write.csv(data.frame(source_data_row = seq_len(nrow(dat)), Subject_ID = dat$Subject_ID,
    Timepoint = dat$Timepoint, Diagnosis = dat$Diagnosis, included = keep,
    analysis_data_row = row_index), file.path(output_dir, "COHORT_SELECTION.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8")
  normalizePath(target, winslash = "/", mustWork = TRUE)
}
