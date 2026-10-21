arg<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE);root<-dirname(dirname(normalizePath(sub('^--file=','',arg[grepl('^--file=',arg)][1]))))
stopifnot(file.exists(file.path(root,'R/config.R')))
source(file.path(root,'R/config.R'))
cfg<-jsonlite::fromJSON(file.path(root,'config.json'))
validate_config(cfg)
expect_config_error <- function(config, fields) {
  error <- tryCatch({ validate_config(config); NULL }, error = identity)
  if (!inherits(error, "error"))
    stop("Invalid configuration was accepted: ", paste(fields, collapse = ", "))
  for (field in fields) {
    if (!grepl(field, conditionMessage(error), fixed = TRUE))
      stop("Configuration error did not identify '", field, "': ", conditionMessage(error))
  }
}

# Misspelled knobs must stop before sampling or scoring can ignore them.
for (case in list(c("sampling", "thin"), c("scoring", "profile"), c("prior", "draw"),
                  c("sampling_recovery", "retry"))) {
  bad <- cfg
  bad[[case[1]]][[case[2]]] <- 5
  expect_config_error(bad, case)
}

# Every documented field is required, except the entire recovery policy.
objects <- list(config = cfg[names(cfg) != "sampling_recovery"],
                sampling = cfg$sampling, scoring = cfg$scoring, prior = cfg$prior,
                sampling_recovery = cfg$sampling_recovery)
for (object in names(objects)) {
  for (field in names(objects[[object]])) {
    bad <- cfg
    if (object == "config") bad[[field]] <- NULL else bad[[object]][[field]] <- NULL
    expect_config_error(bad, c(if (object != "config") object, field))
    bad <- cfg
    if (object == "config") bad[field] <- list(NULL) else bad[[object]][field] <- list(NULL)
    expect_config_error(bad, c(if (object != "config") object, field))
  }
}

# JSON objects must retain unambiguous field names at every level.
for (object in names(objects)) {
  bad <- cfg
  if (object == "config") {
    bad <- c(bad, bad["seed"])
  } else {
    bad[[object]] <- c(bad[[object]], bad[[object]][1L])
  }
  expect_config_error(bad, c(object, "duplicate"))
  for (invalid_name in c("", NA_character_)) {
    bad <- cfg
    if (object == "config") names(bad)[1L] <- invalid_name else
      names(bad[[object]])[1L] <- invalid_name
    expect_config_error(bad, object)
  }
}
duplicate_json <- jsonlite::fromJSON('{"seed":20260922,"seed":1}')
bad <- c(cfg[names(cfg) != "seed"], duplicate_json)
expect_config_error(bad, c("seed", "duplicate"))

numeric_fields <- list(
  config = c("schema_version", "seed", "calibration_bootstraps"),
  sampling = names(cfg$sampling),
  scoring = c("tolerance", "mcse_limit", "workers", "chunk"),
  prior = names(cfg$prior),
  sampling_recovery = c("max_retries", "warmup", "adapt_delta")
)
for (object in names(numeric_fields)) {
  for (field in numeric_fields[[object]]) {
    original <- if (object == "config") cfg[[field]] else cfg[[object]][[field]]
    for (value in list(as.character(original), TRUE, NA_real_, NaN, Inf, -Inf,
                       numeric(), rep(original, 2L), list(original),
                       matrix(original, nrow = 1L))) {
      bad <- cfg
      if (object == "config") bad[[field]] <- value else bad[[object]][[field]] <- value
      expect_config_error(bad, c(if (object != "config") object, field))
    }
  }
}
for (value in list(c("label", "other"), TRUE, NA_character_, character(), 1, list("label"), "")) {
  bad <- cfg; bad$analysis <- value
  expect_config_error(bad, "analysis")
}
for (value in list(c("divergence_only", "divergence_only"), TRUE, NA_character_,
                   character(), 1, list("divergence_only"))) {
  bad <- cfg; bad$sampling_recovery$trigger <- value
  expect_config_error(bad, c("sampling_recovery", "trigger"))
}
for (value in list(as.character(cfg$scoring$nodes), list(32, 64, 128, 256),
                   c(32, 64, NA, 256), matrix(cfg$scoring$nodes, nrow = 1L))) {
  bad <- cfg; bad$scoring$nodes <- value
  expect_config_error(bad, c("scoring", "nodes"))
}
for (object in names(objects)) {
  for (value in list(1, "object", TRUE, list(), c(named = 1), data.frame(field = 1))) {
    bad <- cfg
    if (object == "config") bad <- value else bad[[object]] <- value
    expect_config_error(bad, object)
  }
}
cat("PASS: all configuration objects reject unknown, missing, duplicate, null and malformed fields with field-specific errors.\n")

# JSON integers and doubles are both valid; labels remain descriptive.
double_cfg <- cfg
for (object in names(numeric_fields)) {
  for (field in numeric_fields[[object]]) {
    if (object == "config") double_cfg[[field]] <- as.double(cfg[[field]]) else
      double_cfg[[object]][[field]] <- as.double(cfg[[object]][[field]])
  }
}
double_cfg$scoring$nodes <- as.double(cfg$scoring$nodes)
double_cfg$analysis <- "independent reproducibility run"
validate_config(double_cfg)
integer_nodes <- cfg; integer_nodes$scoring$nodes <- as.integer(cfg$scoring$nodes)
validate_config(integer_nodes)
explicit_null <- jsonlite::fromJSON('{"sampling_recovery":null}')
without_recovery <- cfg
without_recovery["sampling_recovery"] <- explicit_null
validate_config(without_recovery)
cat("PASS: valid numeric JSON types, descriptive labels and explicit null recovery remain supported.\n")

# The schema rejects values outside the established numerical contract.
invalid_ranges <- list(
  list(NULL, "schema_version", c(0, 1.5, 2)),
  list(NULL, "seed", c(0, 1.5, 2e9)),
  list(NULL, "calibration_bootstraps", c(1999, 2000.5, .Machine$integer.max + 1)),
  list("sampling", "chains", c(0, 3, 5)),
  list("sampling", "iter", c(2499, 2500.5, .Machine$integer.max + 1)),
  list("sampling", "warmup", c(0, 1500.5, 4001)),
  list("sampling", "cores", c(0, 1.5, 5)),
  list("sampling", "adapt_delta", c(.989, 1)),
  list("sampling", "max_treedepth", c(14, 15.5)),
  list("sampling", "refresh", c(0, 1.5)),
  list("scoring", "tolerance", c(0, -.001, .0021)),
  list("scoring", "mcse_limit", c(0, -.01, .0301)),
  list("scoring", "chunk", c(0, 1.5)),
  list("prior", "draws", c(39999, 40000.5, 40001)),
  list("prior", "batch_size", c(0, 1.5, 501)),
  list("sampling_recovery", "max_retries", c(0, 1.5, 2)),
  list("sampling_recovery", "warmup", c(0, 1500, 2500.5, .Machine$integer.max)),
  list("sampling_recovery", "adapt_delta", c(.9995, 1))
)
for (case in invalid_ranges) {
  for (value in case[[3]]) {
    bad <- cfg
    if (is.null(case[[1]])) bad[[case[[2]]]] <- value else
      bad[[case[[1]]]][[case[[2]]]] <- value
    # Cross-field constraints may identify either member of the relationship.
    error <- tryCatch({ validate_config(bad); NULL }, error = identity)
    stopifnot(inherits(error, "error"),
              grepl(if (is.null(case[[1]])) case[[2]] else case[[1]],
                    conditionMessage(error), fixed = TRUE))
  }
}
boundary_cfg <- cfg
boundary_cfg$seed <- 1L
boundary_cfg$sampling <- list(chains = 4L, iter = 1001L, warmup = 1L, cores = 1L,
                             adapt_delta = .99, max_treedepth = 15L, refresh = 1L)
boundary_cfg$sampling_recovery <- list(trigger = "divergence_only", max_retries = 1L,
                                      warmup = 2L, adapt_delta = .999)
boundary_cfg$scoring$tolerance <- .001
boundary_cfg$scoring$mcse_limit <- .01
boundary_cfg$scoring$workers <- 1L
boundary_cfg$scoring$chunk <- 1L
boundary_cfg$prior$batch_size <- 1L
validate_config(boundary_cfg)
boundary_cfg$sampling_recovery <- NULL
boundary_cfg$seed <- 1999999999L
boundary_cfg$sampling$iter <- .Machine$integer.max
boundary_cfg$calibration_bootstraps <- .Machine$integer.max
boundary_cfg$prior$batch_size <- 40000L
validate_config(boundary_cfg)
cat("PASS: established numerical bounds, retained-draw requirements and valid boundary settings are preserved.\n")

for(nodes in list(c(2,2),c(32,32),c(64,32,128,256),c(32,64.5,128,256))) {
  bad<-cfg;bad$scoring$nodes<-nodes
  stopifnot(inherits(try(validate_config(bad),silent=TRUE),'try-error'))
}
for(value in c(0,-1,1.5,NA,Inf)) {
  bad<-cfg;bad$scoring$workers<-value
  stopifnot(inherits(try(validate_config(bad),silent=TRUE),'try-error'))
}
cat('PASS: fixed quadrature rules and integer resources enforced.\n')
for(change in list(list(max_retries=2L),list(trigger='any_failure'),list(warmup=1500L),
 list(warmup=2500.5),list(adapt_delta=.9995),list(adapt_delta=1))) {
 bad<-cfg;bad$sampling_recovery<-modifyList(bad$sampling_recovery,change)
 stopifnot(inherits(try(validate_config(bad),silent=TRUE),'try-error'))
}
without_recovery<-cfg;without_recovery$sampling_recovery<-NULL;validate_config(without_recovery)
cat('PASS: at most one explicitly configured divergence-only retry; omitting the recovery policy disables retries.\n')
for(field in c('score_profile','delta','unexpected')) {
  bad<-cfg;bad[[field]]<-'all'
  stopifnot(inherits(try(validate_config(bad),silent=TRUE),'try-error'))
}
cat('PASS: unsupported configuration fields cannot enable other analyses.\n')
local({
  old <- Sys.getenv("CANTAB_JOB_CORES", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("CANTAB_JOB_CORES") else
    Sys.setenv(CANTAB_JOB_CORES = old), add = TRUE)
  Sys.setenv(CANTAB_JOB_CORES = "4")
  capped <- cfg
  capped$sampling$cores <- 2L
  capped$scoring$workers <- 3L
  validate_config(capped)
  adjusted <- apply_runtime_budget(capped)
  stopifnot(adjusted$sampling$cores == 2L, adjusted$scoring$workers == 3L,
            adjusted$sampling$chains == cfg$sampling$chains,
            adjusted$sampling$iter == cfg$sampling$iter)
  Sys.setenv(CANTAB_JOB_CORES = "1")
  serial <- apply_runtime_budget(capped)
  stopifnot(serial$sampling$cores == 1L, serial$scoring$workers == 1L,
            serial$sampling$chains == cfg$sampling$chains)
  cat("PASS: allocation respects configuration caps and the public core budget without changing sampling effort.\n")
})
