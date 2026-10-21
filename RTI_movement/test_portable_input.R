# Synthetic regression for the four normative engines; no study data or MCMC.
# Run from R with source(), or with Rscript. Requires dplyr, loo, digest, jsonlite.
local({
  paths <- as.character(unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE))
  paths <- paths[basename(paths) == 'test_portable_input.R']
  if (!length(paths)) paths <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))
  file <- gsub('~+~', ' ', tail(paths, 1L), fixed = TRUE)
  root <- dirname(dirname(normalizePath(file, mustWork = TRUE)))
  suppressPackageStartupMessages(library(dplyr))
  helper <- new.env(parent = globalenv())
  source(file.path(root, 'R/input_data.R'), local = helper)
  work <- tempfile('cantab_identifier_test_'); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  ids <- c('001', '1', '9007199254740992', '9007199254740993', '02', '2')
  dat <- data.frame(Subject_ID = rep(ids, each = 2L), Diagnosis = 'HC',
    Timepoint = rep(1:2, length(ids)), Platform = rep(c('RS6', 'Connect'), length(ids)),
    Gender = rep(rep(c('F', 'M'), 3L), each = 2L), FSIQ = rep(98:103, each = 2L),
    Age_baseline_years = rep(seq(25, 60, length.out = length(ids)), each = 2L),
    Time_since_baseline_years = rep(c(0, 1.5), length(ids)),
    SWM_Between_errors = rep(c(12, 3), length(ids)),
    SWM_Strategy = rep(c(32, 8), length(ids)), RVP_A_prime = .95,
    RVP_Median_latency = 450, check.names = FALSE)
  dat[['RTI_Median_five-choice_reaction_time']] <- 400
  dat[['RTI_Median_five-choice_movement_time']] <- 360
  metadata_names <- c('SWM_Problem_reached', 'SWM Problem reached', 'SWM_Problem_Reached',
    'SWM_ProblemReached', 'SWM Problem Reached', 'SWM_Problem_reached_QC')
  metadata_values <- rep(c('0', '15', '16', '99', NA, 'not recorded'), each = 2L)
  stopifnot(identical(names(dat), helper$cantab_required_columns), ncol(dat) == 14L)
  input <- file.path(work, 'input with BOM.csv')
  write.csv(dat, input, row.names = FALSE, na = '')
  bytes <- readBin(input, 'raw', n = file.info(input)$size)
  writeBin(c(as.raw(c(239, 187, 191)), bytes), input)
  prepared <- helper$cantab_prepare_input(input, file.path(work, 'prepared'), outcome = 'synthetic')
  stopifnot(identical(helper$read_cantab_input(prepared)$Subject_ID, dat$Subject_ID))

  assignment <- function(expr, name, rhs = NULL) {
    is.call(expr) && identical(expr[[1]], as.name('<-')) && identical(expr[[2]], as.name(name)) &&
      (is.null(rhs) || (is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name(rhs))))
  }
  find_assignments <- function(expr, name, rhs = NULL) {
    if (assignment(expr, name, rhs)) return(list(expr))
    if (!is.call(expr) && !is.expression(expr) && !is.list(expr)) return(list())
    unlist(lapply(as.list(expr), find_assignments, name = name, rhs = rhs), recursive = FALSE)
  }
  fold_functions <- c('make_subject_cv_profile', 'score_grouped_fold_balance',
                      'make_fold_balance_audit', 'make_grouped_folds')
  verify_rows_and_folds <- function(x, e, label) {
    stopifnot(identical(as.character(x$Subject_ID), dat$Subject_ID),
      length(unique(x$Subject_ID)) == length(ids), !anyDuplicated(x[c('Subject_ID', 'Timepoint')]))
    e$cv_balance_include_fsiq <- FALSE
    arguments<-list(model_dat=x,K=5L,seed=1334L,n_candidates=1L)
    if('balanced'%in%names(formals(e$make_grouped_folds)))arguments$balanced<-FALSE
    folds<-do.call(e$make_grouped_folds,arguments)$folds
    map <- unique(data.frame(Subject_ID = as.character(x$Subject_ID), fold = folds))
    stopifnot(nrow(map) == length(ids), !anyDuplicated(map$Subject_ID),
      identical(sort(unique(folds)), 1:5))
    cat(label, 'identifiers, visit keys and grouped folds PASS\n')
  }
  for(outcome in c('RTI_movement','RVP_latency')) {
    e<-new.env(parent=globalenv());source(file.path(root,outcome,'core.R'),local=e)
    ctx<-e$normative_context(prepared,file.path(work,outcome))
    verify_rows_and_folds(ctx$data,ctx,outcome)
    stopifnot(isTRUE(all.equal(ctx$data$Age_years,dat$Age_baseline_years+dat$Time_since_baseline_years)))
    baseline_data<-ctx$data
    for(column in metadata_names) {
      extra<-dat;extra[[column]]<-metadata_values
      extra_file<-file.path(work,'input_extra.csv');write.csv(extra,extra_file,row.names=FALSE,na='')
      observed<-e$normative_context(extra_file,file.path(work,outcome))$data
      stopifnot(identical(observed[names(baseline_data)],baseline_data))
    }
    prefix<-if(outcome=='RTI_movement')'rti'else'rvp'
    units<-get(paste0(prefix,'_input_units'),envir=ctx)
    explicit<-dat;explicit$Age_years<-dat$Age_baseline_years+.25
    stopifnot(identical(units(explicit)$Age,explicit$Age_years*365.25))
    explicit<-dat;explicit$Age<-rep(12345,nrow(dat));stopifnot(identical(units(explicit)$Age,explicit$Age))
    dir.create(file.path(work,outcome));write.csv(ctx$data,file.path(work,outcome,paste0(outcome,'_scores.csv')),row.names=FALSE,na='')
    a<-new.env(parent=globalenv());a$input<-prepared;a$run<-file.path(work,outcome)
    audit<-parse(file.path(root,outcome,'audit_run.R'))
    for(name in c('raw','d')) {
      selected<-find_assignments(audit,name,'read.csv');stopifnot(length(selected)==1L);eval(selected[[1]],a)
    }
    stopifnot(identical(a$raw$Subject_ID,dat$Subject_ID),identical(a$d$Subject_ID,dat$Subject_ID))
    cat(outcome,'14-column derivation, unrelated metadata isolation and audit identifier parity PASS\n')
  }

  e <- new.env(parent = globalenv())
  source(file.path(root, 'SWM_Between_errors/R/data.R'), local = e)
  source(file.path(root, 'SWM_Between_errors/vendor/folds.R'), local = e)
  expressions <- parse(file.path(root, 'SWM_Between_errors/R/pipeline.R'))
  reads <- find_assignments(expressions, 'raw', 'read.csv'); stopifnot(length(reads) == 1L)
  e$ctx <- list(input = prepared); eval(reads[[1]], e)
  cfg <- list(groups = c('HC', 'MDD', 'BD'), trials = c(RS6 = 4, Connect = 1), control_label = 'HC')
  x <- e$prepare_raw_data(e$raw, cfg)
  verify_rows_and_folds(x$hc, e, 'SWM_Between_errors')

  e <- new.env(parent = globalenv())
  for (name in c('swm_discrete_family.R', 'swm_core.R', 'swm_folds.R', 'swm_pipeline.R'))
    source(file.path(root, 'SWM_Strategy/R', name), local = e)
  expressions <- parse(file.path(root, 'SWM_Strategy/R/swm_pipeline.R'))
  reads <- find_assignments(expressions, 'raw', 'read.csv'); stopifnot(length(reads) == 1L)
  e$input <- prepared; eval(reads[[1]], e)
  x <- e$swm_adapt_raw(e$raw)
  verify_rows_and_folds(x, e, 'SWM_Strategy')
  stopifnot(all(x$eligible), all(x$score == 8))
  for (column in metadata_names) {
    extra <- e$raw; extra[[column]] <- metadata_values
    adapted <- e$swm_adapt_raw(extra)
    stopifnot("Unrelated export fields must not exclude valid Strategy rows" =
      identical(adapted[names(x)], x))
  }
  invalid <- e$raw
  invalid$SWM_Strategy[1:5] <- c(0, 15, 57, NA, 32.5)
  invalid[[metadata_names[1]]] <- 16
  adapted <- e$swm_adapt_raw(invalid)
  stopifnot(!any(adapted$eligible[1:5]), all(adapted$eligible[-(1:5)]),
    identical(adapted$SWM_Strategy_QC_reason[1:5],
      c('nonpositive_placeholder', 'outside_harmonized_2_14',
        'outside_harmonized_2_14', 'missing_raw', 'noninteger_raw')))
  cat('PASS: exact numeric-looking identifiers and minimal 14-column input across all four engines; no MCMC.\n')
})
