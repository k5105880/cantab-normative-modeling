# Main model, held-out scoring and diagnostic outputs.
calibrate_normative_scores <- function(dat,spec,selected_family,oof) {
  oof_col <- 'normative_z_cross_validation';rpit_col <- 'PIT_cross_validation';rpit_z_col <- 'PIT_z_cross_validation'
  oof_hc <- dat %>%
    filter(Diagnosis == control_label, !is.na(.data[[oof_col]]))

  moments <- sample_moments(oof_hc[[oof_col]])
  oof_calibration_overall <- data.frame(
    outcome = spec$outcome_prefix,
    family_id = selected_family,
    calibration_method = "predictive_CDF",
    n = nrow(oof_hc),
    n_subjects = n_distinct(oof_hc$Subject_ID),
    mean_z = unname(moments["mean"]),
    sd_z = unname(moments["sd"]),
    skewness_z = unname(moments["skewness"]),
    excess_kurtosis_z = unname(moments["excess_kurtosis"]),
    stringsAsFactors = FALSE
  )

  rpit_moments <- sample_moments(oof_hc[[rpit_z_col]])
  oof_rpit_calibration_overall <- data.frame(
    outcome = spec$outcome_prefix,
    family_id = selected_family,
    calibration_method = "randomized_PIT",
    n = sum(is.finite(oof_hc[[rpit_col]])),
    n_subjects = n_distinct(oof_hc$Subject_ID[is.finite(oof_hc[[rpit_col]])]),
    mean_pit = mean(oof_hc[[rpit_col]], na.rm = TRUE),
    sd_pit = sd(oof_hc[[rpit_col]], na.rm = TRUE),
    target_mean_pit = 0.5,
    target_sd_pit = 1 / sqrt(12),
    mean_z = unname(rpit_moments["mean"]),
    sd_z = unname(rpit_moments["sd"]),
    skewness_z = unname(rpit_moments["skewness"]),
    excess_kurtosis_z = unname(rpit_moments["excess_kurtosis"]),
    stringsAsFactors = FALSE
  )

  calibration_z_col <- oof_col
  calibration_method <- "predictive_CDF_z"

  oof_calibration_by_platform <- oof_hc %>%
    group_by(Platform) %>%
    summarise(
      outcome = spec$outcome_prefix,
      family_id = selected_family,
      calibration_method = calibration_method,
      n = n(),
      n_subjects = n_distinct(Subject_ID),
      mean_z = mean(.data[[calibration_z_col]], na.rm = TRUE),
      sd_z = sd(.data[[calibration_z_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(outcome, family_id, calibration_method)

  oof_calibration_by_practice <- oof_hc %>%
    group_by(Practice_bin) %>%
    summarise(
      outcome = spec$outcome_prefix,
      family_id = selected_family,
      calibration_method = calibration_method,
      n = n(),
      n_subjects = n_distinct(Subject_ID),
      mean_z = mean(.data[[calibration_z_col]], na.rm = TRUE),
      sd_z = sd(.data[[calibration_z_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(outcome, family_id, calibration_method)

  oof_calibration_by_covariate <- make_calibration_by_covariate(
    oof_hc = oof_hc,
    calibration_z_col = calibration_z_col,
    outcome = spec$outcome_prefix,
    family_id = selected_family,
    calibration_method = calibration_method
  )

  oof_coverage <- coverage_from_yrep(oof$yrep, oof$y) %>%
    mutate(
      outcome = spec$outcome_prefix,
      family_id = selected_family,
      .before = 1
    )
  list(OOF_calibration=oof_calibration_overall,OOF_PIT_calibration=oof_rpit_calibration_overall,
    OOF_calibration_by_platform=oof_calibration_by_platform,OOF_calibration_by_practice=oof_calibration_by_practice,
    OOF_calibration_by_covariate=oof_calibration_by_covariate,OOF_interval_coverage=oof_coverage)
}
run_normative_pipeline <- function() {
  spec <- outcome_specs[[1L]];family_id <- spec$primary_family
  model_dat <- prepare_model_data(data,spec,family_id,control_only=TRUE)
  fold_obj <- make_grouped_folds(model_dat,K=cv_K,seed=seed_main+100L,
    balanced=cv_balance_folds,n_candidates=cv_balance_candidates)
  model_dat$CV_fold <- fold_obj$folds
  checkpoint <- file.path(checkpoint_dir,'normative_model.rds')
  if(file.exists(checkpoint)) {
    selected <- readRDS(checkpoint)
    if(!identical(selected$model_dat,model_dat)||!identical(selected$folds,fold_obj$folds))
      stop('Model checkpoint training rows or folds changed')
  } else {
    full <- fit_normative_model(model_dat,spec$outcome_prefix,family_id,seed=seed_main+1000L+1L)
    kf <- run_grouped_kfold(full$fit,fold_obj$folds,save_fits=TRUE,seed=seed_main+10000L+1L)
    selected <- list(fit=full$fit,definition=full$definition,model_dat=model_dat,kfold=kf,folds=fold_obj$folds)
    engine_atomic_save(selected,checkpoint)
  }
  full_diag <- rvp_add_bfmi(diagnose_brms_fit(selected$fit,spec$outcome_prefix,family_id,selected$definition$model_label),selected$fit)
  if(!all(full_diag$validation_pass))stop('Model convergence gate failed')
  rvp_validate_kfold(selected$kfold,model_dat,spec$outcome_prefix,seed_main+10000L+1L,diagnose_brms_fit,max_treedepth_fit)
  cv_summary <- extract_kfold_summary(selected$kfold,spec$outcome_prefix,family_id,prespecified=TRUE)
  if(nrow(cv_summary)!=1L||!is.finite(cv_summary$elpd_kfold))stop('Cross-validation summary is incomplete')
  cv_summary$status <- 'ok'
  pred_dat <- prepare_model_data(data,spec,family_id,control_only=FALSE)
  full_z <- score_reference_population(selected$fit,pred_dat,family_id)
  oof <- score_cross_validation(selected$kfold,family_id,model_dat,rpit_seed=seed_main+31000L)
  dat <- data
  for(column in c('normative_z_full_reference','normative_z_cross_validation','PIT_cross_validation','PIT_z_cross_validation'))dat[[column]] <- NA_real_
  dat$normative_z_full_reference[pred_dat$row_id_global] <- full_z
  dat$normative_z_cross_validation[oof$row_id_global] <- oof$z
  dat$PIT_cross_validation[oof$row_id_global] <- oof$pit_randomized
  dat$PIT_z_cross_validation[oof$row_id_global] <- oof$z_rpit
  hc <- as.character(dat$Diagnosis)==control_label
  dat$normative_z <- dat$normative_z_full_reference
  dat$normative_z[hc] <- dat$normative_z_cross_validation[hc]
  dat$model_eligible <- dat$row_id_global%in%pred_dat$row_id_global
  if(any(dat$model_eligible&!is.finite(dat$normative_z)))stop('An eligible row lacks an accepted normative score')
  dat$precision_status <- ifelse(dat$model_eligible,'ok','ineligible_missing_or_invalid_input')
  dat$score_source <- ifelse(hc,'HC_OOF','full_HC_new_subject')
  calibration <- calibrate_normative_scores(dat,spec,family_id,oof)
  prior <- tryCatch(run_prior_predictive_check(model_dat,selected$definition,spec,family_id,outcome_index=1L),
    error=function(e)data.frame(status='failed',reason=conditionMessage(e)))
  write.csv(prior,'prior_predictive_check.csv',row.names=FALSE,na='')
  if(nrow(prior)!=1L||!identical(prior$status,'ok'))warning('Prior predictive diagnostic is incomplete; inspect prior_predictive_check.csv')
  for(name in names(calibration))write.csv(calibration[[name]],paste0(name,'.csv'),row.names=FALSE,na='')
  write.csv(full_diag,'full_fit_diagnostics.csv',row.names=FALSE,na='')
  write.csv(cv_summary,'cross_validation_summary.csv',row.names=FALSE,na='')
  write.csv(model_dat[c('row_id_global','Subject_ID','CV_fold')],'subject_folds.csv',row.names=FALSE)
  write.csv(fold_obj$audit,'fold_balance_audit.csv',row.names=FALSE)
  # Keep source cells exact in the CSV; model and calibration calculations above remain numeric.
  export_rows <- dat
  raw_cells <- read.csv(infile,check.names=FALSE,colClasses='character',na.strings=NULL,
    fileEncoding='UTF-8-BOM',stringsAsFactors=FALSE)
  if(nrow(raw_cells)!=nrow(dat))stop('Source row count changed before export')
  generated_columns <- c('row_id_global','normative_z','normative_z_full_reference',
    'normative_z_cross_validation','PIT_cross_validation','PIT_z_cross_validation',
    'model_eligible','precision_status','score_source')
  for(column in setdiff(names(raw_cells),generated_columns))export_rows[[column]] <- raw_cells[[column]]
  write.csv(export_rows,'RVP_latency_scores.csv',row.names=FALSE,na='')
  engine_atomic_save(list(spec=spec,selected_fit=selected$fit,selected_definition=selected$definition,
    selected_model_dat=model_dat,selected_kfold=selected$kfold,folds=fold_obj$folds,oof=oof,
    prior_predictive_check=prior,state_cdf_settings=state_cdf_settings),'normative_objects.rds')
  writeLines('COMPLETE','model_run_status.txt')
  invisible(dat)
}
