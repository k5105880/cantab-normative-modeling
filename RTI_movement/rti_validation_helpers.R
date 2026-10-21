# Validation only: these helpers do not change the fitted likelihood or priors.
rti_strict_numeric <- function(x,column='numeric input') {
  missing <- is.na(x) | trimws(as.character(x))==''
  value <- suppressWarnings(as.numeric(x))
  if(any(!missing & !is.finite(value))) stop('Invalid numeric value in column: ',column)
  value[missing] <- NA_real_
  value
}

rti_oof_index <- function(y, response) {
  nm <- names(y)
  if(is.null(nm) || length(y)!=length(response) || anyNA(nm) || any(!grepl('^[0-9]+$',nm)))
    stop('OOF predictions require explicit complete model-row indices.')
  idx <- as.integer(nm)
  if(anyNA(idx) || anyDuplicated(idx) || !identical(sort(idx),seq_along(response)))
    stop('OOF indices are not a unique complete partition.')
  if(!isTRUE(all.equal(as.numeric(y),as.numeric(response[idx]),tolerance=1e-12)))
    stop('OOF response/index alignment failed.')
  idx
}

rti_input_units <- function(dat) {
  aliases <- c(Age='Age_years',Age_baseline='Age_baseline_years',Time_since_baseline='Time_since_baseline_years')
  for(day in names(aliases)) {
    year <- aliases[[day]]
    if(!year %in% names(dat)) next
    y <- suppressWarnings(as.numeric(dat[[year]]))
    if(any(!is.na(dat[[year]]) & !is.finite(y))) stop('Invalid years: ',year)
    if(day %in% names(dat)) {
      v <- suppressWarnings(as.numeric(dat[[day]]))
      if(any(!is.na(dat[[day]]) & !is.finite(v))) stop('Invalid days: ',day)
      ok <- is.finite(v) & is.finite(y)
      if(any(abs(v[ok]/365.25-y[ok])>1e-6)) stop('Conflicting day/year units: ',day)
    } else dat[[day]] <- y*365.25
  }
  # The documented 14-column format supplies baseline age and elapsed years.
  # Derive current age only when neither current-age column was supplied.
  if (!'Age' %in% names(dat) &&
      all(c('Age_baseline', 'Time_since_baseline') %in% names(dat))) {
    dat$Age <- suppressWarnings(as.numeric(as.character(dat$Age_baseline))) +
      suppressWarnings(as.numeric(as.character(dat$Time_since_baseline)))
  }
  dat
}

rti_true_baseline <- function(dat,tol=1e-6) {
  label <- !is.na(dat$Timepoint) & dat$Timepoint==1
  zero <- is.finite(dat$Time_since_baseline_years) & abs(dat$Time_since_baseline_years)<=tol
  if(any(label & !zero)) stop('Timepoint 1 must have time since baseline zero.')
  idx <- which(label | zero)
  if(length(idx)>1) stop('Multiple baseline rows for one participant.')
  if(!length(idx)) NA_integer_ else idx
}

rti_diagnostics_pass <- function(d) {
  cols <- c('n_divergent','n_max_treedepth','max_rhat','min_ess_bulk','min_ess_tail','min_ebfmi')
  if(!all(cols %in% names(d))) return(rep(FALSE,nrow(d)))
  finite <- apply(d[cols],1,function(v) all(is.finite(v)))
  finite & with(d,n_divergent==0 & n_max_treedepth==0 & max_rhat<1.01 &
    min_ess_bulk>=400 & min_ess_tail>=400 & min_ebfmi>=.3)
}

rti_validate_input <- function(dat) {
  key <- c('Subject_ID','Timepoint')
  if(!all(key %in% names(dat))) stop('Missing participant/timepoint keys')
  if(anyNA(dat[key]) || any(!nzchar(trimws(as.character(dat$Subject_ID))))) stop('Missing participant/timepoint key')
  if(anyDuplicated(dat[key])) stop('Duplicate participant/timepoint rows')
  if(any(!as.character(dat$Platform) %in% c('RS6','Connect'))) stop('Unknown or missing platform')
  groups <- unique(as.character(dat$Diagnosis))
  if(anyNA(dat$Diagnosis) || any(!groups %in% c('HC','BD','MDD'))) stop('Unknown or missing diagnosis')
  if(any(is.finite(dat$Time_since_baseline_years) & dat$Time_since_baseline_years < 0)) stop('Negative time since baseline')
  rt <- dat[['RTI_Median_five-choice_movement_time']]
  if(any(!is.na(rt) & (!is.finite(rt) | rt<=0))) stop('RTI movement time must be positive finite milliseconds')
  if(any(!is.finite(dat$Timepoint) | dat$Timepoint<1 | dat$Timepoint!=floor(dat$Timepoint))) stop('Invalid timepoint')
  for(a in split(dat,dat$Subject_ID,drop=TRUE)) rti_true_baseline(a)
  invisible(TRUE)
}

rti_add_bfmi <- function(diag,fit) {
  np <- brms::nuts_params(fit)
  energy <- split(np$Value[np$Parameter=='energy__'],np$Chain[np$Parameter=='energy__'])
  diag$min_ebfmi <- if(length(energy)) min(vapply(energy,function(e)mean(diff(e)^2)/var(e),numeric(1))) else NA_real_
  diag$validation_pass <- rti_diagnostics_pass(diag)
  diag
}

rti_validate_kfold <- function(kf,model_dat,role,seed,diagnose,max_treedepth) {
  if(is.null(kf$fits)) stop('Fold fits required for validation')
  omitted_all <- integer(); rows <- list()
  for(k in seq_len(nrow(kf$fits))) {
    f <- kf$fits[[k,'fit']]; idx <- as.integer(kf$fits[[k,'omitted']])
    omitted_all <- c(omitted_all,idx)
    if(any(as.character(model_dat$Subject_ID[idx]) %in% as.character(f$data$Subject_ID))) stop('Subject leakage in CV fold ',k)
    d <- diagnose(f,role,as.character(f$family$family),paste0('fold_',k))
    d <- rti_add_bfmi(d,f);d$fold <- k;d$seed <- seed
    rows[[k]] <- d
  }
  if(!identical(sort(omitted_all),seq_len(nrow(model_dat)))) stop('Invalid CV partition')
  d <- do.call(rbind,rows)
  name <- paste0('cv_diagnostics_',gsub('[^A-Za-z0-9_-]','_',role),'_',seed,'.csv')
  write.csv(d,name,row.names=FALSE,na='')
  if(!all(d$validation_pass)) stop('CV convergence gate failed; inspect ',name)
  invisible(d)
}

