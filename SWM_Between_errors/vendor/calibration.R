# Calibration, platform comparisons and prior predictive summaries.
bootstrap_platform_gap <- function(d, person_equal = FALSE, B = 2000, seed = 1) {
    if (!all(c("RS6", "Connect") %in% as.character(d$Platform))) 
        return(NULL)
    a <- aggregate(d$pit_mean, list(id = as.character(d$Subject_ID), platform = as.character(d$Platform)), function(x) c(sum = sum(x), 
        n = length(x)))
    ss <- a$x[, "sum"]
    nn <- a$x[, "n"]
    if (person_equal) {
        ss <- ss/nn
        nn[] <- 1
    }
    ids <- unique(a$id)
    sR <- sC <- nR <- nC <- numeric(length(ids))
    for (i in seq_len(nrow(a))) {
        j <- match(a$id[i], ids)
        if (a$platform[i] == "RS6") {
            sR[j] <- ss[i]
            nR[j] <- nn[i]
        }
        else {
            sC[j] <- ss[i]
            nC[j] <- nn[i]
        }
    }
    effect <- function(ix) if (sum(nC[ix]) > 0 && sum(nR[ix]) > 0) 
        sum(sC[ix])/sum(nC[ix]) - sum(sR[ix])/sum(nR[ix])
    else NA_real_
    set.seed(seed)
    boot <- replicate(B, effect(sample.int(length(ids), length(ids), replace = TRUE)))
    ci <- quantile(boot, c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE, names = FALSE)
    out <- data.frame(weighting = if (person_equal) 
        "person_equal_within_platform"
    else "all_visits", n_RS6 = sum(d$Platform == "RS6"), n_Connect = sum(d$Platform == "Connect"), subjects_RS6 = sum(nR > 
        0), subjects_Connect = sum(nC > 0), gap_Connect_minus_RS6 = effect(seq_along(ids)), lo95 = ci[1], lo90 = ci[2], 
        hi90 = ci[3], hi95 = ci[4], bootstrap_finite = sum(is.finite(boot)))
    for (margin in c(0.1, 0.2, 0.3)) out[[paste0("within_illustrative_margin_", margin)]] <- ci[2] > -margin && 
        ci[3] < margin
    out
}

summarize_calibration <- function(d, person_equal = FALSE) {
    w <- rep(1, nrow(d))
    if (person_equal) 
        w <- 1/as.numeric(table(as.character(d$Subject_ID))[as.character(d$Subject_ID)])
    w <- w/sum(w)
    av <- function(x) sum(w * x)
    z <- data.frame(n = nrow(d), subjects = length(unique(d$Subject_ID)), weighting = if (person_equal) 
        "person_equal"
    else "all_visits", PIT_mean = av(d$pit_mean), PIT_SD = sqrt(max(0, av(d$pit_second) - av(d$pit_mean)^2)), PIT_low_1.96 = av(d$pit_low), 
        PIT_high_1.96 = av(d$pit_high), observed_zero = av(d$NM_response == 0), predicted_zero = av(d$pred_zero), 
        mean_log_density = av(d$log_pred_density))
    for (level in c(50, 80, 95)) {
        z[[paste0("observed_coverage_", level)]] <- av(d[[paste0("coverage_", level)]])
        z[[paste0("predictive_coverage_", level)]] <- av(d[[paste0("expected_coverage_", level)]])
    }
    z
}

summarize_prior_predictions <- function(fit, d, seed) {
    set.seed(seed)
    nd <- d
    nd$Subject_ID <- factor(paste0("prior_new_", as.character(d$Subject_ID)))
    prep <- prepare_predictions(fit, newdata = nd, re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian", 
        ndraws = min(1200, posterior::ndraws(fit)))
    pp <- posterior_predict(prep)
    bind_rows(lapply(c("RS6", "Connect"), function(pl) {
        values <- as.numeric(pp[, d$Platform == pl, drop = FALSE])
        qq <- quantile(values[is.finite(values)], c(0.5, 0.95, 0.99), names = FALSE)
        data.frame(Platform = pl, zero_fraction = mean(values == 0, na.rm = TRUE), q50 = qq[1], q95 = qq[2], q99 = qq[3], 
            fraction_gt1000 = mean(values > 1000, na.rm = TRUE), nonfinite_fraction = mean(!is.finite(values)))
    }))
}

