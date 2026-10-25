# Normative score prediction and predictive-interval coverage.
pit_moments <- function(lower, mass, upper) {
    stopifnot(length(lower) == length(mass), length(mass) == length(upper), all(is.finite(c(lower, mass, upper))), 
        all(lower >= 0), all(mass > 0), all(upper >= 0), all(abs(lower + mass + upper - 1) < 1e-08))
    a <- b <- numeric(length(lower))
    lo <- lower <= 0.5
    up <- upper <= 0.5
    a[lo] <- qnorm(lower[lo])
    a[!lo] <- qnorm(pmin(1, upper[!lo] + mass[!lo]), lower.tail = FALSE)
    b[up] <- qnorm(upper[up], lower.tail = FALSE)
    b[!up] <- qnorm(pmin(1, lower[!up] + mass[!up]))
    pa <- dnorm(a)
    pb <- dnorm(b)
    ap <- a * pa
    bp <- b * pb
    ap[!is.finite(a)] <- 0
    bp[!is.finite(b)] <- 0
    m <- (pa - pb)/mass
    e2 <- 1 + (ap - bp)/mass
    small <- mass < 1e-05 & pmin(lower, upper) > 10 * mass
    if (any(small)) {
        nodes <- (c(-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664) + 1)/2
        weights <- c(0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189)/2
        for (i in which(small)) {
            if (mass[i] < 1e-09 * min(lower[i], upper[i])) {
                if (lower[i] <= 0.5) {
                  anchor <- qnorm(lower[i])
                  delta <- mass[i] * nodes
                }
                else {
                  anchor <- qnorm(upper[i], lower.tail = FALSE)
                  delta <- -mass[i] * (1 - nodes)
                }
                step <- delta/dnorm(anchor)
                q <- anchor + step + anchor * step^2/2
            }
            else q <- if (lower[i] <= 0.5) 
                qnorm(lower[i] + mass[i] * nodes)
            else qnorm(upper[i] + mass[i] * (1 - nodes), lower.tail = FALSE)
            m[i] <- sum(weights * q)
            e2[i] <- sum(weights * q * q)
        }
    }
    bad <- !is.finite(m) | !is.finite(e2) | e2 < m * m - 1e-06
    if (any(bad)) {
        mid <- ifelse(lower <= 0.5, qnorm(lower + mass/2), qnorm(upper + mass/2, lower.tail = FALSE))
        m[bad] <- mid[bad]
        e2[bad] <- mid[bad]^2
    }
    data.frame(mean = m, second = pmax(e2, m * m), low = pmin(1, pmax(0, (pnorm(-1.96) - lower)/mass)), high = pmin(1, 
        pmax(0, (pnorm(1.96, lower.tail = FALSE) - upper)/mass)))
}

new_subject_ids <- function(fit, ids, prefix) {
    known <- union(as.character(fit$data$Subject_ID), levels(fit$data$Subject_ID))
    while (any(paste0(prefix, ids) %in% known)) prefix <- paste0(prefix, "_")
    factor(paste0(prefix, ids))
}

predict_normative_coverage <- function(fit, data, seed, ndraws = 8000, coverage = FALSE) {
    set.seed(seed)
    nd <- data
    nd$Subject_ID <- new_subject_ids(fit, as.character(nd$Subject_ID), "new_subject_")
    out <- list()
    gd <- getFromNamespace("get_dpar", "brms")
    for (ix in split(seq_len(nrow(nd)), ceiling(seq_len(nrow(nd))/120))) {
        model_ix <- if (length(ix) == 1L) 
            rep(ix, 2L)
        else ix
        keep <- seq_along(ix)
        prep <- prepare_predictions(fit, newdata = nd[model_ix, , drop = FALSE], re_formula = NULL, allow_new_levels = TRUE, 
            sample_new_levels = "gaussian", ndraws = min(ndraws, posterior::ndraws(fit)))
        mu <- gd(prep, "mu", inv_link = TRUE)[, keep, drop = FALSE]
        sigma <- gd(prep, "sigma", inv_link = TRUE)[, keep, drop = FALSE]
        hu <- gd(prep, "hu", inv_link = TRUE)[, keep, drop = FALSE]
        p <- family_probabilities(fit, nd$NM_response[ix], mu, sigma, hu)
        sc <- score_family_probabilities(p)
        pm <- pit_from_score(sc)
        sc$row_id_global <- nd$row_id_global[ix]
        sc$pred_zero <- colMeans(hu)
        sc$log_pred_density <- sc$log_mass
        sc$pit_mean <- pm$mean
        sc$pit_second <- pm$second
        sc$pit_low <- pm$low
        sc$pit_high <- pm$high
        if (coverage) {
            pp <- posterior_predict(prep)[, keep, drop = FALSE]
            for (level in c(0.5, 0.8, 0.95)) {
                b <- apply(pp, 2, quantile, probs = c((1 - level)/2, (1 + level)/2), type = 1)
                lp <- family_probabilities(fit, b[1, ], mu, sigma, hu)
                hp <- family_probabilities(fit, b[2, ], mu, sigma, hu)
                sc[[paste0("coverage_", 100 * level)]] <- nd$NM_response[ix] >= b[1, ] & nd$NM_response[ix] <= 
                  b[2, ]
                sc[[paste0("expected_coverage_", 100 * level)]] <- 1 - colMeans(lp$lower) - colMeans(hp$upper)
            }
        }
        out[[length(out) + 1]] <- sc
    }
    bind_rows(out)
}


# Coverage uses predictive simulation;
# state probability/PIT moments use the deterministic integrated likelihood.
predict_normative <- function(fit,data,seed,ndraws=8000,coverage=FALSE,
 diagnostics_file=tempfile('swm_state_precision_',fileext='.csv')) {
 sc<-swm_state_predict_release(fit,data,'hdg',diagnostics_file)
 if(coverage) {
  sampled<-predict_normative_coverage(fit,data,seed,ndraws,coverage=TRUE)
  keys<-grep('^(coverage_|expected_coverage_)',names(sampled),value=TRUE)
  sc[keys]<-sampled[match(sc$row_id_global,sampled$row_id_global),keys,drop=FALSE]
  sc$coverage_method<-'sampled_new_person_predictive_distribution'
 }
 pm<-pit_from_score(sc)
 sc$log_pred_density<-sc$log_mass
 sc$pit_mean<-pm$mean;sc$pit_second<-pm$second;sc$pit_low<-pm$low;sc$pit_high<-pm$high
 sc
}
