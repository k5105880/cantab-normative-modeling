# Convert predictive probabilities to normative scores and calibration moments.
log_add <- function(a, b) {
    m <- pmax(a, b)
    z <- m + log1p(exp(pmin(a, b) - m))
    z[is.infinite(m) & m < 0] <- -Inf
    z
}

log_mean_columns <- function(x) apply(x, 2, function(v) {
    m <- max(v)
    if (is.infinite(m) && m < 0) 
        -Inf
    else m + log(mean(exp(v - m)))
})

score_family_probabilities <- function(p) {
    if (is.null(p$log_mass)) {
        p$log_lower <- log(p$lower)
        p$log_mass <- log(p$mass)
        p$log_upper <- log(p$upper)
    }
    ll <- log_mean_columns(p$log_lower)
    lm <- log_mean_columns(p$log_mass)
    lu <- log_mean_columns(p$log_upper)
    if (any(!is.finite(lm))) 
        stop("Predictive log mass is nonfinite; inspect fitted parameters and input support.")
    pl <- pmin(0, log_add(ll, lm - log(2)))
    pu <- pmin(0, log_add(lu, lm - log(2)))
    z <- ifelse(pl <= log(0.5), qnorm(pl, log.p = TRUE), qnorm(pu, lower.tail = FALSE, log.p = TRUE))
    data.frame(lower = exp(ll), mass = exp(lm), upper = exp(lu), mid_p = exp(pl), z = z, log_lower = ll, log_mass = lm, 
        log_upper = lu)
}

pit_from_score <- function(sc) {
    good <- sc$mass > 1e-250
    out <- data.frame(mean = rep(NA_real_, nrow(sc)), second = NA_real_, low = NA_real_, high = NA_real_)
    if (any(good)) 
        out[good, ] <- pit_moments(sc$lower[good], sc$mass[good], sc$upper[good])
    for (i in which(!good)) {
        ll <- sc$log_lower[i]
        lm <- sc$log_mass[i]
        lu <- sc$log_upper[i]
        if (lm < min(ll, lu) - log(1e+09)) {
            m <- sc$z[i]
            e2 <- m * m
        }
        else {
            a <- if (ll <= log(0.5)) 
                qnorm(ll, log.p = TRUE)
            else qnorm(pmin(0, log_add(lu, lm)), lower.tail = FALSE, log.p = TRUE)
            b <- if (lu <= log(0.5)) 
                qnorm(lu, lower.tail = FALSE, log.p = TRUE)
            else qnorm(pmin(0, log_add(ll, lm)), log.p = TRUE)
            ra <- if (is.finite(a)) 
                exp(dnorm(a, log = TRUE) - lm)
            else 0
            rb <- if (is.finite(b)) 
                exp(dnorm(b, log = TRUE) - lm)
            else 0
            m <- ra - rb
            e2 <- 1 + if (is.finite(a)) 
                a * ra
            else 0
            e2 <- e2 - if (is.finite(b)) 
                b * rb
            else 0
        }
        out[i, ] <- c(m, max(e2, m * m), as.numeric(sc$z[i] < -1.96), as.numeric(sc$z[i] > 1.96))
    }
    stopifnot(all(is.finite(as.matrix(out))))
    out
}

