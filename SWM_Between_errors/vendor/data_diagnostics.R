# Sampling diagnostics and acceptance thresholds.
between_parameter_diagnostics <- function(draws) {
    vars <- posterior::variables(draws)
    structural <- invalid_structure <- rep(FALSE, length(vars))
    for (k in seq_along(vars)) {
        hit <- regmatches(vars[k], regexec('^L_[0-9]+\\[([0-9]+),([0-9]+)\\]$', vars[k]))[[1]]
        if (length(hit) != 3L) next
        i <- as.integer(hit[2]); j <- as.integer(hit[3])
        expected <- if (i == 1L && j == 1L) 1 else if (i < j) 0 else next
        x <- as.numeric(draws[,,vars[k]])
        structural[k] <- all(is.finite(x)) && all(x == expected)
        invalid_structure[k] <- !structural[k]
    }
    checked <- draws[,,!structural,drop=FALSE]
    if (!length(posterior::variables(checked)))
        return(data.frame(max_rhat=Inf,min_ess_bulk=0,min_ess_tail=0,
            n_parameters_without_finite_diagnostics=1L,
            parameters_without_finite_diagnostics='no sampled parameters',
            n_structural_constants=sum(structural)))
    s <- posterior::summarise_draws(checked, 'rhat', 'ess_bulk', 'ess_tail')
    finite <- vapply(s$variable, function(v) all(is.finite(checked[,,v])), logical(1)) &
        is.finite(s$rhat) & is.finite(s$ess_bulk) & is.finite(s$ess_tail) &
        !s$variable %in% vars[invalid_structure]
    data.frame(max_rhat=if(any(finite))max(s$rhat[finite])else Inf,
        min_ess_bulk=if(any(finite))min(s$ess_bulk[finite])else 0,
        min_ess_tail=if(any(finite))min(s$ess_tail[finite])else 0,
        n_parameters_without_finite_diagnostics=sum(!finite),
        parameters_without_finite_diagnostics=paste(s$variable[!finite],collapse=';'),
        n_structural_constants=sum(structural))
}

diagnostics <- function(fit, max_treedepth=fit$stan_args$control$max_treedepth) {
    draws <- as_draws_array(fit)
    s <- between_parameter_diagnostics(draws)
    np <- nuts_params(fit)
    pars <- c('energy__','divergent__','treedepth__')
    valid <- length(max_treedepth)==1L && is.finite(max_treedepth) &&
        max_treedepth>=1 && all(c('Parameter','Value','Chain','Iteration') %in% names(np))
    if (valid) {
        z <- np[np$Parameter %in% pars,,drop=FALSE]
        valid <- nrow(z)==prod(dim(draws)[1:2])*length(pars) &&
            all(is.finite(z$Value)) && !anyNA(z[c('Chain','Iteration')]) &&
            !anyDuplicated(z[c('Parameter','Chain','Iteration')]) &&
            length(unique(z$Chain))==dim(draws)[2]
        if (valid) valid <- all(table(z$Chain,factor(z$Parameter,levels=pars))==dim(draws)[1])
    }
    if (!valid) return(cbind(s,data.frame(divergences=NA_real_,treedepth_hits=NA_real_,min_EBFMI=NA_real_)))
    z <- z[order(z$Chain,z$Iteration),,drop=FALSE]
    energies <- split(z$Value[z$Parameter=='energy__'], z$Chain[z$Parameter=='energy__'])
    bfmi <- vapply(energies,function(e)if(length(e)<2L || !is.finite(var(e)) || var(e)<=0)
        NA_real_ else mean(diff(e)^2)/var(e),numeric(1))
    cbind(s,data.frame(divergences=sum(z$Value[z$Parameter=='divergent__']),
        treedepth_hits=sum(z$Value[z$Parameter=='treedepth__']>=max_treedepth),
        min_EBFMI=min(bfmi)))
}

diagnostic_pass <- function(d) {
    required <- c('max_rhat','min_ess_bulk','min_ess_tail','divergences',
        'treedepth_hits','min_EBFMI','n_parameters_without_finite_diagnostics')
    if (!is.data.frame(d) || nrow(d)!=1L || !all(required %in% names(d)) ||
        !all(vapply(d[required],is.numeric,logical(1))) ||
        !all(is.finite(unlist(d[required])))) return(FALSE)
    d$max_rhat<1.01 && d$min_ess_bulk>=400 && d$min_ess_tail>=400 &&
        d$divergences==0 && d$treedepth_hits==0 && d$min_EBFMI>0.3 &&
        d$n_parameters_without_finite_diagnostics==0
}
