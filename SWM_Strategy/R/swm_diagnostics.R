swm_rows <- function(pmf,d,model,uniforms,log_pmf=NULL) {
  stopifnot(length(pmf)==nrow(d),length(uniforms)==nrow(d))
  scores<-do.call(rbind,lapply(seq_len(nrow(d)),function(i)
    swm_score_pmf(pmf[[i]],d$index[i],d$grid[i],uniforms[i],
      if(is.null(log_pmf))NULL else log_pmf[[i]])))
  cbind(model=model,d[,c('row_id_global','Subject_ID','Platform','CV_fold','Age_band',
    'Age_baseline_years','Practice_bin','score','index','grid')],scores)
}

swm_parameter_diagnostics <- function(draws) {
  variables<-posterior::variables(draws)
  # brms saves the full Cholesky correlation matrix when save_pars(all=TRUE).
  # Its first diagonal entry is exactly 1 and its upper triangle exactly 0.
  # Verify both the structural identity and all draws; never exempt arbitrary
  # constant sampled parameters, which may instead indicate a stuck chain.
  structural<-vapply(variables,function(v) {
    hit<-regmatches(v,regexec('^L_[0-9]+\\[([0-9]+),([0-9]+)\\]$',v))[[1]]
    if(length(hit)!=3L)return(FALSE)
    i<-as.integer(hit[2]);j<-as.integer(hit[3])
    expected<-if(i==1L && j==1L)1 else if(i<j)0 else return(FALSE)
    x<-as.numeric(draws[,,v])
    all(is.finite(x)) && all(x==expected)
  },logical(1))
  if(all(structural))stop('No sampled quantities available for convergence diagnostics')
  checked<-draws[,,!structural,drop=FALSE]
  s<-posterior::summarise_draws(checked,'rhat','ess_bulk','ess_tail')
  finite_draws<-vapply(s$variable,function(v)all(is.finite(checked[,,v])),logical(1))
  finite<-finite_draws & is.finite(s$rhat)&is.finite(s$ess_bulk)&is.finite(s$ess_tail)
  if(!any(finite))stop('No usable convergence diagnostics')
  data.frame(max_rhat=max(s$rhat[finite]),min_ess_bulk=min(s$ess_bulk[finite]),
    min_ess_tail=min(s$ess_tail[finite]),n_parameters_without_finite_diagnostics=sum(!finite),
    n_structural_constants=sum(structural),
    structural_constants=paste(variables[structural],collapse=';'),
    parameters_without_finite_diagnostics=paste(s$variable[!finite],collapse=';'))
}

swm_diagnostics <- function(fit,label,max_treedepth=20) {
  n<-brms::nuts_params(fit)
  s<-swm_parameter_diagnostics(posterior::as_draws_array(fit))
  cbind(data.frame(label=label,n_divergent=sum(n$Parameter=='divergent__' & n$Value==1),
    treedepth_hits=sum(n$Parameter=='treedepth__' & n$Value>=max_treedepth)),s)
}
