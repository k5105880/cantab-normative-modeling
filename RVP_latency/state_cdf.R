# Marginal normative scoring. Shared numerical method for latency and movement.
# Conditional likelihood CDF -> Gaussian new-person integration -> posterior
# probability average -> inverse normal. No simulated residual/new-person draws.
state_cdf_version <- 'analytic_state_cdf_v1'
state_cdf_settings <- list(quad_nodes=32L, check_nodes=64L, max_nodes=1024L,
  z_tolerance=1e-4, tail_relative_tolerance=1e-4, draw_probability_tolerance=1e-6,
  posterior_mcse_z_tolerance=.02, chunk_size=25L)

.state_integer <- function(x, name, lower=1L) {
  if(length(x)!=1L || !is.numeric(x) || !is.finite(x) || x!=floor(x) || x<lower)
    stop(name,' must be an integer >= ',lower,call.=FALSE)
  as.integer(x)
}
state_gauss_normal <- local({
  cache<-new.env(parent=emptyenv())
  function(n) {
    n<-.state_integer(n,'quadrature nodes',2L);key<-as.character(n)
    if(exists(key,cache,inherits=FALSE))return(cache[[key]])
    J<-matrix(0,n,n);J[cbind(1:(n-1L),2:n)]<-sqrt(1:(n-1L))
    e<-eigen(J+t(J),symmetric=TRUE);ix<-order(e$values)
    w<-e$vectors[1,ix]^2;keep<-w>0
    cache[[key]]<-list(nodes=e$values[ix][keep],weights=w[keep]/sum(w[keep]))
    cache[[key]]
  }
})
.state_logadd <- function(a,b) {
  m<-pmax(a,b);out<-m+log1p(exp(-abs(a-b)))
  out[is.infinite(m)&m<0]<--Inf;out
}
.state_logmean <- function(x) {
  m<-apply(x,2,max);ans<-m+log(colMeans(exp(sweep(x,2,m,'-'))))
  ans[is.infinite(m)&m<0]<--Inf;pmin(0,ans)
}
.z_from_log_probabilities <- function(ll,lu) {
  ifelse(ll<=lu,qnorm(ll,log.p=TRUE),qnorm(lu,lower.tail=FALSE,log.p=TRUE))
}

# aux is shape for Gamma/log mean, sigma for Student/identity mean.
state_cdf_components <- function(response, eta, aux, nu, re_variance, family_id,nodes) {
  if(!is.matrix(eta)||!is.numeric(eta)||any(dim(eta)<1L)||
     !is.matrix(aux)||!identical(dim(aux),dim(eta))||
     !is.matrix(re_variance)||!identical(dim(re_variance),dim(eta))||
     length(response)!=ncol(eta))stop('Incompatible draws-by-rows dimensions',call.=FALSE)
  if(any(!is.finite(response))||any(!is.finite(eta))||any(!is.finite(aux))||
     any(!is.finite(re_variance)))stop('All state CDF inputs must be finite',call.=FALSE)
  if(any(aux<=0))stop('Likelihood scale/shape must be positive',call.=FALSE)
  if(any(re_variance<0))stop('Random-effect variance must be nonnegative',call.=FALSE)
  if(!family_id %in% c('gamma_raw','student_log1p'))stop('Unsupported analytic state family',call.=FALSE)
  if(family_id=='gamma_raw'&&any(response<=0))stop('Gamma response must be positive',call.=FALSE)
  if(family_id=='student_log1p'&&(length(nu)!=nrow(eta)||any(!is.finite(nu))||any(nu<=0)))
    stop('Student nu must be finite and positive per draw',call.=FALSE)
  g<-state_gauss_normal(nodes);sd_re<-sqrt(re_variance)
  y<-matrix(response,nrow(eta),ncol(eta),byrow=TRUE)
  ll<-lu<-matrix(-Inf,nrow(eta),ncol(eta))
  for(k in seq_along(g$nodes)) {
    loc<-eta+sd_re*g$nodes[k]
    if(family_id=='gamma_raw') {
      # Unit-rate argument avoids overflow from separately computing exp(eta).
      q<-exp(log(y)+log(aux)-loc)
      low<-pgamma(q,shape=aux,log.p=TRUE)
      up<-pgamma(q,shape=aux,lower.tail=FALSE,log.p=TRUE)
    } else {
      q<-(y-loc)/aux
      low<-pt(q,df=nu,log.p=TRUE);up<-pt(q,df=nu,lower.tail=FALSE,log.p=TRUE)
    }
    ll<-.state_logadd(ll,matrix(low,nrow(eta),ncol(eta))+log(g$weights[k]))
    lu<-.state_logadd(lu,matrix(up,nrow(eta),ncol(eta))+log(g$weights[k]))
  }
  list(lower=ll,upper=lu,log_lower=.state_logmean(ll),log_upper=.state_logmean(lu))
}

# Means are diagnosed with original chains/iteration order; no IID assumption or
# pseudochains from arbitrarily batched/thinned predictions. Scaling by each
# row's maximum log tail avoids underflow in extreme tails. MCSE is a first-order
# delta-method estimate for z; it is NOT quadrature error or posterior SD of z.
state_posterior_precision <- function(log_tail,metadata,log_mean_tail,z) {
  if(!requireNamespace('posterior',quietly=TRUE))stop('posterior package required')
  if(!is.matrix(log_tail)||!is.data.frame(metadata)||nrow(metadata)!=nrow(log_tail)||
     !all(c('.chain','.iteration','.draw') %in% names(metadata))||
     anyNA(metadata)||anyDuplicated(metadata$.draw)||
     anyDuplicated(metadata[c('.chain','.iteration')]))stop('Invalid posterior chain metadata',call.=FALSE)
  chains<-sort(unique(metadata$.chain));ids<-lapply(chains,function(k) {
    i<-which(metadata$.chain==k);i[order(metadata$.iteration[i])]
  })
  balanced<-length(unique(lengths(ids)))==1L
  contiguous<-all(vapply(ids,function(i) length(i)>=4L && all(diff(metadata$.iteration[i])==1L),logical(1)))
  available<-balanced&&contiguous
  rows<-lapply(seq_len(ncol(log_tail)),function(j) {
    x<-exp(log_tail[,j]-max(log_tail[,j]));scale<-exp(max(log_tail[,j]))
    if(any(!is.finite(x))||!is.finite(log_mean_tail[j])||!is.finite(z[j]))stop('Nonfinite tail precision inputs')
    se<-ess<-rh<-NA_real_;method<-'unavailable_short_or_incomplete_chains'
    if(available) {
      a<-vapply(ids,function(i)x[i],numeric(length(ids[[1]])))
      if(diff(range(x))==0) {se<-0;ess<-length(x);method<-'constant_draw_contribution'} else {
        se<-as.numeric(posterior::mcse_mean(a));ess<-as.numeric(posterior::ess_mean(a));rh<-as.numeric(posterior::rhat(a))
        method<-'chain_aware_spectral_ess_delta_method'
      }
    }
    rel<-se/mean(x)
    zse<-rel*exp(log_mean_tail[j]-dnorm(z[j],log=TRUE))
    good<-is.finite(zse)&&zse<=state_cdf_settings$posterior_mcse_z_tolerance
    data.frame(posterior_draws=nrow(log_tail),posterior_chains=length(chains),
      posterior_mcse_p=se*scale,posterior_relative_mcse_tail=rel,
      posterior_mcse_z=zse,posterior_ess_mean=ess,posterior_rhat=rh,
      posterior_precision_ok=good,posterior_mcse_method=method)
  })
  do.call(rbind,rows)
}

state_cdf_checked <- function(response,eta,aux,nu,re_variance,family_id,metadata,
  quad_nodes=state_cdf_settings$quad_nodes,check_nodes=state_cdf_settings$check_nodes,
  max_nodes=state_cdf_settings$max_nodes,z_tolerance=state_cdf_settings$z_tolerance,
  tail_relative_tolerance=state_cdf_settings$tail_relative_tolerance,
  draw_probability_tolerance=state_cdf_settings$draw_probability_tolerance) {
  quad_nodes<-.state_integer(quad_nodes,'quad_nodes',2);check_nodes<-.state_integer(check_nodes,'check_nodes',2)
  max_nodes<-.state_integer(max_nodes,'max_nodes',2)
  if(quad_nodes>=check_nodes||check_nodes>max_nodes)stop('Require quad_nodes < check_nodes <= max_nodes')
  for(tol in list(z_tolerance,tail_relative_tolerance,draw_probability_tolerance))
    if(length(tol)!=1L||!is.finite(tol)||tol<=0)stop('Invalid quadrature tolerance')
  calc<-function(ix,n) state_cdf_components(response[ix],eta[,ix,drop=FALSE],aux[,ix,drop=FALSE],nu,re_variance[,ix,drop=FALSE],family_id,n)
  compare<-function(a,b,ix,n) {
    ratio<-if(family_id=='gamma_raw')sqrt(re_variance[,ix,drop=FALSE]*aux[,ix,drop=FALSE]) else sqrt(re_variance[,ix,drop=FALSE])/aux[,ix,drop=FALSE]
    resolution<-apply(ratio,2,max)*pi/sqrt(n)
    zd<-abs(.z_from_log_probabilities(a$log_lower,a$log_upper)-.z_from_log_probabilities(b$log_lower,b$log_upper))
    taildiff<-pmax(abs(expm1(a$log_lower-b$log_lower)),abs(expm1(a$log_upper-b$log_upper)))
    drawdiff<-pmax(apply(abs(exp(a$lower)-exp(b$lower)),2,max),apply(abs(exp(a$upper)-exp(b$upper)),2,max))
    good<-is.finite(resolution)&resolution<=1&is.finite(zd)&zd<=z_tolerance&is.finite(taildiff)&taildiff<=tail_relative_tolerance&
      is.finite(drawdiff)&drawdiff<=draw_probability_tolerance&
      abs(exp(b$log_lower)+exp(b$log_upper)-1)<1e-10
    data.frame(quad_z_diff=zd,quad_relative_tail_diff=taildiff,
      quad_max_draw_probability_diff=drawdiff,quad_step_to_scale=resolution,numerical_ok=good)
  }
  ix<-seq_len(ncol(eta));coarse<-calc(ix,quad_nodes);fine<-calc(ix,check_nodes)
  diag<-compare(coarse,fine,ix,check_nodes);used<-rep(check_nodes,length(ix));order<-check_nodes
  while(any(!diag$numerical_ok)&&order<max_nodes) {
    bad<-which(!diag$numerical_ok);order<-min(order*2L,max_nodes)
    nextfine<-calc(bad,order)
    previous<-lapply(fine,function(x)if(is.matrix(x))x[,bad,drop=FALSE] else x[bad])
    diag[bad,]<-compare(previous,nextfine,bad,order)
    for(nm in names(fine)) {
      if(is.matrix(fine[[nm]]))fine[[nm]][,bad]<-nextfine[[nm]] else fine[[nm]][bad]<-nextfine[[nm]]
    }
    used[bad]<-order
  }
  if(any(!diag$numerical_ok))stop('State CDF quadrature did not converge for rows ',paste(which(!diag$numerical_ok),collapse=','),
    ' at max_nodes=',max_nodes,'; no empirical-CDF fallback is permitted',call.=FALSE)
  z<-.z_from_log_probabilities(fine$log_lower,fine$log_upper)
  if(any(!is.finite(z)))stop('Nonfinite analytic state z; scores withheld',call.=FALSE)
  lower_tail<-fine$log_lower<=fine$log_upper
  tails<-fine$lower;tails[,!lower_tail]<-fine$upper[,!lower_tail,drop=FALSE]
  precision<-state_posterior_precision(tails,metadata,pmin(fine$log_lower,fine$log_upper),z)
  cbind(data.frame(p_lower=exp(fine$log_lower),p_upper=exp(fine$log_upper),
    log_p_lower=fine$log_lower,log_p_upper=fine$log_upper,z=z,
    clipped=FALSE,probability_rounded_boundary=exp(fine$log_lower)==0|exp(fine$log_lower)==1,
    quad_nodes_used=used),diag,precision)
}

state_re_variance <- function(cov_draws,newdata) {
  d<-dim(cov_draws);cn<-dimnames(cov_draws)[[2]]
  if(length(d)!=3L||d[1]<1L||d[2]!=d[3]||any(!is.finite(cov_draws))||
     is.null(cn)||!identical(cn,dimnames(cov_draws)[[3]])||anyDuplicated(cn)||
     !(setequal(cn,'Intercept')||setequal(cn,c('Intercept','Time_since_baseline_1y'))))
    stop('Unsupported Subject_ID covariance dimensions or coefficients',call.=FALSE)
  tol<-1e-10*max(1,abs(cov_draws))
  if(max(abs(cov_draws-aperm(cov_draws,c(1,3,2))))>tol)stop('Covariance must be symmetric')
  for(i in seq_len(d[1])) {
    if(min(eigen(matrix(cov_draws[i,,],d[2],d[3]),symmetric=TRUE,only.values=TRUE)$values) < -tol)
      stop('Covariance must be positive semidefinite',call.=FALSE)
  }
  z<-matrix(1,nrow(newdata),length(cn),dimnames=list(NULL,cn))
  if('Time_since_baseline_1y'%in%cn) {
    t<-newdata$Time_since_baseline_1y
    if(length(t)!=nrow(newdata)||any(!is.finite(t)))stop('Finite time required for random slope')
    z[,'Time_since_baseline_1y']<-t
  }
  v<-matrix(0,d[1],nrow(newdata))
  for(j in seq_along(cn))for(k in seq_along(cn))v<-v+outer(cov_draws[,j,k],z[,j]*z[,k])
  if(any(!is.finite(v))||any(v< -tol))stop('Invalid projected new-person variance')
  pmax(v,0)
}

state_validate_fit <- function(fit,family_id) {
  if(!inherits(fit,'brmsfit'))stop('Expected brmsfit')
  expected<-switch(family_id,gamma_raw=c('gamma','log'),student_log1p=c('student','identity'),stop('Unsupported analytic state family'))
  if(!identical(tolower(fit$family$family),expected[1])||!identical(fit$family$link,expected[2]))stop('Family/link does not match analytic state scorer')
  tt<-brms::brmsterms(fit$formula)
  aux<-if(family_id=='gamma_raw')'shape' else 'sigma'
  if(!identical(fit$formula$formula[[2]],as.name('NM_response'))||length(tt$nlpars)||length(fit$formula$pfix)||
     !all(names(tt$dpars)%in%c('mu',aux)))stop('Unsupported response/distributional structure')
  re<-tt$dpars$mu$re
  if(is.null(re)||nrow(re)!=1L||!identical(as.character(re$group),'Subject_ID'))stop('Exactly one Subject_ID term required')
  gc<-re$gcall[[1]]
  if(!identical(gc$dist,'gaussian')||any(vapply(c('by','pw','cov','type'),function(x)!identical(gc[[x]],''),logical(1))))stop('Only ordinary Gaussian subject effects supported')
  allowed_mu<-c('formula','fe','re','sm','allvars','family','dpar','resp','respform','transform')
  allowed_aux<-c('formula','fe','allvars','family','dpar','resp','respform','transform')
  if(!all(names(tt$dpars$mu)%in%allowed_mu)||!all(names(tt$dpars[[aux]])%in%allowed_aux))stop('Unsupported special/random auxiliary effects')
  invisible(TRUE)
}

state_predict_marginal <- function(fit,newdata,family_id,chunk_size=state_cdf_settings$chunk_size) {
  state_validate_fit(fit,family_id);chunk_size<-.state_integer(chunk_size,'chunk_size')
  if(!is.data.frame(newdata)||!nrow(newdata)||!all(c('NM_response','row_id_global','Subject_ID')%in%names(newdata))||
     anyNA(newdata$row_id_global)||anyDuplicated(newdata$row_id_global))stop('Scoring data require unique row IDs, response and Subject_ID')
  # ALL retained posterior draws, in canonical brms ordering, with chain identity.
  parameter<-if(family_id=='student_log1p')'nu' else 'lp__'
  draws<-posterior::as_draws_df(fit,variable=parameter)
  metadata<-as.data.frame(draws)[c('.chain','.iteration','.draw')]
  if(!identical(as.integer(metadata$.draw),seq_len(nrow(metadata))))stop('Unexpected posterior draw ordering')
  nu<-if(family_id=='student_log1p')as.numeric(draws$nu) else NULL
  vc<-brms::VarCorr(fit,summary=FALSE)
  if(!identical(names(vc),'Subject_ID'))stop('Exactly one Subject_ID covariance block required')
  covariance<-vc$Subject_ID$cov
  # brms returns only sd (no cov member) for a one-coefficient random intercept.
  if(is.null(covariance)) {
    sd<-vc$Subject_ID$sd
    if(!is.matrix(sd)||ncol(sd)!=1L||!identical(colnames(sd),'Intercept'))
      stop('Unsupported scalar Subject_ID random effect')
    covariance<-array(sd^2,c(nrow(sd),1L,1L),dimnames=list(NULL,'Intercept','Intercept'))
  }
  if(dim(covariance)[1]!=nrow(metadata))stop('Covariance and posterior draws are not aligned')
  variance<-state_re_variance(covariance,newdata)
  result<-list()
  for(start in seq.int(1L,nrow(newdata),by=chunk_size)) {
    idx<-seq.int(start,min(start+chunk_size-1L,nrow(newdata)));nd<-newdata[idx,,drop=FALSE]
    # mgcv 1.9-4 drops matrix dimensions for a singleton tp/m=c(2,0)
    # prediction. Duplicate only for basis evaluation, then discard padding.
    prediction_data<-if(nrow(nd)==1L)nd[c(1L,1L),,drop=FALSE] else nd
    eta<-brms::posterior_linpred(fit,newdata=prediction_data,re_formula=NA,allow_new_levels=TRUE)
    aux<-brms::posterior_linpred(fit,newdata=prediction_data,re_formula=NA,allow_new_levels=TRUE,
      dpar=if(family_id=='gamma_raw')'shape' else 'sigma',transform=TRUE)
    eta<-eta[,seq_len(nrow(nd)),drop=FALSE];aux<-aux[,seq_len(nrow(nd)),drop=FALSE]
    if(nrow(eta)!=nrow(metadata))stop('Linear predictions and chain metadata are not aligned')
    out<-state_cdf_checked(nd$NM_response,eta,aux,nu,variance[,idx,drop=FALSE],family_id,metadata)
    result[[length(result)+1L]]<-cbind(data.frame(row_id_global=nd$row_id_global),out)
  }
  out<-do.call(rbind,result);rownames(out)<-NULL;out
}

state_predict_oof <- function(kf_obj,model_dat,family_id) {
  if(is.null(kf_obj$fits)||!all(c('fit','omitted')%in%colnames(kf_obj$fits)))stop('OOF saved fits and omitted indices required')
  result<-list();all_idx<-integer()
  for(k in seq_len(nrow(kf_obj$fits))) {
    fit<-kf_obj$fits[[k,'fit']];idx<-kf_obj$fits[[k,'omitted']]
    if(!is.numeric(idx)||anyNA(idx)||any(idx!=floor(idx))||any(idx<1|idx>nrow(model_dat))||anyDuplicated(idx))stop('Invalid OOF omitted rows')
    idx<-as.integer(idx);all_idx<-c(all_idx,idx)
    nd<-model_dat[idx,,drop=FALSE]
    if(any(as.character(nd$Subject_ID)%in%as.character(fit$data$Subject_ID)))stop('Subject leakage in analytic OOF scoring')
    out<-state_predict_marginal(fit,nd,family_id);out$fold<-k;result[[k]]<-out
  }
  if(!identical(sort(all_idx),seq_len(nrow(model_dat))))stop('OOF omitted rows must form a unique complete partition')
  out<-do.call(rbind,result);rownames(out)<-NULL;out
}

state_write_diagnostics <- function(result,prefix,family_id,label) {
  result$scoring_method<-state_cdf_version
  write.csv(result,paste0(prefix,'_',family_id,'_',label,'_state_cdf_diagnostics.csv'),row.names=FALSE,na='')
  if(any(!result$posterior_precision_ok))stop(sum(!result$posterior_precision_ok),
    ' state rows have posterior MCSE(z) > ',state_cdf_settings$posterior_mcse_z_tolerance,
    ' or unavailable MCSE; inspect state_cdf_diagnostics.csv. Scores withheld; increase posterior sampling. Quadrature convergence is separate.',call.=FALSE)
  invisible(result)
}

# Preserve typed numerical scoring failures through orchestration.
state_scoring_guard <- function(expr) {
  tryCatch(force(expr),error=function(e) {
    class(e)<-unique(c('state_cdf_error',class(e)));stop(e)
  })
}
