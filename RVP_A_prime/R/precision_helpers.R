# Numerical scoring helpers. Fitted distributions and their priors are unchanged.
normal_quadrature <- function(n) {
  stopifnot(length(n)==1L,is.finite(n),n>=2,n==as.integer(n))
  J <- matrix(0,n,n);J[cbind(seq_len(n-1),2:n)]<-sqrt(seq_len(n-1));J<-J+t(J)
  e <- eigen(J,symmetric=TRUE)
  weight <- e$vectors[1,]^2
  list(node=e$values,weight=weight/sum(weight))
}

aprime_logadd <- function(a,b) {
  m <- pmax(a,b)
  ans <- m+log1p(exp(-abs(a-b)))
  ans[is.infinite(m)&m<0] <- -Inf
  ans
}
aprime_logmean <- function(x) {
  m <- max(x)
  if(is.infinite(m)&&m<0)return(-Inf)
  m+log(mean(exp(x-m)))
}
aprime_probability <- function(log_lower,log_upper) {
  ifelse(log_lower<=log(.5),exp(log_lower),-expm1(log_upper))
}
aprime_normal_quantile <- function(log_lower,log_upper) {
  # Only evaluate the smaller tail; its complement can round just above one.
  z <- qnorm(pmin(log_lower,log_upper),log.p=TRUE)
  ifelse(log_lower<=log_upper,z,-z)
}

random_effect_variance <- function(sd_intercept,sd_time,rho,time) {
  v <- sd_intercept^2+time^2*sd_time^2+2*time*rho*sd_intercept*sd_time
  if(any(!is.finite(v)|v < -1e-12))stop('Invalid random-effect covariance')
  pmax(v,0)
}

xbeta_bounds_logit <- function(y,eta,phi,kappa) {
  lens<-lengths(list(y,eta,phi,kappa));n<-max(lens)
  stopifnot(all(lens %in% c(1L,n)))
  y<-rep_len(y,n);eta<-rep_len(eta,n);phi<-rep_len(phi,n);kappa<-rep_len(kappa,n)
  stopifnot(all(is.finite(c(y,eta,phi,kappa))),all(y>=0&y<=1),all(phi>0),all(kappa>=0))
  # plogis(-eta) avoids cancellation in 1-plogis(eta) for positive tails.
  a<-plogis(eta)*phi;b<-plogis(-eta)*phi
  stopifnot(all(a>0),all(b>0))
  # Evaluate both tails using the smaller beta argument, which cannot round
  # to one when a tiny boundary width still carries substantial probability.
  lower_half<-y<=.5
  x<-(ifelse(lower_half,y,1-y)+kappa)/(1+2*kappa)
  shape1<-ifelse(lower_half,a,b);shape2<-ifelse(lower_half,b,a)
  log_cdf<-pbeta(x,shape1,shape2,log.p=TRUE)
  log_survival<-pbeta(x,shape1,shape2,lower.tail=FALSE,log.p=TRUE)
  ll<-ifelse(lower_half,log_cdf,log_survival)
  lu<-ifelse(lower_half,log_survival,log_cdf)
  lm<-rep(-Inf,n)
  lm[y==0]<-ll[y==0];ll[y==0]<--Inf
  lm[y==1]<-lu[y==1];lu[y==1]<--Inf
  list(lower=aprime_probability(ll,aprime_logadd(lm,lu)),
    upper=aprime_probability(aprime_logadd(ll,lm),lu),
    log_lower=ll,log_mass=lm,log_upper=lu)
}

aprime_cdf_map <- function(chunks,work,workers,os_type=.Platform$OS.type) {
  if(!length(chunks))return(list())
  # Windows R has no fork workers. Serial integration retains the same arithmetic.
  if(os_type=='windows'||workers<=1L)lapply(chunks,work) else
    parallel::mclapply(chunks,work,mc.cores=min(workers,length(chunks)),mc.preschedule=TRUE)
}

score_cdf_quadrature <- function(fit,data,role,nodes=c(32L,64L,128L,256L),
                                 tolerance=.002,mcse_limit=.03,workers=1L,chunk=8L) {
  stopifnot(nrow(data)>0,!anyDuplicated(data$row_id_global),all(is.finite(data$NM_response)))
  draws<-posterior::as_draws_matrix(fit);N<-nrow(draws);chains<-posterior::nchains(fit)
  re_names<-grep('^(sd|cor)_Subject_ID__',colnames(draws),value=TRUE)
  stopifnot(setequal(re_names,c('sd_Subject_ID__Intercept','sd_Subject_ID__Time_since_baseline_1y',
    'cor_Subject_ID__Intercept__Time_since_baseline_1y')),N%%chains==0)
  si<-draws[,'sd_Subject_ID__Intercept'];st<-draws[,'sd_Subject_ID__Time_since_baseline_1y']
  rho<-draws[,'cor_Subject_ID__Intercept__Time_since_baseline_1y']
  quadratures<-lapply(nodes,normal_quadrature)
  work<-function(ix) {
    nd<-data[ix,,drop=FALSE]
    pr<-brms::prepare_predictions(fit,newdata=nd,re_formula=NA,draw_ids=seq_len(N))
    rows<-lapply(seq_len(nrow(nd)),function(j) {
      lp<-brms::get_dpar(pr,'mu',i=j,inv_link=FALSE)
      ph<-brms:::get_dpar(pr,'phi',i=j);ka<-brms:::get_dpar(pr,'kappa',i=j)
      ss<-sqrt(random_effect_variance(si,st,rho,nd$Time_since_baseline_1y[j]))
      previous<-NA_real_;error<-Inf
      for(q in seq_along(nodes)) {
        g<-quadratures[[q]];ll<-lm<-lu<-rep(-Inf,N)
        keep<-which(g$weight>1e-20)
        lw<-log(g$weight[keep]/sum(g$weight[keep]))
        for(h in seq_along(keep)) {
          k<-keep[h]
          bounds<-xbeta_bounds_logit(nd$NM_response[j],lp+ss*g$node[k],ph,ka)
          ll<-aprime_logadd(ll,lw[h]+bounds$log_lower)
          lm<-aprime_logadd(lm,lw[h]+bounds$log_mass)
          lu<-aprime_logadd(lu,lw[h]+bounds$log_upper)
        }
        # Normalize the three disjoint events without reconstructing a tiny
        # tail from its complement. This removes quadrature roundoff at one.
        total<-aprime_logadd(aprime_logadd(ll,lm),lu)
        ll<-ll-total;lm<-lm-total;lu<-lu-total
        log_p<-aprime_logadd(ll,lm-log(2))
        log_q<-aprime_logadd(lu,lm-log(2))
        mean_p<-aprime_logmean(log_p);mean_q<-aprime_logmean(log_q)
        z<-aprime_normal_quantile(mean_p,mean_q)
        if(!is.finite(z))stop('Unresolved marginal probability: row ',nd$row_id_global[j])
        if(q>1L){error<-abs(z-previous);if(error<=tolerance)break}
        previous<-z
      }
      # This includes posterior-chain Monte Carlo variation. It is not a subject-level CI.
      tail<-if(mean_p<=mean_q)log_p else log_q
      scale<-max(tail);centered<-expm1(tail-scale);spread<-max(abs(centered))
      # Center/rescale before estimating MCSE so nearly constant probabilities
      # do not disappear in a subtraction or underflow in variance calculations.
      se<-if(spread==0)0 else
        as.numeric(posterior::mcse_mean(matrix(centered/spread,ncol=chains)))*spread
      mcse<-if(isTRUE(se==0))0 else exp(log(se)+scale-dnorm(z,log=TRUE))
      lower<-aprime_logmean(ll);mass<-aprime_logmean(lm);upper<-aprime_logmean(lu)
      data.frame(row_id_global=nd$row_id_global[j],role=role,
        p_lower=aprime_probability(lower,aprime_logadd(mass,upper)),
        p_upper=aprime_probability(aprime_logadd(lower,mass),upper),
        log_lower=lower,log_mass=mass,log_upper=upper,
        z=z,n_posterior_draws=N,quadrature_nodes=nodes[q],quadrature_z_difference=error,
        posterior_mcse_z=mcse,precision_status=if(error>tolerance)'quadrature_warning'else
          if(!is.finite(mcse)||mcse>mcse_limit)'posterior_MC_warning'else'ok')
    })
    do.call(rbind,rows)
  }
  chunks<-split(seq_len(nrow(data)),ceiling(seq_len(nrow(data))/chunk))
  ans<-aprime_cdf_map(chunks,work,workers)
  if(any(vapply(ans,inherits,logical(1),'try-error')))stop('Quadrature worker failed')
  out<-do.call(rbind,ans);out<-out[match(data$row_id_global,out$row_id_global),]
  stopifnot(identical(out$row_id_global,data$row_id_global))
  out
}
