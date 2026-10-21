# Numerical scoring helpers. Fitted distributions and their priors are unchanged.
normal_quadrature <- function(n) {
  stopifnot(length(n)==1L,is.finite(n),n>=2,n==as.integer(n))
  J <- matrix(0,n,n);J[cbind(seq_len(n-1),2:n)]<-sqrt(seq_len(n-1));J<-J+t(J)
  e <- eigen(J,symmetric=TRUE)
  list(node=e$values,weight=e$vectors[1,]^2)
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
  f<-pbeta((y+kappa)/(1+2*kappa),a,b)
  low<-high<-f;low[y==0]<-0;high[y==1]<-1
  list(lower=low,upper=high)
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
      mu<-brms:::get_dpar(pr,'mu',i=j);lp<-qlogis(mu)
      ph<-brms:::get_dpar(pr,'phi',i=j);ka<-brms:::get_dpar(pr,'kappa',i=j)
      ss<-sqrt(random_effect_variance(si,st,rho,nd$Time_since_baseline_1y[j]))
      previous<-NA_real_;error<-Inf
      for(q in seq_along(nodes)) {
        g<-quadratures[[q]];low<-high<-numeric(N)
        for(k in which(g$weight>1e-20)) {
          bounds<-xbeta_bounds_logit(nd$NM_response[j],lp+ss*g$node[k],ph,ka)
          low<-low+g$weight[k]*bounds$lower;high<-high+g$weight[k]*bounds$upper
        }
        pd<-(low+high)/2;p<-mean(pd);z<-qnorm(p)
        if(!is.finite(z))stop('Unresolved marginal probability: row ',nd$row_id_global[j])
        if(q>1L){error<-abs(z-previous);if(error<=tolerance)break}
        previous<-z
      }
      # This includes posterior-chain Monte Carlo variation. It is not a subject-level CI.
      mcse<-as.numeric(posterior::mcse_mean(matrix(pd,ncol=chains)))/dnorm(z)
      data.frame(row_id_global=nd$row_id_global[j],role=role,p_lower=mean(low),p_upper=mean(high),
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
