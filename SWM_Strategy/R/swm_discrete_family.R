# Category indices are a likelihood encoding, NOT binomial counts.
# The scientific score remains RS6 raw/4 or Connect raw on 2--14.
swm_encode <- function(raw, platform) {
  if (length(raw) != length(platform) || any(!is.finite(raw)) ||
      any(!platform %in% c('RS6','Connect'))) stop('Invalid raw score/platform')
  divisor <- ifelse(platform == 'RS6',4,1)
  score <- raw/divisor
  if(any(abs(raw-round(raw))>1e-8) || any(score<2 | score>14))
    stop('Scores must be integer raw values within RS6 8--56 / Connect 2--14')
  data.frame(score=score,index=as.integer(round((score-2)*divisor)),
             grid=as.integer(12*divisor))
}

swm_logdiff <- function(a,b) {
  out <- rep(-Inf,length(a)); ok <- is.finite(a) & (b<a)
  out[ok] <- a[ok] + log(-expm1(b[ok]-a[ok])); out
}

swm_logpmf <- function(k, mu, phi, grid, zi=0) {
  values<-list(k,mu,phi,grid,zi);n<-max(lengths(values))
  if(any(!lengths(values)%in%c(1L,n))) stop('Incompatible probability vector lengths')
  values<-lapply(values,rep_len,length.out=n)
  k<-values[[1]];mu<-values[[2]];phi<-values[[3]];grid<-values[[4]];zi<-values[[5]]
  if(any(!is.finite(c(k,mu,phi,grid,zi))) || any(mu<=0|mu>=1) ||
     any(phi<=0) || any(grid<1|grid!=round(grid)) || any(zi<0|zi>1))
    stop('Invalid distribution parameters')
  valid<-k>=0 & k<=grid & k==round(k)
  lo<-pmax(0,(k-.5)/grid);hi<-pmin(1,(k+.5)/grid)
  a<-mu*phi;b<-(1-mu)*phi;lp<-rep(-Inf,n)
  left<-valid & k==0;right<-valid & k==grid;middle<-valid & !left & !right
  lp[left]<-pbeta(hi[left],a[left],b[left],log.p=TRUE)
  lp[right]<-pbeta(lo[right],a[right],b[right],lower.tail=FALSE,log.p=TRUE)
  j<-which(middle)
  if(length(j)) {
    upper<-pbeta(hi[j],a[j],b[j],log.p=TRUE)
    low<-pbeta(lo[j],a[j],b[j],log.p=TRUE)
    lp[j]<-swm_logdiff(upper,low)
    tail<-j[upper > log(.5)]
    if(length(tail)) lp[tail]<-swm_logdiff(
      pbeta(lo[tail],a[tail],b[tail],lower.tail=FALSE,log.p=TRUE),
      pbeta(hi[tail],a[tail],b[tail],lower.tail=FALSE,log.p=TRUE))
  }
  out<-log1p(-zi)+lp
  j<-which(left)
  if(length(j)) {
    x<-log(zi[j]);y<-out[j];m<-pmax(x,y)
    out[j]<-m+log(exp(x-m)+exp(y-m))
  }
  out
}

swm_rng <- function(mu,phi,grid,zi=0) {
  n<-max(length(mu),length(phi),length(grid),length(zi))
  mu<-rep_len(mu,n);phi<-rep_len(phi,n);grid<-rep_len(grid,n);zi<-rep_len(zi,n)
  z<-floor(rbeta(n,mu*phi,(1-mu)*phi)*grid+.5)
  z[runif(n)<zi]<-0L
  as.integer(pmin(grid,pmax(0,z)))
}

swm_check_folds <- function(d) {
  if(anyNA(d[,c('Subject_ID','CV_fold')])) stop('Missing participant/fold')
  if(any(vapply(split(d$CV_fold,d$Subject_ID),function(x)length(unique(x)),integer(1))!=1L))
    stop('Participant appears in multiple folds')
  invisible(TRUE)
}

swm_logsumexp <- function(x) {
  if(!length(x)||all(x==-Inf))return(-Inf)
  m<-max(x);m+log(sum(exp(x-m)))
}

swm_average_logpmf <- function(logP) {
  if(!is.matrix(logP)||any(dim(logP)<1L)||anyNA(logP)||any(logP==Inf))
    stop('Invalid draw-by-category log probabilities')
  totals<-apply(logP,1,swm_logsumexp)
  if(any(!is.finite(totals))||any(abs(totals)>1e-6))stop('Predictive log mass invalid')
  normalized<-sweep(logP,1,totals,'-')
  apply(normalized,2,swm_logsumexp)-log(nrow(logP))
}

swm_score_pmf <- function(p,k,grid,u,log_p=NULL) {
  if(length(p)!=grid+1 || any(!is.finite(p)) || any(p<0) || abs(sum(p)-1)>1e-7 ||
     length(k)!=1 || k<0 || k>grid || k!=round(k) || u<0 || u>1)
    stop('Invalid PMF or observation for scoring')
  p<-p/sum(p)
  if(is.null(log_p))log_p<-log(p)
  if(length(log_p)!=length(p)||anyNA(log_p)||any(log_p==Inf)||
     any(abs(exp(log_p)-p)>1e-7))stop('Invalid log PMF for scoring')
  total<-swm_logsumexp(log_p)
  if(!is.finite(total)||abs(total)>1e-7)stop('Log PMF is not normalized')
  log_p<-log_p-total
  # Sum each tail directly. Never recover the small upper tail as 1-CDF.
  lower<-swm_logsumexp(c(log_p[seq_len(k)],log_p[k+1]-log(2)))
  after<-if(k<grid)log_p[seq.int(k+2,grid+1)]else numeric()
  upper<-swm_logsumexp(c(after,log_p[k+1]-log(2)))
  z<-if(lower<=upper)qnorm(lower,log.p=TRUE)else qnorm(upper,lower.tail=FALSE,log.p=TRUE)
  cdf<-cumsum(p);cdf[length(cdf)]<-1
  before<-if(k==0)0 else cdf[k];mass<-p[k+1]
  pit<-before+u*mass;mid<-before+.5*mass
  qindex<-function(q)which(cdf>=q)[1]-1L
  out<-list(crps=(12/grid)*sum((cdf[seq_len(grid)]-as.numeric(k<=0:(grid-1)))^2),
    log_score=log_p[k+1],zero_probability=as.integer(log_p[k+1]==-Inf),
    probability_underflow=as.integer(mass==0 && is.finite(log_p[k+1])),p_floor=p[1],
    observed_floor=as.integer(k==0),p_ceiling=p[length(p)],
    pit=pit,mid_p=mid,z_mid=z,log_mid_p_lower=lower,log_mid_p_upper=upper,
    numeric_tail_status=if(is.finite(z))'ok'else'numeric_tail_limit',
    z_rpit=qnorm(pmin(1-1e-6,pmax(1e-6,pit))),
    mean_score=2+12*sum((0:grid)*p)/grid)
  for(level in c(50,80,90,95)) {
    alpha<-(1-level/100)/2;lo<-qindex(alpha);hi<-qindex(1-alpha)
    out[[paste0('coverage',level)]]<-as.integer(k>=lo & k<=hi)
    out[[paste0('below',level)]]<-as.integer(k<lo)
    out[[paste0('above',level)]]<-as.integer(k>hi)
  }
  as.data.frame(out)
}

# Official brms custom-family interface; no package is loaded for base-R tests.
swm_parameters <- function(i,prep) {
  list(mu=brms::get_dpar(prep,'mu',i=i),phi=brms::get_dpar(prep,'phi',i=i),
       zi=brms::get_dpar(prep,'zi',i=i),grid=prep$data$vint1[i])
}
log_lik_swm_floor_beta <- function(i,prep) {
  p<-swm_parameters(i,prep);swm_logpmf(prep$data$Y[i],p$mu,p$phi,p$grid,p$zi)
}
posterior_predict_swm_floor_beta <- function(i,prep,...) {
  p<-swm_parameters(i,prep);swm_rng(p$mu,p$phi,p$grid,p$zi)
}
posterior_epred_swm_floor_beta <- function(prep) {
  vapply(seq_along(prep$data$Y),function(i) {
    p<-swm_parameters(i,prep);res<-rep(0,length(p$mu))
    for(k in seq_len(p$grid)) res<-res+k*exp(swm_logpmf(k,p$mu,p$phi,p$grid,p$zi))
    res
  },numeric(prep$ndraws))
}
swm_family <- function() {
  brms::custom_family('swm_floor_beta',dpars=c('mu','phi','zi'),
    links=c('logit','log','logit'),lb=c(0,0,0),ub=c(1,NA,1),
    type='int',vars='vint1[n]')
}
