# Hurdle discretized lower-truncated Gaussian, integer support 0,1,2,...
# Standalone numerical functions require no brms installation.
hdg_log_interval <- function(a,b,mu,sigma){
 za<-(a-mu)/sigma;zb<-(b-mu)/sigma
 right<-za>=0
 hi<-ifelse(right,pnorm(za,lower.tail=FALSE,log.p=TRUE),pnorm(zb,log.p=TRUE))
 lo<-ifelse(right,pnorm(zb,lower.tail=FALSE,log.p=TRUE),pnorm(za,log.p=TRUE))
 ans<-hi+log(-expm1(pmin(0,lo-hi)))
 width<-(b-a)/sigma;mid<-((a+b)/2-mu)/sigma
 narrow<-width*(1+abs(mid))<.01 & width>0
 # Fourth-order local integral; avoids subtracting indistinguishable CDFs.
 if(any(narrow)){
  w<-width[narrow];z<-mid[narrow]
  ans[narrow]<-dnorm(z,log=TRUE)+log(w)+log1p((z*z-1)*w*w/24+(z^4-6*z*z+3)*w^4/1920)
 }
 ans
}
hdg_log_add <- function(a,b){m<-pmax(a,b);z<-m+log1p(exp(pmin(a,b)-m));z[is.infinite(m)&m<0]<--Inf;z}
hdg_probabilities <- function(y,mu,sigma,hu){
 mu<-as.matrix(mu);sigma<-as.matrix(sigma);hu<-as.matrix(hu)
 if((is.matrix(y)&&!identical(dim(y),dim(mu)))||(!is.matrix(y)&&length(y)!=ncol(mu))||!identical(dim(mu),dim(sigma))||!identical(dim(mu),dim(hu)))stop('HDG dimensions disagree')
 if(any(!is.finite(y)|y<0|y!=round(y))||any(!is.finite(mu)|mu<=0)||any(!is.finite(sigma)|sigma<=0)||any(!is.finite(hu)|hu<0|hu>1))stop('Invalid HDG support or parameters')
 yy<-if(is.matrix(y))y else matrix(rep(y,each=nrow(mu)),nrow(mu));zero<-yy==0
 den<-pnorm(.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)
 lp<-hdg_log_interval(pmax(.5,yy-.5),yy+.5,mu,sigma)-den
 mass<-(1-hu)*exp(pmin(0,lp))
 prev<-pnorm(pmax(.5,yy-.5),mu,sigma,lower.tail=FALSE,log.p=TRUE)
 last<-pnorm(yy+.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)
 lower<-hu+(1-hu)*(-expm1(pmin(0,prev-den)))
 upper<-(1-hu)*exp(pmin(0,last-den))
 mass[zero]<-hu[zero];lower[zero]<-0;upper[zero]<-1-hu[zero]
 logmass<-log1p(-hu)+pmin(0,lp);logupper<-log1p(-hu)+pmin(0,last-den)
 lpos<-hdg_log_interval(.5,pmax(.5,yy-.5),mu,sigma)-den
 loglower<-hdg_log_add(log(hu),log1p(-hu)+pmin(0,lpos))
 lower<-exp(loglower);lower[zero]<-0
 logmass[zero]<-log(hu[zero]);loglower[zero]<--Inf;logupper[zero]<-log1p(-hu[zero])
 list(lower=lower,mass=mass,upper=upper,log_lower=loglower,log_mass=logmass,log_upper=logupper)
}
hdg_rng <- function(mu,sigma,hu){
 stopifnot(length(mu)==length(sigma),length(mu)==length(hu),all(is.finite(mu)&mu>0),all(is.finite(sigma)&sigma>0),all(is.finite(hu)&hu>=0&hu<=1))
 lt<-pnorm(.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)+log(runif(length(mu)))
 x<-mu+sigma*qnorm(lt,lower.tail=FALSE,log.p=TRUE)
 y<-pmax(1,floor(x+.5));y[runif(length(mu))<hu]<-0;y
}
HDG_STAN <- '
  real hdg_log_tail(real z) {
    real q;
    real az = abs(z);
    if (az > 20) {
      real u = inv_square(az);
      real series = 1 + u * (-1 + u * (3 + u * (-15 + u * (105 + u * (-945 + u * (10395 + u * (-135135 + u * 2027025)))))));
      q = -0.5 * square(az) - log(az) - 0.5 * log(2 * pi()) + log(series);
    } else {
      // std_normal_lccdf underflows already around z=8 in current Stan Math.
      q = log(erfc(az / sqrt(2.0))) - log(2.0);
    }
    if (z >= 0) return q;
    return log1m_exp(q);
  }
  real hdg_log_interval(real a, real b, real mu, real sigma) {
    real za = (a - mu) / sigma;
    real zb = (b - mu) / sigma;
    real w = (b - a) / sigma;
    real z = (0.5 * (a + b) - mu) / sigma;
    if (w * (1 + abs(z)) < 0.01)
      return std_normal_lpdf(z) + log(w) + log1p((square(z)-1)*square(w)/24 + (pow(z,4)-6*square(z)+3)*pow(w,4)/1920);
    if (za >= 0) return log_diff_exp(hdg_log_tail(za), hdg_log_tail(zb));
    return log_diff_exp(hdg_log_tail(-zb), hdg_log_tail(-za));
  }
  real hurdle_discrete_gaussian_lpmf(int y, real mu, real sigma, real hu) {
    if (y == 0) return bernoulli_lpmf(1 | hu);
    if (y < 0) return negative_infinity();
    return bernoulli_lpmf(0 | hu) + hdg_log_interval(y-0.5,y+0.5,mu,sigma) - hdg_log_tail((0.5-mu)/sigma);
  }
'
# brms stores these functions on the custom family, keeping cached fits portable.
log_lik_hurdle_discrete_gaussian <- function(i,prep){
 gd<-getFromNamespace('get_dpar','brms')
 mu<-gd(prep,'mu',i=i);sigma<-gd(prep,'sigma',i=i);hu<-gd(prep,'hu',i=i);y<-prep$data$Y[i]
 if(y==0)return(log(hu))
 log1p(-hu)+hdg_log_interval(y-.5,y+.5,mu,sigma)-pnorm(.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)
}
posterior_predict_hurdle_discrete_gaussian <- function(i,prep,...){
 gd<-getFromNamespace('get_dpar','brms')
 hdg_rng(gd(prep,'mu',i=i),gd(prep,'sigma',i=i),gd(prep,'hu',i=i))
}
hdg_family <- function()brms::custom_family('hurdle_discrete_gaussian',dpars=c('mu','sigma','hu'),links=c('log','log','logit'),
 lb=c(0,0,0),ub=c(NA,NA,1),type='int',
 log_lik=log_lik_hurdle_discrete_gaussian,posterior_predict=posterior_predict_hurdle_discrete_gaussian)
hdg_stanvars <- function()brms::stanvar(scode=HDG_STAN,block='functions')
