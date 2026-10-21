# Primary visit/state only: Gaussian new-person integration of conditional
# probabilities, then posterior probability average, then inverse normal.
# Standalone identical copies in both SWM bundles. Longitudinal samplers remain
# separate. Full retained chains, including iteration order, are never thinned.
swm_state_settings<-list(quad_nodes=16L,check_nodes=32L,max_nodes=1024L,
  z_tolerance=1e-4,tail_relative_tolerance=1e-4,draw_probability_tolerance=1e-6,
  posterior_mcse_z_tolerance=.03,chunk_size=20L)
swm_state_gauss_normal<-local({
 cache<-new.env(parent=emptyenv())
 function(n) {
  if(length(n)!=1L||!is.finite(n)||n<2||n!=floor(n))stop('Invalid quadrature order')
  key<-as.character(n);if(exists(key,cache,inherits=FALSE))return(cache[[key]])
  J<-matrix(0,n,n);J[cbind(1:(n-1L),2:n)]<-sqrt(1:(n-1L))
  e<-eigen(J+t(J),symmetric=TRUE);i<-order(e$values);w<-e$vectors[1,i]^2;k<-w>0
  cache[[key]]<-list(nodes=e$values[i][k],weights=w[k]/sum(w[k]));cache[[key]]
 }
})
swm_state_logadd<-function(a,b) {
 m<-pmax(a,b);x<-m+log1p(exp(-abs(a-b)));x[is.infinite(m)&m<0]<--Inf;x
}
swm_state_logdiff<-function(a,b) {
 x<-a+log(-expm1(pmin(0,b-a)));x[is.infinite(a)&a<0]<--Inf;x
}
swm_state_logmean<-function(x) {
 m<-apply(x,2,max);v<-m+log(colMeans(exp(sweep(x,2,m,'-'))));v[is.infinite(m)&m<0]<--Inf;pmin(0,v)
}
swm_z_from_log_probabilities<-function(lo,up)ifelse(lo<=up,qnorm(lo,log.p=TRUE),qnorm(up,lower.tail=FALSE,log.p=TRUE))

swm_state_re_variance<-function(covariance,newdata) {
 d<-dim(covariance);cn<-dimnames(covariance)[[2]]
 if(length(d)!=3L||d[1]<1L||d[2]!=d[3]||any(!is.finite(covariance))||
   is.null(cn)||!identical(cn,dimnames(covariance)[[3]])||
   !(setequal(cn,'Intercept')||setequal(cn,c('Intercept','Time_since_baseline_1y'))))stop('Unsupported subject covariance')
 tol<-1e-10*max(1,abs(covariance))
 if(max(abs(covariance-aperm(covariance,c(1,3,2))))>tol)stop('Asymmetric subject covariance')
 if(any(covariance[,1,1]< -tol))stop('Nonpositive subject covariance')
 if(d[2]==2L&&any(covariance[,2,2]< -tol|covariance[,1,1]*covariance[,2,2]-covariance[,1,2]^2 < -tol))stop('Nonpositive subject covariance')
 z<-matrix(1,nrow(newdata),length(cn),dimnames=list(NULL,cn))
 if('Time_since_baseline_1y'%in%cn) {
  t<-newdata$Time_since_baseline_1y
  if(length(t)!=nrow(newdata)||any(!is.finite(t)))stop('Finite random-slope time required')
  z[,'Time_since_baseline_1y']<-t
 }
 v<-matrix(0,d[1],nrow(newdata))
 for(j in seq_along(cn))for(k in seq_along(cn))v<-v+outer(covariance[,j,k],z[,j]*z[,k])
 if(any(!is.finite(v))||any(v< -tol))stop('Invalid new-person variance')
 pmax(v,0)
}

# Likelihood probability below, at, and above the observed category. Each tail
# is evaluated directly in log space (never as one minus a rounded CDF).
swm_state_log_probabilities<-function(response,location,aux,zero_eta,family_id,grid=NULL) {
 nr<-nrow(location);nc<-ncol(location);y<-matrix(response,nr,nc,byrow=TRUE)
 lzero<-plogis(zero_eta,log.p=TRUE);lpositive<-plogis(-zero_eta,log.p=TRUE)
 if(family_id=='hdg') {
  mu<-exp(location);sigma<-aux
  den<-pnorm(.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)
  log_interval<-function(a,b) {
   za<-(a-mu)/sigma;zb<-(b-mu)/sigma
   hi<-ifelse(za>=0,pnorm(za,lower.tail=FALSE,log.p=TRUE),pnorm(zb,log.p=TRUE))
   lo<-ifelse(za>=0,pnorm(zb,lower.tail=FALSE,log.p=TRUE),pnorm(za,log.p=TRUE))
   out<-swm_state_logdiff(hi,lo);w<-(b-a)/sigma;mid<-((a+b)/2-mu)/sigma
   narrow<-w*(1+abs(mid))<.01&w>0
   out[narrow]<-dnorm(mid[narrow],log=TRUE)+log(w[narrow])+log1p((mid[narrow]^2-1)*w[narrow]^2/24+(mid[narrow]^4-6*mid[narrow]^2+3)*w[narrow]^4/1920)
   out
  }
  mass<-lpositive+pmin(0,log_interval(pmax(.5,y-.5),y+.5)-den)
  lower<-swm_state_logadd(lzero,lpositive+pmin(0,log_interval(.5,pmax(.5,y-.5))-den))
  upper<-lpositive+pmin(0,pnorm(y+.5,mu,sigma,lower.tail=FALSE,log.p=TRUE)-den)
  mass[y==0]<-lzero[y==0];lower[y==0]<--Inf;upper[y==0]<-lpositive[y==0]
 } else if(family_id=='beta') {
  g<-matrix(grid,nr,nc,byrow=TRUE);lo<-pmax(0,(y-.5)/g);hi<-pmin(1,(y+.5)/g)
  # Separate logit tails preserve both shape parameters when plogis rounds to1.
  a<-exp(log(aux)+plogis(location,log.p=TRUE));b<-exp(log(aux)+plogis(-location,log.p=TRUE))
  cl<-matrix(pbeta(lo,a,b,log.p=TRUE),nr,nc);ch<-matrix(pbeta(hi,a,b,log.p=TRUE),nr,nc)
  sl<-matrix(pbeta(lo,a,b,lower.tail=FALSE,log.p=TRUE),nr,nc);sh<-matrix(pbeta(hi,a,b,lower.tail=FALSE,log.p=TRUE),nr,nc)
  lp<-swm_state_logdiff(ch,cl);right<-ch>log(.5);lp[right]<-swm_state_logdiff(sl[right],sh[right])
  mass<-lpositive+lp;lower<-swm_state_logadd(lzero,lpositive+cl);upper<-lpositive+sh
  mass[y==0]<-swm_state_logadd(lzero[y==0],mass[y==0]);lower[y==0]<--Inf;upper[y==g]<--Inf
 } else stop('Unknown SWM state family')
 list(lower=lower,mass=mass,upper=upper)
}

swm_state_components<-function(response,eta,aux,zero_eta,variance,family_id,grid,nodes) {
 g<-if(all(variance==0))list(nodes=0,weights=1)else swm_state_gauss_normal(nodes);sd<-sqrt(variance)
 logs<-list(lower=matrix(-Inf,nrow(eta),ncol(eta)),mass=matrix(-Inf,nrow(eta),ncol(eta)),upper=matrix(-Inf,nrow(eta),ncol(eta)))
 for(k in seq_along(g$nodes)) {
  p<-swm_state_log_probabilities(response,eta+sd*g$nodes[k],aux,zero_eta,family_id,grid)
  for(n in names(logs))logs[[n]]<-swm_state_logadd(logs[[n]],p[[n]]+log(g$weights[k]))
 }
 if(any(vapply(logs,function(x)anyNA(x)||any(x==Inf),logical(1))))stop('Nonfinite integrated likelihood')
 logs$mid_lower<-swm_state_logadd(logs$lower,logs$mass-log(2));logs$mid_upper<-swm_state_logadd(logs$upper,logs$mass-log(2))
 logs
}

swm_state_posterior_precision<-function(log_tail,metadata,log_mean_tail,z) {
 if(!is.matrix(log_tail)||nrow(metadata)!=nrow(log_tail)||!all(c('.chain','.iteration','.draw')%in%names(metadata))||
    anyNA(metadata)||anyDuplicated(metadata$.draw)||anyDuplicated(metadata[c('.chain','.iteration')]))stop('Invalid state chain metadata')
 chains<-sort(unique(metadata$.chain));ids<-lapply(chains,function(k){i<-which(metadata$.chain==k);i[order(metadata$.iteration[i])]})
 available<-length(unique(lengths(ids)))==1L&&all(vapply(ids,function(i)length(i)>=4&&all(diff(metadata$.iteration[i])==1),logical(1)))
 do.call(rbind,lapply(seq_len(ncol(log_tail)),function(j) {
  scale<-max(log_tail[,j]);x<-exp(log_tail[,j]-scale);se<-ess<-rh<-NA_real_;method<-'unavailable_short_or_incomplete_chains'
  if(available&&all(is.finite(x))&&is.finite(z[j])) {
   if(diff(range(x))==0){se<-0;ess<-length(x);method<-'constant_draw_contribution'}else{
    a<-vapply(ids,function(i)x[i],numeric(length(ids[[1]])))
    se<-as.numeric(posterior::mcse_mean(a));ess<-as.numeric(posterior::ess_mean(a));rh<-as.numeric(posterior::rhat(a));method<-'chain_aware_spectral_ess_delta_method'
   }
  }
  rel<-se/mean(x);zse<-rel*exp(log_mean_tail[j]-dnorm(z[j],log=TRUE))
  data.frame(posterior_draws=nrow(log_tail),posterior_chains=length(chains),posterior_mcse_p=se*exp(scale),
    posterior_relative_mcse_tail=rel,posterior_mcse_z=zse,posterior_ess_mean=ess,posterior_rhat=rh,
    posterior_precision_ok=is.finite(zse)&&zse<=swm_state_settings$posterior_mcse_z_tolerance,posterior_mcse_method=method)
 }))
}

swm_state_checked<-function(response,eta,aux,zero_eta,variance,metadata,family_id,grid=NULL,
 quad_nodes=swm_state_settings$quad_nodes,check_nodes=swm_state_settings$check_nodes,max_nodes=swm_state_settings$max_nodes) {
 if(!is.matrix(eta)||any(dim(eta)<1)||any(vapply(list(aux,zero_eta,variance),function(x)!is.matrix(x)||!identical(dim(x),dim(eta)),logical(1)))||
    length(response)!=ncol(eta)||any(!is.finite(c(response,eta,aux,zero_eta,variance)))||any(aux<=0)||any(variance<0)||
    any(response<0|response!=round(response)))stop('Invalid state integration inputs')
 if(family_id=='beta'&&(length(grid)!=length(response)||any(!is.finite(grid)|grid<1|grid!=round(grid)|response>grid)))stop('Invalid state category grid')
 if(any(!is.finite(c(quad_nodes,check_nodes,max_nodes)))||any(c(quad_nodes,check_nodes,max_nodes)!=floor(c(quad_nodes,check_nodes,max_nodes)))||quad_nodes<2||quad_nodes>=check_nodes||check_nodes>max_nodes)stop('Invalid adaptive quadrature orders')
 calc<-function(ix,n)swm_state_components(response[ix],eta[,ix,drop=FALSE],aux[,ix,drop=FALSE],zero_eta[,ix,drop=FALSE],variance[,ix,drop=FALSE],family_id,if(is.null(grid))NULL else grid[ix],n)
 compare<-function(a,b,ix,n) {
  ratio<-if(family_id=='hdg')sqrt(variance[,ix,drop=FALSE])*matrix(pmax(response[ix],.5),nrow(eta),length(ix),byrow=TRUE)/aux[,ix,drop=FALSE]else sqrt(variance[,ix,drop=FALSE]*aux[,ix,drop=FALSE])*.5
  if(family_id=='hdg')ratio[,response[ix]==0]<-0
  resolution<-apply(ratio,2,max)*pi/sqrt(n)
  al<-swm_state_logmean(a$mid_lower);au<-swm_state_logmean(a$mid_upper);bl<-swm_state_logmean(b$mid_lower);bu<-swm_state_logmean(b$mid_upper)
  zd<-abs(swm_z_from_log_probabilities(al,au)-swm_z_from_log_probabilities(bl,bu));td<-pmax(abs(expm1(al-bl)),abs(expm1(au-bu)))
  dd<-pmax(apply(abs(exp(a$mid_lower)-exp(b$mid_lower)),2,max),apply(abs(exp(a$mid_upper)-exp(b$mid_upper)),2,max))
  norm<-apply(abs(exp(b$lower)+exp(b$mass)+exp(b$upper)-1),2,max)
  data.frame(quad_z_diff=zd,quad_relative_tail_diff=td,quad_max_draw_probability_diff=dd,quad_step_to_scale=resolution,
    numerical_ok=is.finite(resolution)&resolution<=1&is.finite(zd)&zd<=swm_state_settings$z_tolerance&is.finite(td)&td<=swm_state_settings$tail_relative_tolerance&is.finite(dd)&dd<=swm_state_settings$draw_probability_tolerance&norm<1e-9)
 }
 ix<-seq_len(ncol(eta));a<-calc(ix,quad_nodes);b<-calc(ix,check_nodes);diag<-compare(a,b,ix,check_nodes);used<-rep(check_nodes,length(ix));n<-check_nodes
 while(any(!diag$numerical_ok)&&n<max_nodes) {
  ix<-which(!diag$numerical_ok);n<-min(2L*n,max_nodes);finer<-calc(ix,n)
  diag[ix,]<-compare(lapply(b,function(x)x[,ix,drop=FALSE]),finer,ix,n)
  for(key in names(b))b[[key]][,ix]<-finer[[key]]
  used[ix]<-n
 }
 ll<-swm_state_logmean(b$lower);lm<-swm_state_logmean(b$mass);lu<-swm_state_logmean(b$upper)
 ml<-swm_state_logmean(b$mid_lower);mu<-swm_state_logmean(b$mid_upper);z<-swm_z_from_log_probabilities(ml,mu)
 tail<-b$mid_lower;i<-ml>mu;tail[,i]<-b$mid_upper[,i,drop=FALSE]
 precision<-swm_state_posterior_precision(tail,metadata,pmin(ml,mu),z)
 cbind(data.frame(lower=exp(ll),mass=exp(lm),upper=exp(lu),log_lower=ll,log_mass=lm,log_upper=lu,
  mid_p=exp(ml),z=z,pred_zero=colMeans(plogis(zero_eta)),log_mid_p_lower=ml,log_mid_p_upper=mu,quad_nodes_used=used,
  scoring_method='conditional_probability_gaussian_quadrature_all_posterior_draws'),diag,precision)
}

swm_state_validate_fit<-function(fit,family_id) {
 fam<-if(family_id=='hdg')'hurdle_discrete_gaussian'else if(family_id=='beta')'swm_floor_beta'else stop('Unsupported state family')
 if(!identical(fit$family$name,fam)||!identical(fit$family$link,if(family_id=='hdg')'log'else'logit')||
  !identical(fit$family[[if(family_id=='hdg')'link_sigma'else'link_phi']],'log')||
  !identical(fit$family[[if(family_id=='hdg')'link_hu'else'link_zi']],'logit'))stop('State family/link mismatch')
 tt<-brms::brmsterms(fit$formula);aux<-if(family_id=='hdg')'sigma'else'phi';zero<-if(family_id=='hdg')'hu'else'zi'
 if(length(tt$nlpars)||length(fit$formula$pfix)||!setequal(names(tt$dpars),c('mu',aux,zero)))stop('Unsupported state distributional structure')
 for(dp in c(aux,zero))if(!is.null(tt$dpars[[dp]][['re',exact=TRUE]]))stop('Random auxiliary effects are not supported by one-dimensional state integration')
 re<-tt$dpars$mu$re
 if(is.null(re)||nrow(re)!=1||!identical(as.character(re$group),'Subject_ID'))stop('One Gaussian Subject_ID term required')
 gc<-re$gcall[[1]]
 if(!identical(gc$dist,'gaussian')||any(vapply(c('by','pw','cov','type'),function(x)!identical(gc[[x]],''),logical(1))))stop('Unsupported subject effects')
 if(family_id=='beta'&&!identical(names(tt$dpars$mu$adforms),'vint'))stop('Expected the Strategy category grid')
 # Frozen HDG exposure offsets are already included by brms get_dpar below.
 allowed<-c(if(family_id=='beta')'adforms',if(family_id=='hdg')'offset','formula','fe','re','sm','allvars','family','dpar','resp','respform','transform')
 if(any(vapply(tt$dpars,function(dp)!all(names(dp)%in%allowed),logical(1))))stop('Unsupported special state predictor')
 invisible(TRUE)
}

swm_state_predict<-function(fit,newdata,family_id,chunk_size=swm_state_settings$chunk_size) {
 swm_state_validate_fit(fit,family_id)
 if(!is.data.frame(newdata)||!nrow(newdata)||anyNA(newdata$row_id_global)||anyDuplicated(newdata$row_id_global))stop('State data require unique row IDs')
 if(length(chunk_size)!=1||!is.finite(chunk_size)||chunk_size<1||chunk_size!=floor(chunk_size))stop('Invalid state chunk size')
 metadata<-as.data.frame(posterior::as_draws_df(fit,variable='lp__'))[c('.chain','.iteration','.draw')]
 if(!identical(as.integer(metadata$.draw),seq_len(nrow(metadata))))stop('Unexpected posterior draw order')
 vc<-brms::VarCorr(fit,summary=FALSE)
 if(!identical(names(vc),'Subject_ID'))stop('Unsupported covariance blocks')
 covariance<-vc$Subject_ID$cov
 if(is.null(covariance)) {
  sd<-vc$Subject_ID$sd
  if(!is.matrix(sd)||ncol(sd)!=1||!identical(colnames(sd),'Intercept'))stop('Unsupported scalar covariance')
  covariance<-array(sd^2,c(nrow(sd),1L,1L),dimnames=list(NULL,'Intercept','Intercept'))
 }
 if(dim(covariance)[1]!=nrow(metadata))stop('Posterior covariance/draw misalignment')
 variance<-swm_state_re_variance(covariance,newdata)
 chunks<-split(seq_len(nrow(newdata)),ceiling(seq_len(nrow(newdata))/chunk_size))
 evaluate<-function(ix) {
  nd<-newdata[ix,,drop=FALSE];model_nd<-if(nrow(nd)==1L)nd[c(1L,1L),,drop=FALSE]else nd
  # One expensive brms/mgcv adapter per chunk, never per quadrature node.
  prep<-brms::prepare_predictions(fit,newdata=model_nd,re_formula=NA,allow_new_levels=TRUE,draw_ids=seq_len(nrow(metadata)))
  gd<-getFromNamespace('get_dpar','brms')
  get<-function(dp,inv=FALSE) {
   x<-gd(prep,dp,inv_link=inv)
   if(is.null(dim(x)))x<-matrix(x,nrow(metadata),nrow(model_nd))
   if(nrow(x)!=nrow(metadata))stop('Prediction/draw misalignment')
   x[,seq_len(nrow(nd)),drop=FALSE]
  }
  result<-swm_state_checked(if(family_id=='hdg')nd$NM_response else nd$index,
   get('mu'),get(if(family_id=='hdg')'sigma'else'phi',TRUE),get(if(family_id=='hdg')'hu'else'zi'),
   variance[,ix,drop=FALSE],metadata,family_id,if(family_id=='beta')nd$grid else NULL)
  cbind(row_id_global=nd$row_id_global,result)
 }
 budget<-suppressWarnings(as.numeric(Sys.getenv('CANTAB_JOB_CORES','1')))
 if(length(budget)!=1L||!is.finite(budget)||budget!=floor(budget)||budget<1)stop('CANTAB_JOB_CORES must be a positive integer')
 workers<-min(4L,budget,length(chunks));out<-if(workers>1L&&.Platform$OS.type!='windows')parallel::mclapply(chunks,evaluate,mc.cores=workers,mc.set.seed=FALSE)else lapply(chunks,evaluate)
 if(any(vapply(out,inherits,logical(1),'try-error')))stop('State chunk failed: ',paste(out[vapply(out,inherits,logical(1),'try-error')],collapse='; '))
 out<-do.call(rbind,out);rownames(out)<-NULL;out
}

swm_state_write_diagnostics<-function(result,path) {
 if(!is.character(path)||length(path)!=1||!nzchar(path))stop('State diagnostic path required')
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 write.csv(result,path,row.names=FALSE,na='')
 if(any(!result$numerical_ok|!result$posterior_precision_ok|!is.finite(result$z)))stop(
  'Primary state precision unresolved; diagnostics saved to ',path,
  '. Require converged quadrature and chain-aware MCSE(z) <= ',swm_state_settings$posterior_mcse_z_tolerance,
  '; scores withheld. Increase posterior sampling for posterior error; quadrature error is separate.',call.=FALSE)
 invisible(result)
}

# Even an exceptional numerical/adapter failure is persisted before refusing
# release, with the affected row IDs. This does not hide or substitute a score.
swm_state_predict_release<-function(fit,newdata,family_id,diagnostics_file) {
 result<-tryCatch(swm_state_predict(fit,newdata,family_id),error=function(e)
  data.frame(row_id_global=newdata$row_id_global,z=NA_real_,numerical_ok=FALSE,
    posterior_precision_ok=FALSE,state_error=conditionMessage(e)))
 swm_state_write_diagnostics(result,diagnostics_file)
 result
}
