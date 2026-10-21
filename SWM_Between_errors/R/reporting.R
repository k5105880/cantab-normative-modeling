report_normative <- function(p,jobs,scores,ctx) {
 out<-file.path(ctx$run,'aggregate');save<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
 save(dplyr::bind_rows(lapply(jobs,`[[`,'diagnostics')),'all_model_diagnostics')
 save(dplyr::bind_rows(lapply(jobs,`[[`,'prior')),'prior_predictive')
 cal<-gaps<-list();counter<-0L
 for(id in names(jobs)){
  d<-jobs[[id]]$oof;age<-cut(d$Age_baseline_years,c(-Inf,40,60,Inf),labels=c('age_le40','age_40to60','age_gt60'))
  strata<-list(all=rep(TRUE,nrow(d)),baseline=d$Practice_bin=='0_baseline',followup=d$Practice_bin!='0_baseline')
  for(a in levels(age))strata[[a]]<-age==a
  for(v in unique(as.character(d$Practice_bin)))strata[[v]]<-d$Practice_bin==v
  for(st in names(strata)){
   g<-d[strata[[st]],,drop=FALSE];if(!nrow(g))next
   for(eq in c(FALSE,TRUE)){
    counter<-counter+1L
    gg<-bootstrap_platform_gap(g,eq,B=ctx$config$bootstrap_reps,seed=ctx$config$seed+counter)
    if(!is.null(gg)){gg$model<-id;gg$stratum<-st;gaps[[length(gaps)+1]]<-gg}
    for(pl in c('RS6','Connect')){
     s<-g[g$Platform==pl,,drop=FALSE];if(!nrow(s))next
     cc<-summarize_calibration(s,eq);cc$model<-id;cc$stratum<-st;cc$Platform<-pl;cal[[length(cal)+1]]<-cc
    }
   }
  }
 }
 save(dplyr::bind_rows(cal),'calibration_by_stratum');save(dplyr::bind_rows(gaps),'platform_gaps_bootstrap')
 cols<-'normative_z'
 counts<-dplyr::bind_rows(lapply(cols,function(n)dplyr::bind_rows(lapply(c('all',ctx$config$groups),function(g){x<-scores[is.finite(scores[[n]])&(g=='all'|scores$Diagnosis==g),,drop=FALSE];data.frame(column=n,group=g,rows=nrow(x),subjects=length(unique(x$Subject_ID)))}))))
 save(counts,'score_availability')
 dictionary<-data.frame(column=cols,role='Normative model',positive_direction='More errors / worse relative to the HC reference distribution',HC='Subject-wise held-out',nonHC='Full HC fit')
 save(dictionary,'score_dictionary')
 counts<-data.frame(rows_raw=nrow(p$raw),rows_state=nrow(p$valid),HC_rows=nrow(p$hc),HC_subjects=length(unique(p$hc$Subject_ID)))
 save(counts,'sample_counts')
 writeLines(c('Normative model for SWM error counts.',
  'Healthy controls use participant-held-out fits; other participants use the full healthy-control fit.',
  'Prior predictive, calibration, coverage and platform-gap diagnostics are included.',
  'Normative z scores use deterministic mid-P, all posterior draws and Gaussian quadrature.',
  'Exposure is set to RS6=4 / Connect=1; confirm that these values match your test protocol.',
  'Bootstrap intervals condition on saved predictions.'),file.path(out,'RESULTS_GUIDE.txt'))
 invisible(NULL)
}
