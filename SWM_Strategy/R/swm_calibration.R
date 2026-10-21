swm_calibration<-function(pred) {
  r<-pred$rows;groups<-list(ALL=rep(TRUE,nrow(r)))
  for(p in c('RS6','Connect')) {
    groups[[p]]<-r$Platform==p
    for(a in c('under40','40_to_under65','65plus'))groups[[paste(p,a,sep=':')]]<-r$Platform==p & r$Age_band==a
  }
  do.call(rbind,lapply(names(groups),function(g) {
    x<-r[groups[[g]],]
    data.frame(model=pred$model,group=g,n=nrow(x),subjects=length(unique(x$Subject_ID)),
      mean_z=mean(x$z_mid),sd_z=sd(x$z_mid),
      sampled_randomized_pit_mean_z=mean(x$ez),sampled_randomized_pit_sd_z=sqrt(mean(x$ez2)-mean(x$ez)^2),
      observed_floor=mean(x$observed_floor),predicted_floor=mean(x$p_floor),
      crps=mean(x$crps),coverage95=mean(x$coverage95))
  }))
}
