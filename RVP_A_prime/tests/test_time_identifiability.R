args<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
for(f in c('model_spec.R','prepare_raw.R','prepare_data.R','folds.R','fit_models.R'))source(file.path(root,'R',f))
source(file.path(root,'tests/synthetic_input.R'))
tmp<-tempfile('aprime_identifiability_');dir.create(tmp)
write.csv(aprime_test_input(),file.path(tmp,'raw.csv'),row.names=FALSE)
aprime_prepare_raw(file.path(tmp,'raw.csv'),tmp)
p<-prepare_data(file.path(tmp,'aligned_raw.csv'),NULL,tmp)
cfg<-jsonlite::fromJSON(file.path(root,'config.json'))
audits<-list()
for(t in fit_tasks(p,cfg,tmp)) {
 model<-build_model(t$train)
 formula_text<-paste(deparse(model$formula$formula,width.cutoff=500L),collapse=' ')
 stopifnot(grepl('s(Time_since_baseline_1y, k = 4)',formula_text,fixed=TRUE),
   grepl('PlatformConnect_indicator:Time_since_baseline_1y',formula_text,fixed=TRUE),
   !grepl('Platform:Time_since_baseline_1y',formula_text,fixed=TRUE))
 a<-model_design_audit(model,t$train)
 audits[[length(audits)+1L]]<-cbind(kind=t$kind,fold=t$fold,attr(a,'full_audit'))
}
a<-do.call(rbind,audits)
stopifnot(nrow(a)==12L,all(a$full_rank))
unlink(tmp,recursive=TRUE)
cat('PASS: primary common s(t) + Connect indicator:t; all six synthetic normative designs full rank, including penalized bases and phi.\n')
