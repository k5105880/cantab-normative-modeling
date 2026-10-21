# Synthetic design regression; no private input or model fitting is needed.
args<-gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE)
here<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
source(file.path(here,'rti_core.R'))
local({
  root<-tempfile('reaction_design_');dir.create(root)
  on.exit(unlink(root,recursive=TRUE))
  d<-expand.grid(Timepoint=1:4,subject=1:40)
  d$Subject_ID<-sprintf('%03d',d$subject)
  d$Diagnosis<-'HC';d$Gender<-ifelse(d$subject%%2,'F','M')
  d$Platform<-ifelse(d$subject<=20,'RS6','Connect')
  d$FSIQ<-85+d$subject%%17
  d$Age_baseline_years<-25+d$subject
  d$Time_since_baseline_years<-(d$Timepoint-1)*(0.7+(d$subject%%7)/10)
  d[['RTI_Median_five-choice_reaction_time']]<-300+d$subject+d$Timepoint
  input<-file.path(root,'input.csv');write.csv(d,input,row.names=FALSE)
  folds<-file.path(root,'folds.csv')
  write.csv(data.frame(Subject_ID=sprintf('%03d',1:40),CV_fold=rep(1:5,8)),
            folds,row.names=FALSE)
  data<-rti_prepare(input,folds)
  stopifnot(identical(names(data),c('all','population')))
  def<-rti_definition()
  for(fold in 0:5) {
    hc<-data$population
    if(fold)hc<-hc[hc$CV_fold!=fold,,drop=FALSE]
    sd<-rti_standata(hc,def)
    a<-cantab_audit_standata(sd,paste('normative fold',fold))
    stopifnot(all(a$full_rank))
  }
  duplicated<-def
  duplicated$formula$formula<-as.formula(gsub(
    's(Time_since_baseline_1y, k = 4, bs = "tp", m = c(2, 0))',
    's(Time_since_baseline_1y, k=4)',
    paste(deparse(def$formula$formula,width.cutoff=500L),collapse=' '),fixed=TRUE))
  stopifnot(!identical(deparse(duplicated$formula$formula),deparse(def$formula$formula)))
  sd<-brms::make_standata(duplicated$formula,data=data$population,
                        family=duplicated$family,prior=duplicated$priors)
  error<-tryCatch(cantab_audit_standata(sd,'duplicated linear time term'),error=identity)
  stopifnot(inherits(error,'error'),grepl('Rank-deficient population design',conditionMessage(error)))
  cat('PASS: six normative training designs are identifiable; duplicated time slope is rejected.\n')
})
