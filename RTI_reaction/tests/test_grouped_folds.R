# Participant IDs may retain factor levels from people outside the training set.
args<-gsub('~+~', ' ', commandArgs(FALSE), fixed=TRUE)
here<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
local({
  e<-new.env(parent=asNamespace('dplyr'))
  e$cv_balance_include_fsiq<-FALSE
  sys.source(file.path(here,'grouped_folds.R'),envir=e)
  d<-expand.grid(Timepoint=1:4,subject=seq(2,20,2))
  d$Subject_ID<-factor(sprintf('%03d',d$subject),levels=sprintf('%03d',1:20))
  d$Platform<-ifelse(d$subject<=10,'RS6','Connect')
  d$Gender<-ifelse(d$subject%%4,'F','M')
  d$Age_baseline_10y<-(30+d$subject-50)/10
  d$FSIQ_15<-(90+d$subject-100)/15
  a<-e$make_grouped_folds(d,5L,1334L,TRUE,8L)
  b<-e$make_grouped_folds(droplevels(d),5L,1334L,TRUE,8L)
  stopifnot(identical(a,b),!anyNA(a$folds),setequal(a$folds,1:5),
    all(vapply(split(a$folds,d$Subject_ID,drop=TRUE),
      function(x)length(unique(x))==1L,logical(1))))
  d$Subject_ID<-factor(as.character(d$Subject_ID),levels=rev(levels(d$Subject_ID)))
  custom<-e$make_grouped_folds(d,5L,1334L,FALSE,1L)
  set.seed(1334L)
  expected<-loo::kfold_split_grouped(5L,droplevels(d$Subject_ID))
  stopifnot(identical(custom$folds,as.integer(expected)))
  cat('PASS: unused participant factor levels cannot cause missing folds; custom level order, partition, balance, and audit are preserved.\n')
})
