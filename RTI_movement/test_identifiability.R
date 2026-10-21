# Regression against the actual production model constructors; no MCMC.
args<-commandArgs(TRUE)
here<-dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE)))
scratch<-tempfile('normative_regression_');dir.create(scratch)
source(file.path(here,'core.R'));source(file.path(here,'synthetic_input.R'))
input<-if(length(args))normalizePath(args[1],mustWork=TRUE)else file.path(scratch,'synthetic.csv')
if(!length(args))write.csv(synthetic_normative_input(),input,row.names=FALSE)
e<-normative_context(input,scratch)
setwd(scratch)
spec<-e$outcome_specs[[1]];family<-spec$primary_family
md<-e$prepare_model_data(e$data,spec,family)
make<-function(f,d)brms::make_standata(f$formula,data=d,family=f$family,prior=f$priors)
fixed_time<-e$cantab_time_term
for(role in 'normative') {
  d<-md;builder<-e$build_model_definition
  fixed<-builder(d,family,outcome_prefix=spec$outcome_prefix)
  sd<-make(fixed,d)
  e$cantab_audit_standata(sd,paste(role,'fixed'))
  # Two unrestricted platform slopes remain and no smooth linear time null basis remains.
  stopifnot(sum(grepl('Platform.*:Time_since_baseline_1y',colnames(sd$X)))==2L)
  e$cantab_time_term<-function(smooth,platform_slopes)if(smooth)'s(Time_since_baseline_1y, k = 4)' else 'Time_since_baseline_1y'
  unconstrained<-builder(d,family,outcome_prefix=spec$outcome_prefix)
  old<-make(unconstrained,d)
  rejected<-tryCatch({e$cantab_audit_standata(old,paste(role,'unconstrained'));FALSE},error=function(err)grepl('Rank-deficient',conditionMessage(err)))
  stopifnot(rejected)
  cat('UNCONSTRAINED_TIME_REJECTED',role,'\n')
  e$cantab_time_term<-fixed_time
  e$time_term_type<-'linear'
  linear<-make(builder(d,family,outcome_prefix=spec$outcome_prefix),d)
  e$cantab_audit_standata(linear,paste(role,'linear fallback'))
  stopifnot('Time_since_baseline_1y' %in% colnames(linear$X))
  e$time_term_type<-'smooth';e$include_platform_time_interaction<-FALSE
  no_interaction<-make(builder(d,family,outcome_prefix=spec$outcome_prefix),d)
  e$cantab_audit_standata(no_interaction,paste(role,'no platform slopes'))
  stopifnot(ncol(no_interaction$Xs)>ncol(sd$Xs))
  e$time_term_type<-'linear'
  plain<-make(builder(d,family,outcome_prefix=spec$outcome_prefix),d)
  e$cantab_audit_standata(plain,paste(role,'linear no platform slopes'))
  stopifnot('Time_since_baseline_1y' %in% colnames(plain$X))
  e$time_term_type<-'smooth';e$include_platform_time_interaction<-TRUE
}
cat('IDENTIFIABILITY_REGRESSION_PASS: unconstrained failure, corrected designs, linear fallback and common time effect.\n')
# Unused non-control levels must never produce missing grouped folds.
original<-e$make_grouped_folds(md,5L,1334L)
with_unused<-md
with_unused$Subject_ID<-factor(as.character(md$Subject_ID),levels=c('A_unused',levels(md$Subject_ID),'Z_unused'))
stopifnot(identical(original,e$make_grouped_folds(with_unused,5L,1334L)))
custom<-md;custom$Subject_ID<-factor(as.character(md$Subject_ID),levels=rev(levels(md$Subject_ID)))
custom_unused<-custom;custom_unused$Subject_ID<-factor(as.character(custom$Subject_ID),levels=c('A_unused',levels(custom$Subject_ID),'Z_unused'))
stopifnot(identical(e$make_grouped_folds(custom,5L,1334L),e$make_grouped_folds(custom_unused,5L,1334L)))
cat('PASS: grouped fold assignment ignores unused levels and preserves custom level order.\n')
