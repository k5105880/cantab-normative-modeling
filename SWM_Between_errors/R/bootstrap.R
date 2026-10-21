# Load definitions into a dedicated environment without attaching packages or
# changing the caller's options, environment variables, or R library paths.
load_runtime <- function(root, envir=new.env(parent=parent.frame())) {
 packages<-c('brms','posterior','dplyr','mgcv','cmdstanr','digest','loo','zip')
 missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '),'. Complete the shared environment/setup.R setup first.')
 # Match the namespace precedence of brms, posterior, then dplyr, while keeping
 # these imports local to the analysis runtime.
 for(package in c('brms','posterior','dplyr')) {
  for(name in getNamespaceExports(package))
   assign(name,getExportedValue(package,name),envir=envir)
 }
 files<-c('cantab_design_audit.R','config.R','R/state_precision.R',
  file.path('vendor',c('count_distribution.R','folds.R','data_diagnostics.R','probability_scores.R','normative_prediction.R','calibration.R')),
  'R/data.R','vendor/common_model_definitions.R','vendor/final_model_definition.R',
  'R/models.R','R/reporting.R','R/pipeline.R')
 for(file in files)sys.source(file.path(root,file),envir=envir)
 helper<-new.env(parent=envir)
 sys.source(file.path(dirname(root),'R','input_data.R'),envir=helper)
 envir$cantab_output_directory<-helper$cantab_output_directory
 envir
}
