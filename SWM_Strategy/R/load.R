swm_root <- dirname(dirname(normalizePath(sys.frame(1)$ofile)))
source(file.path(swm_root,'cantab_design_audit.R'))
for (f in c('state_precision.R','swm_discrete_family.R','swm_diagnostics.R','swm_core.R','swm_folds.R',
            'swm_runtime.R','swm_calibration.R','swm_pipeline.R')) source(file.path(swm_root,'R',f))
