# Exercise reference selection and snapshots before dispatching any MCMC.
args<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
local({
  tmp<-tempfile('aprime references ');dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE))
  e<-new.env(parent=globalenv())
  for(file in c('prepare_raw.R','prepare_data.R','folds.R','portable_pipeline.R'))
    sys.source(file.path(root,'R',file),envir=e)
  sys.source(file.path(root,'tests/synthetic_input.R'),envir=e)
  # The real launcher writes its source snapshot before this external-process boundary.
  e$system2<-function(command,args,...)stop('TEST_BEFORE_STAGES',call.=FALSE)
  source(file.path(root,'run.R'),local=e)
  raw<-e$aprime_test_input()
  raw$SWM_Between_errors<-4L;raw$SWM_Strategy<-5L
  raw$RVP_Median_latency<-400
  raw[['RTI_Median_five-choice_reaction_time']]<-320
  raw[['RTI_Median_five-choice_movement_time']]<-250
  raw$Diagnosis2<-rep(c('subtype A','','subtype B','unrecorded'),length.out=nrow(raw))
  raw$source_note<-rep(c('comma, retained','quote " retained','  padded  ','note'),length.out=nrow(raw))
  input<-file.path(tmp,'raw.csv');write.csv(raw,input,row.names=FALSE)
  refs<-file.path(tmp,'references');dir.create(file.path(refs,'reference'),recursive=TRUE)
  cohort<-unique(raw[c('Subject_ID','Diagnosis')])
  cohort<-cohort[cohort$Diagnosis!='MDD',]
  write.csv(cohort,file.path(refs,'diagnosis_reference.csv'),row.names=FALSE)
  fold_map<-data.frame(Subject_ID=sprintf('synthetic_%03d',1:60),CV_fold=rep(5:1,12))
  root_folds<-file.path(refs,'primary_subject_folds.csv')
  nested_folds<-file.path(refs,'reference/primary_subject_folds.csv')
  write.csv(fold_map,root_folds,row.names=FALSE)
  # Historical files have no authority over cohort labels, fold choice or identity.
  historical<-raw[c('Subject_ID','Timepoint','Platform','Diagnosis')]
  historical$Diagnosis<-'MDD'
  write.csv(historical,file.path(refs,'reference/diagnosis_key.csv'),row.names=FALSE)
  writeLines('unused obsolete provenance',file.path(refs,'provenance.json'))
  stopped<-function(expr,pattern) {
    err<-tryCatch(force(expr),error=identity)
    if(!inherits(err,'error'))stop('Expected failure containing: ',pattern)
    if(!grepl(pattern,conditionMessage(err),fixed=TRUE))
      stop('Expected ',pattern,'; got: ',conditionMessage(err))
  }
  root_run<-file.path(tmp,'root folds run')
  stopped(e$run_outcome(input,root_run,'check',reference_dir=refs),'TEST_BEFORE_STAGES')
  snapshot<-file.path(root_run,'reference/primary_subject_folds.csv')
  stopifnot(file.exists(snapshot),identical(e$aprime_hash(snapshot),e$aprime_hash(root_folds)),
    identical(list.files(file.path(root_run,'reference')),'primary_subject_folds.csv'))
  aligned<-e$aprime_read_cells(file.path(root_run,'input/aligned_raw.csv'))
  expected<-e$aprime_read_cells(input)
  expected<-expected[expected$Diagnosis!='MDD',,drop=FALSE];rownames(expected)<-NULL
  stopifnot(identical(aligned,expected),
    setequal(list.files(file.path(root_run,'input')),c('aligned_raw.csv','source_manifest.json')))
  p<-e$prepare_data(file.path(root_run,'input/aligned_raw.csv'),dirname(snapshot),file.path(tmp,'prepared'))
  stopifnot(all(p$folds$CV_fold==fold_map$CV_fold[match(p$folds$Subject_ID,fold_map$Subject_ID)]))
  provenance<-jsonlite::fromJSON(file.path(tmp,'prepared/fold_provenance.json'))
  stopifnot(identical(provenance$source,'supplied_frozen_reference'),
    identical(provenance$sha256,e$aprime_hash(root_folds)))
  inventory<-e$aprime_reference_inventory(refs)
  stopifnot(setequal(basename(names(inventory)),c('diagnosis_reference.csv','primary_subject_folds.csv')))
  writeLines('another unused obsolete map',file.path(refs,'reference/diagnosis_key.csv'))
  stopifnot(identical(inventory,e$aprime_reference_inventory(refs)))
  manifest<-jsonlite::fromJSON(file.path(root_run,'identity.json'))
  e$aprime_assert_identity(manifest,root,normalizePath(input,winslash='/'),refs)
  changed_map<-fold_map;changed_map$CV_fold<-6L-changed_map$CV_fold
  write.csv(changed_map,root_folds,row.names=FALSE)
  stopped(e$aprime_assert_identity(manifest,root,normalizePath(input,winslash='/'),refs),'References changed')
  write.csv(fold_map,root_folds,row.names=FALSE)
  jsonlite::write_json(list(note='cohort metadata'),file.path(refs,'diagnosis_reference.json'))
  stopifnot('diagnosis_reference.json'%in%basename(names(e$aprime_reference_inventory(refs))))
  unlink(file.path(refs,'diagnosis_reference.json'))
  stopifnot(file.rename(root_folds,nested_folds))
  nested_run<-file.path(tmp,'nested folds run')
  stopped(e$run_outcome(input,nested_run,'check',reference_dir=refs),'TEST_BEFORE_STAGES')
  stopifnot(identical(e$aprime_hash(file.path(nested_run,'reference/primary_subject_folds.csv')),
    e$aprime_hash(nested_folds)))
  # Two candidates are ambiguous even when their contents are identical.
  stopifnot(file.copy(nested_folds,root_folds))
  stopped(e$run_outcome(input,file.path(tmp,'ambiguous run'),'check',reference_dir=refs),
    'Ambiguous reference')
  unlink(c(root_folds,nested_folds))
  no_folds_run<-file.path(tmp,'cohort only run')
  stopped(e$run_outcome(input,no_folds_run,'check',reference_dir=refs),'TEST_BEFORE_STAGES')
  stopifnot(!dir.exists(file.path(no_folds_run,'reference')),
    identical(basename(names(e$aprime_reference_inventory(refs))),'diagnosis_reference.csv'))
  cat('PASS: shared cohort selection, raw cells, root/nested fold snapshots, supplied-fold hashes, ambiguity rejection and ignored historical files; MCMC not dispatched.\n')
})
