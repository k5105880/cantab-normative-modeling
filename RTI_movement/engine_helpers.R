# Run identity, checkpoint integrity and dependency checks.
engine_required_packages <- c('dplyr','brms','posterior','loo','cmdstanr','digest','mgcv','future')
engine_check_packages <- function() {
  missing <- engine_required_packages[!vapply(engine_required_packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing))stop('Missing R packages: ',paste(missing,collapse=', '))
  invisible(TRUE)
}
engine_atomic_save <- function(value,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temporary<-tempfile('checkpoint_',tmpdir=dirname(path));on.exit(unlink(temporary),add=TRUE)
  saveRDS(value,temporary)
  if(!file.rename(temporary,path))stop('Cannot save checkpoint: ',path)
}
engine_file_manifest <- function(files) {
  if(any(!file.exists(files)))stop('Missing run artifact')
  data.frame(path=normalizePath(files,mustWork=TRUE),sha256=unname(vapply(files,function(p)digest::digest(file=p,algo='sha256'),character(1))))
}
engine_guard <- function(root,input,run_dir) {
  code<-sort(list.files(root,pattern='[.]R$',full.names=TRUE))
  identity<-list(files=engine_file_manifest(c(code,input)),R=R.version.string,
    packages=vapply(engine_required_packages,function(p)as.character(utils::packageVersion(p)),character(1)),
    cmdstan=as.character(cmdstanr::cmdstan_version()))
  path<-file.path(run_dir,'engine_identity.rds')
  if(file.exists(path)) {
    if(!identical(readRDS(path),identity))stop('Input, code or software changed; use a fresh output directory')
  } else {
    stale<-list.files(run_dir,pattern='^(normative_objects|normative_complete|brms_cache|model_checkpoints|cv_checkpoints)')
    if(length(stale))stop('Existing model artifacts without a run identity')
    engine_atomic_save(identity,path)
    capture.output(sessionInfo(),file=file.path(run_dir,'session_info.txt'))
  }
  invisible(identity)
}
engine_complete <- function(run_dir) {
  path<-file.path(run_dir,'normative_complete.rds');if(!file.exists(path))return(FALSE)
  seal<-readRDS(path)
  if(!identical(seal,engine_file_manifest(seal$path)))stop('Completed output files changed; resume refused')
  TRUE
}
engine_mark_complete <- function(run_dir) {
  files<-c(list.files(run_dir,pattern='[.]csv$',full.names=TRUE),file.path(run_dir,'normative_objects.rds'),
    list.files(file.path(run_dir,'audit'),full.names=TRUE))
  engine_atomic_save(engine_file_manifest(files),file.path(run_dir,'normative_complete.rds'))
}
