# Internal process entry point. Use run_outcome() from run.R.
main <- function(args=commandArgs(TRUE)) {
  mode <- match.arg(if(length(args))args[1L]else'--check',c('--check','--fresh','--resume'))
  input <- Sys.getenv('NM_INPUT_FILE','');run <- Sys.getenv('NM_RUN_DIR','')
  if(!nzchar(input)||!nzchar(run))stop('Use run_outcome(input_file,output_dir)')
  source(file.path(here,'core.R'),local=TRUE)
  if(mode=='--check') {
    source(file.path(here,'preflight.R'),local=TRUE)
    return(preflight_models(input,file.path(run,'preflight')))
  }
  run_normative_engine(input,run)
}
here <- local({
  paths<-unlist(lapply(sys.frames(),function(f)f$ofile),use.names=FALSE)
  if(!length(paths))paths<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
  dirname(normalizePath(gsub('~+~',' ',tail(paths,1L),fixed=TRUE),mustWork=TRUE))
})
if(sys.nframe()==0L)main()
