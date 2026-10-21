args<-gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
root<-dirname(dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))))
local({
  scratch<-tempfile('aprime_source_');dir.create(scratch)
  previous<-setwd(scratch);on.exit({setwd(previous);unlink(scratch,recursive=TRUE)},add=TRUE)
  env_before<-Sys.getenv();options_before<-options();search_before<-search()
  globals_before<-ls(.GlobalEnv,all.names=TRUE)
  source_env<-new.env(parent=.GlobalEnv)
  files<-c(list.files(file.path(root,'R'),pattern='[.]R$',full.names=TRUE),file.path(root,'run.R'))
  for(file in files)source(file,local=source_env)
  stopifnot(identical(env_before,Sys.getenv()),identical(options_before,options()),
    identical(search_before,search()),identical(globals_before,ls(.GlobalEnv,all.names=TRUE)),
    !length(list.files(scratch,all.files=TRUE,no..=TRUE)),
    is.function(source_env$run_outcome),is.function(source_env$aprime_driver),
    is.function(source_env$cantab_audit_standata),
    identical(names(formals(source_env$build_model)),'data'))
})
cat('PASS: every implementation file and entry point sources without changing environment, options, search path, global objects or working-directory files.\n')
