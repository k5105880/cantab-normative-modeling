here<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
local({
  before_env<-Sys.getenv();before_opts<-options();before_search<-search();before_global<-ls(.GlobalEnv,all.names=TRUE)
  scratch<-tempfile('source_isolation_');dir.create(scratch);old<-setwd(scratch)
  on.exit({setwd(old);unlink(scratch,recursive=TRUE)},add=TRUE)
  e<-new.env(parent=.GlobalEnv)
  for(file in list.files(here,pattern='[.]R$',full.names=TRUE))
    if(!grepl('/test_',file))source(file,local=e)
  stopifnot(identical(before_env,Sys.getenv()),identical(before_opts,options()),identical(before_search,search()),
    identical(before_global,ls(.GlobalEnv,all.names=TRUE)),!length(list.files(scratch,all.files=TRUE,no..=TRUE)))
})
cat('PASS: source isolation for all implementation and entry-point files.\n')
