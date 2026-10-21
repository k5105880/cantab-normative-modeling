# CANTAB retry2: audit the actual brms population-level likelihood design.
# Include penalized smooth bases Zs_*, not participant random-effect Z_*.
# Separate response/distributional predictors must be audited separately.
cantab_audit_standata <- function(standata, label, output_file=NULL, tolerance=1e-9) {
  if (!is.list(standata)) stop('standata must be a named list')
  selected <- names(standata)[grepl('^(X($|_)|Xs($|_)|Zs_)', names(standata)) &
                              vapply(standata, is.matrix, logical(1))]
  if (!length(selected)) stop('No population design matrices found: ',label)
  parameter <- function(nm) {
    if (nm %in% c('X','Xs')) return('mu')
    if (startsWith(nm,'Zs_')) {
      stem <- sub('_[0-9]+_[0-9]+$', '', nm)
      if (identical(stem,'Zs')) return('mu')
      return(sub('^Zs_', '', stem))
    }
    sub('^Xs?_', '', nm)
  }
  groups <- split(selected,vapply(selected,parameter,character(1)))
  rows <- lapply(names(groups),function(p) {
    keys <- groups[[p]]
    blocks <- standata[keys]
    nr <- vapply(blocks,nrow,integer(1))
    if (length(unique(nr))!=1L || !length(nr) || nr[1]<1L)
      stop('Inconsistent or empty design row counts: ',label,' / ',p)
    x <- do.call(cbind,blocks)
    if (!is.numeric(x) || any(!is.finite(x)))
      stop('Nonfinite design: ',label,' / ',p)
    if (!ncol(x)) stop('Empty design: ',label,' / ',p)
    norms <- sqrt(colSums(x*x))
    if (any(!is.finite(norms))) stop('Nonfinite column norm: ',label,' / ',p)
    scaled <- sweep(x,2,ifelse(norms>0,norms,1),'/')
    rank <- qr(scaled,tol=tolerance,LAPACK=FALSE)$rank
    singular <- svd(scaled,nu=0,nv=0)$d
    ratio <- if (rank<ncol(x)) 0 else min(singular)/max(singular)
    data.frame(label=label,parameter=p,n=nrow(x),columns=ncol(x),rank=rank,
               full_rank=rank==ncol(x),relative_min_singular=ratio,
               matrices=paste(keys,collapse=';'),stringsAsFactors=FALSE)
  })
  result <- do.call(rbind,rows)
  if (!is.null(output_file)) {
    dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
    present <- file.exists(output_file) && file.info(output_file)$size>0
    write.table(result,output_file,sep=',',row.names=FALSE,
                col.names=!present,append=present,quote=TRUE)
  }
  if (any(!result$full_rank)) {
    bad <- result[!result$full_rank,]
    stop('Rank-deficient population design: ',label,' / ',
         paste(sprintf('%s %d/%d',bad$parameter,bad$rank,bad$columns),collapse=', '))
  }
  result
}
