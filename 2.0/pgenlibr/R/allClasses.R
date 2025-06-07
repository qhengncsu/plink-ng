setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("charOrNULL", c("character", "NULL"))

setClass("PlinkMatrix", representation(ptr = "externalptr", samples = "integer",
                                       variants = "integer", fname = "character", xim = "numeric", Dim = "integer", ncov = "integer", covs="matrixOrNULL", colname="charOrNULL"))

#' @export
PlinkMatrixbase <- function(fname, samples, variants, covs=NULL) {
  if (is.unsorted(samples, strictly = TRUE)) {
    stop("The subset indices must be sorted in strictly increasing order.")
  }

  if (length(variants) != length(unique(variants))) {
    stop("Must not have duplicated variant index.")
  }
  dimensions <- c(length(samples), length(variants))
  ncov = 0L
  if(!is.null(covs)) {
    if(nrow(covs) != dimensions[1]) {
      stop("Covariates must have same number of rows as SNP Matrix")
    }
    ncov = ncol(covs)
    dimensions[2] = dimensions[2] + ncov
  }

  samples <- as.integer(samples)
  variants <- as.integer(variants)


  new("PlinkMatrix", ptr = new("externalptr"), samples = samples, variants = variants, fname = fname,
      Dim = dimensions, xim = 0, ncov = ncov, covs = covs, colname=NULL
  )
}
setGeneric("actualize", function(x) standardGeneric("actualize"))
setGeneric("setcovs", function(x, newcov) standardGeneric("setcovs"))


#' @export
setMethod(
  "actualize",
  "PlinkMatrix",
  function(x) {
    if(x@fname == "fromPgen"){
      warning("already actualized")
      return(x)
    }
    x@xim = double(length(x@variants))
    x@ptr = getcompactptr(x@fname, x@variants,x@samples, x@xim)
    return(x)
  }
)

#' @export
PlinkMatrixFromPgen <- function(pgen, variants, covs=NULL, variant_names=NULL){
  if (length(variants) != length(unique(variants))) {
    stop("Must not have duplicated variant index.")
  }
  sample_size = integer(1)
  xim = double(length(variants))
  xptr = getcompactptrfromPgen(pgen, variants, xim, sample_size)
  dimensions <- c(sample_size, length(variants))
  ncov <- integer(1)
  if(!is.null(covs)) {
    if(nrow(covs) != dimensions[1]) {
      stop("Covariates must have same number of rows as SNP Matrix")
    }
    ncov = ncol(covs)
    dimensions[2] = dimensions[2] + ncov
  }
  colname = NULL
  if(!is.null(variant_names)){
    cov_names = colnames(covs)
    if((!is.null(covs)) && is.null(cov_names)){
      cov_names = paste0("V", 1:ncol(covs))
    }
    colname = c(cov_names, variant_names)

  }
  new("PlinkMatrix", ptr = xptr, samples = integer(1), variants = variants, fname = "fromPgen",
      Dim = dimensions, xim = xim, ncov = ncov, covs = covs, colname=colname
  )
}

#' @export
setMethod(
  "setcovs",
  c("PlinkMatrix", "matrix"),
  function(x, newcov) {
    if(!is.null(x@covs)){
      stop("Already has covariates")
    }
    if(nrow(newcov) != x@Dim[1]){
      stop("covariates dimension does not match that of the snp matrix")
    }
    x@covs = newcov
    x@Dim[2] = x@Dim[2] + ncol(newcov)
    x@ncov = ncol(newcov)
    cov_names = colnames(newcov)
    if(is.null(cov_names)){
      cov_names = paste0("V", 1:ncol(newcov))
    }
    x@colname = c(cov_names, x@colname)
    return(x)
  }
)

#' @export
setMethod(
  "dim",
  "PlinkMatrix",
  function(x) {
    return(x@Dim)
  }
)

#' @export
prepareFeatures <- function(pgen, vars, names){
  plinkmatrix <- pgenlibr::PlinkMatrixFromPgen(pgen, match(names, vars), variant_names=names)
  plinkmatrix
}
