#' rMSI2
#' 
#'  MSI data handling and visualization in R
#' 
#' @docType package
#' @author Pere Rafols
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib rMSI2
#' @name rMSI2
NULL  

.onUnload <- function (libpath)
{
  library.dynam.unload("rMSI2", libpath)
}