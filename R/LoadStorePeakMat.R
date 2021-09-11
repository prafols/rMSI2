#' LoadPeakMatrix.
#' 
#' Loads a binned peaks matrix from HDD.
#'
#' @param data_path full path to .pkmat file where data is stored.
#'
#' @return  a rMSI2 peak list object of the rMSIprocPeakMatrix class.
#' @export
#'
LoadPeakMatrix <- function( data_path )
{
  #Remove extension if supplied
  fname <- basename(data_path)
  fname <- (unlist(strsplit(fname, split = ".", fixed = T)))[1]
  fname <- paste0(fname, ".pkmat")
  data <- readRDS(file = file.path(dirname(data_path), fname))
  
  if( class(data) != "rMSIprocPeakMatrix")
  {
    rm(data)
    gc()
    stop("The provided files does not contain a valid rMSIprocPeakMatrix objext\n")
  }
  
  return(data)
}

#' StorePeakMatrix.
#'
#' Stores a binned peaks matrix to HDD.
#' Data is stored compressed using RData format with .pkmat extension.
#'
#' @param data_path full path including filename where data must be stored.
#' @param data a rMSI2 peak list object of the rMSIprocPeakMatrix class. 
#' @export
#'
StorePeakMatrix <- function( data_path, data )
{
  if( class(data) != "rMSIprocPeakMatrix")
  {
    stop("The provided peak matrix is not in rMSIprocPeakMatrix format\n")
  }
  
  #Remove extension if supplied
  fname <- basename(data_path)
  fname <- (unlist(strsplit(fname, split = ".", fixed = T)))[1]
  fname <- paste0(fname, ".pkmat")
  
  dir.create(dirname(data_path), recursive = T, showWarnings = F)
  saveRDS(data, file= file.path( dirname(data_path), fname ))
}

#' buildImgIdVectorFromPeakMatrix.
#' 
#' Builds a integer vector containing all rMSI objects ID's accroding the peak matrix row order.
#' The resulting ID vector can be used to locate a spectrum ID in a peak matrix of multiple data sets.
#'
#' @param pkMat an rMSIproc peak matrix object.
#'
#' @return an integer vector with all ID's of each image in the pkMat object.
#' @export
#'
buildImgIdVectorFromPeakMatrix <- function(pkMat)
{
  imgIDs <- rep(NA, sum(pkMat$numPixels))
  istart <- 1
  for( i in 1:length(pkMat$numPixels))
  {
    istop <- istart + pkMat$numPixels[i] - 1
    imgIDs[istart:istop] <- 1:pkMat$numPixels[i]
    istart <- istop + 1
  }
  return(as.integer(imgIDs))
}

#' getImgIdsFromPeakMatrixRows.
#' 
#' Calculate the rMSI object ID's corresponding to a subset of row of a rMSIproc peak matrix.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param rows a vector of peak matrix rows.
#'
#' @return a list with the rMSI object ID's and image that correpond the specified rows.
#' @export
#'
getImgIdsFromPeakMatrixRows <- function(pkMat, rows)
{
  #Sort rows assending and remove duplicates for fater performance
  rows <- sort(unique(rows), decreasing = F)
  
  id_lst <- list()
  allIds <- buildImgIdVectorFromPeakMatrix(pkMat)
  iRowStart <- 1
  for( i in 1:length(pkMat$names))
  {
    iRowStop <- iRowStart + pkMat$numPixels[i] - 1
    iSubRows<- which(rows %in% iRowStart:iRowStop )
    if( length(iSubRows) > 0)
    {
      id_lst[[i]] <- list( name = pkMat$names[i], uuid = pkMat$uuid[i], id = allIds[rows[iSubRows]] , pkMatRow = rows[iSubRows])
    }
    else
    {
      id_lst[[i]] <- list( name = pkMat$names[i], uuid = pkMat$uuid[i], id = c() , pkMatRow = c())
    }
    iRowStart <- iRowStop + 1
  }
  return(id_lst)
}

#' getPeakMatrixRowsFromImgIds.
#' 
#' Obtains a vector of rMSIproc peak matrix rows corresponding to a vector of rMSI obj ID's.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param img_num the number of the data set object to select in pkMat.
#' @param ids a vector of rMSI obj ID's.
#'
#' @return a vector containing the rows of peak matrix that correspond to the selected rMSI object ID's.
#' @export
#'
getPeakMatrixRowsFromImgIds <- function( pkMat, img_num, ids )
{
  if(length(pkMat$numPixels) < img_num)
  {
    stop(paste("Error: pkMat does not containg", img_num, "images.\n"))
  }
  
  if(img_num > 1)
  {
    startRow <- sum(pkMat$numPixels[1:(img_num-1)]) + 1 
  }
  else
  {
    startRow <- 1
  }
  stopRow <- startRow + pkMat$numPixels[img_num] - 1
  selRows <- ids + startRow - 1
  
  if(max(selRows) > stopRow || min(selRows) < startRow)
  {
    stop("Error: Id selection is out of range.\n")
  }
  
  return(selRows)
}
