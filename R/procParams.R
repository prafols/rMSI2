#########################################################################
#     rMSI2 - R package for MSI data handling and visualization
#     Copyright (C) 2021 Pere Rafols Soler
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

###############################################################
###            DATA DESCRIPTION R5 OBJECTS                  ###
###############################################################

DataInfo <- setRefClass("DataInfo", 
                         fields = list(
                          version = "character",
                          raw_data_path = "data.frame",
                          data_is_peaklist = "logical",
                          roi_list = "list",
                          outputpath = "character",
                          fixBrokenUUID = "logical"
                          ),
                         
                         #Constructor
                         method = list(
                           initialize = function(..., version = as.character(packageVersion("rMSI2")))
                           {
                             raw_data_path <<- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("imzML", "subimage_roi_xml")) 
                             data_is_peaklist <<- F
                             fixBrokenUUID <<- F
                             callSuper(..., raw_data_path = raw_data_path,  version = version)
                           })
                        )
DataInfo$methods(
  
  #' getNumberOfImages.
  #'
  #' get the number of subimages listed in a DataInfo object.
  #' The parseROIs() method must be run befor in order to fill the subimage list.
  #' 
  #' @return the number of images.
  #' 
  #'
  getNumberOfImages = function()
  {
    if(length(roi_list) == 0)
    {
      stop("ERROR: empty roi_list, call parseROIs() first!")
    }
    
    return (length(roi_list))
  },
  
  
  #' getImgPathPos.
  #' 
  #' Get a single image information: name, path and roi positions.
  #'
  #' @param index image index.
  #'
  #' @return a list with name, path and roi position.
  #'
  getImgPathPos = function(index)
  {
    if(index < 1 || index > length(roi_list))
    {
      stop("ERROR: index out of roi_list bounds")
    }
    
    if(length(roi_list) == 0)
    {
      stop("ERROR: empty roi_list, call parseROIs() first!")
    }
    
    return(roi_list[[index]])
  },
  
  #' appendImzMLDataPath.
  #' 
  #' Append a imzML raw data path to the list of images to process. 
  #' A given imzML file (specified in the path argument) may contain multiple MS images (in example: for an acqusition with multiple tissue sections).
  #' The subimage_roi_xml XML file must contain the pixel positions of each pixel in each image ROI. Then, the imzML will be splitted to multiple images.
  #'
  #' @param path_imzML complete path to the data imzML file.
  #' @param subimage_roi_xml optional path to the XML file containing the ROIs to split the imzML file in multiple images.
  #'
  #'
  appendImzMLDataPath = function(path_imzML, subimage_roi_xml = NA)
  {
    aux_df <- setNames( data.frame( path_imzML, subimage_roi_xml), names(raw_data_path))
    raw_data_path <<- rbind(raw_data_path, aux_df)
  },
  
  #' parseROIs.
  #' 
  #' Parses the ROI information previously loaded in the raw_data_path data.frame with the appendImzMLDataPath() method.
  #' The results are sotred in the roi_lst list. Each item in the roi_lst corresponds to a single MS image to process.
  #'
  parseROIs = function()
  {
    roi_list <<- list() #Clear previous roi list
    for(i in 1:nrow(raw_data_path))
    {
      cat(paste0("Parsing ROI info of imzML ", i ,  " of ", nrow(raw_data_path), "\n"))
      
      if(!is.na(raw_data_path$subimage_roi_xml[i]))
      {
        rois <- ParseBrukerXML(raw_data_path$subimage_roi_xml[i])  
        for( roi in rois)
        {
          roi_list[[length(roi_list)+1]] <<- list( name = paste0(unlist(strsplit(basename( as.character(raw_data_path$imzML[i])), split = "\\."))[1] , "_ROI_", roi$name ), 
                                                imzML =  as.character(raw_data_path$imzML[i]), 
                                                ROIpos = roi$pos)
        }
        
      }
      else
      {
        #No ROI file provided so just a single image
        roi_list[[length(roi_list)+1]] <<- list( name = unlist(strsplit(basename(  as.character(raw_data_path$imzML[i])), split = "\\."))[1], 
                                               imzML =  as.character(raw_data_path$imzML[i]), 
                                               ROIpos = NULL) #NULL pos means all pixels in the imzML
      }
    }
  },
  
  #' setOutputPath.
  #' set the output directory to store the results of data processing.
  #' Setting the output path is mandatory.
  #' 
  setOutputPath = function(datapath)
  {
    outputpath <<- path.expand(datapath)
  },
  
  #' setImzMLIsPeakList.
  #' set the input data imzML files as peak lists instead of spectral data.
  #' 
  setImzMLIsPeakList = function(bIsPeakList)
  {
    data_is_peaklist <<- bIsPeakList
  },
  
  
 #' setFixBrokenUUID.
 #' set to true if UUID in the ibd files must be fixed in case of a mismatch.
 #'
  setFixBrokenUUID = function(bFixBrokenUUID)
  {
    fixBrokenUUID <<- bFixBrokenUUID
  }
)

#TODO document me!
#' ImzMLDataDescription.
#' 
#' Create an empty DataInfo object to latter fill it.
#'
#' @return
#' @export
#'
#' @examples
ImzMLDataDescription <- function()
{
  data_info <- DataInfo()
  return (data_info)
}

###############################################################
###            PROCESSING PARAMETERS R5 OBJECTS             ###
###############################################################

SmoothingParams <- setRefClass("SmoothingParams", 
                               fields = list(
                                 enable = "logical",
                                 kernelSize = "integer"
                               ),
                               
                               #Constructor
                               method = list(
                                 initialize = function(...,
                                                       enable = T,
                                                       kernelSize = as.integer(7)
                                 )
                                 {
                                   callSuper(..., enable = enable, kernelSize = kernelSize)
                                 })
)

AlignmentParams <- setRefClass("AlignmentParams", 
                               fields = list(
                                 enable = "logical",
                                 bilinear = "logical",
                                 iterations = "integer",
                                 maxShiftppm = "numeric",
                                 refLow = "numeric",
                                 refMid = "numeric",
                                 refHigh = "numeric",
                                 overSampling = "integer",
                                 winSizeRelative = "numeric"
                               ),
                               
                               #Constructor
                               method = list(
                                 initialize = function(...,
                                                       enable = T,
                                                       bilinear = F,
                                                       iterations = as.integer(1),
                                                       maxShiftppm = 200,
                                                       refLow = 0.1,
                                                       refMid = 0.5,
                                                       refHigh = 0.9,
                                                       overSampling = as.integer(2),
                                                       winSizeRelative = 0.6
                                                       )
                                 {
                                   callSuper(..., enable = enable, bilinear = bilinear, iterations = iterations, maxShiftppm = maxShiftppm, 
                                             refLow = refLow, refMid = refMid, refHigh = refHigh, overSampling = overSampling, winSizeRelative = winSizeRelative)
                                 })
)


PeakPickingParams <- setRefClass("PeakPickingParams", 
                               fields = list(
                                 enable = "logical",
                                 SNR = "numeric",
                                 WinSize = "integer",
                                 overSampling = "integer"
                               ),
                               
                               #Constructor
                               method = list(
                                 initialize = function(...,
                                                       enable = T,
                                                       SNR = 5,
                                                       WinSize = as.integer(20),
                                                       overSampling = as.integer(10)
                                                       )
                                 {
                                   callSuper(..., enable = enable, SNR = SNR, WinSize = WinSize, overSampling = overSampling)
                                 })
)

PeakBinningParams <- setRefClass("PeakBinningParams", 
                                 fields = list(
                                   enable = "logical",
                                   tolerance = "numeric",
                                   tolerance_in_ppm = "logical",
                                   binFilter = "numeric"
                                 ),
                                 
                                 #Constructor
                                 method = list(
                                   initialize = function(...,
                                                         enable = T,
                                                         tolerance = 6,
                                                         tolerance_in_ppm = F,
                                                         binFilter = 0.05,
                                                         fillpeaks = T
                                   )
                                   {
                                     callSuper(..., enable = enable, tolerance = tolerance, tolerance_in_ppm = tolerance_in_ppm, binFilter = binFilter)
                                   })
)

PreProcParams <- setRefClass("PreProcParams", 
                             fields = list(
                              merge = "logical", #TRUE to process multiple images using a common mass axis
                              smoothing = "SmoothingParams",
                              alignment = "AlignmentParams",
                              massCalibration = "logical",
                              peakpicking = "PeakPickingParams",
                              peakbinning = "PeakBinningParams"
                              ),
                             
                             #Constructor
                             method = list(
                               initialize = function(...,
                                                     merge = T,
                                                     massCalibration = T
                                                     )
                               {
                                 callSuper(..., merge = merge, massCalibration = massCalibration)
                               })
                            )

PeakAnnotationParams <- setRefClass("PeakAnnotationParams", 
                                    fields = list(
                                      ppmMassTolerance = "numeric", 
                                      isotopeLikelihoodScoreThreshold = "numeric",
                                      adductElementsTable = "data.frame"
                                    ),
                                    
                                    #Constructor
                                    method = list(
                                      initialize = function(...,
                                                            ppmMassTolerance = as.integer(1),
                                                            isotopeLikelihoodScoreThreshold = 0.8,
                                                            adductElementsTable = data.frame(name = c("+K","+Na","+H"),
                                                                                             mass = c(38.963706,22.98976,1.007825),
                                                                                             priority = c(0,0,0))
                                      )
                                      {
                                        callSuper(..., ppmMassTolerance = ppmMassTolerance, 
                                                  isotopeLikelihoodScoreThreshold = isotopeLikelihoodScoreThreshold,
                                                  adductElementsTable = adductElementsTable)
                                      })
)

ProcParams <- setRefClass("ProcParams",
                          fields = list(
                           version = "character",
                           preprocessing = "PreProcParams",
                           peakAnnotation = "PeakAnnotationParams"
                           ),
                          
                          #Constructor
                          method = list(
                            initialize = function(..., version = as.character(packageVersion("rMSI2")))
                                          {
                                            callSuper(..., version = version)
                                          })
                         )

#Set read-only fields
ref <- getRefClass("ProcParams")
ref$lock("version")

#Adding methods to the ProcParams class
ProcParams$methods(
                    getMergedProcessing = function()
                    {
                      return(preprocessing$merge)
                    },
                    
                    setMergedProcessing = function(bMerge)
                    {
                      preprocessing$merge <<- bMerge
                    }
                    
                   )

#TODO document me!
#' ProcessingParameters.
#' 
#' Create an empty ProcParams object to latter fill it.
#' 
#' @return
#' @export
#'
#' @examples
ProcessingParameters <- function()
{
  proc_params <- ProcParams()
  return (proc_params)
}

#' LoadProcParams.
#' 
#' Loads processing parameters from the HDD.
#'
#' @param filename full path to the .params file where processing paramers are stored.
#'
#' @return  a rMSI2 processing parameters object.
#' @export
#'
LoadProcParams <- function( filename )
{
  #Remove extension if supplied
  fname <- basename(filename)
  fname <- sub("\\.[^.]*$", "", fname)
  fname <- paste0(fname, ".params")
  data <- readRDS(file = file.path(dirname(filename), fname))
  
  if( class(data) != "ProcParams")
  {
    rm(data)
    gc()
    stop("The provided files does not contain a valid ProcParams object\n")
  }
  
  return(data)
}

#' StoreProcParams.
#'
#' Stores all processing parameters to HDD.
#' Parameters arestored compressed using RData format with .params extension.
#'
#' @param filename full path to the .params file where processing paramers are stored.
#' @param params processing parameters object of the ProcParams class. 
#' @export
#'
StoreProcParams <- function( filename, params )
{
  if( class(params) != "ProcParams")
  {
    stop("The provided parameters are not in ProcParams format\n")
  }
  
  #Remove extension if supplied
  fname <- basename(filename)
  fname <- sub("\\.[^.]*$", "", fname)
  fname <- paste0(fname, ".params")
  
  dir.create(dirname(filename), recursive = T, showWarnings = F)
  saveRDS(params, file= file.path( dirname(filename), fname ))
}
