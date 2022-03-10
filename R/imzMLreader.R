#########################################################################
#     rMSI - R package for MSI data handling and visualization
#     Copyright (C) 2014 Pere Rafols Soler
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

#' fixImzMLDuplicates.
#' Delete duplicates and possible zero-drop errors (fixing Bruker's bugs in imzML).
#' 
#' @param mass spectrum mass axis.
#' @param intensity spectrum intensity.
#'
#' @return a list with non-duplicated mass and intensity vectors.
#'
fixImzMLDuplicates <- function(mass, intensity)
{
  idup <- which(duplicated(mass))
  if(length(idup) > 0)
  {
    for( i in 1:length(idup))
    {
      intensity[idup[i] - 1] <- max( intensity[idup[i]], intensity[idup[i] - 1])
    }
    mass <- mass[-idup]
    intensity <- intensity[-idup]
  }
  return(list(mass=mass, intensity=intensity))
}

#' import_imzML.
#'
#' @param imzML_File full path to .imzML file 
#' @param ibd_File path to the binary file (default the same as imzML file but with .ibd extension)
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param fun_text This is a callback function to update the label widget of loading data. See details for more information.
#' @param close_signal function to be called if loading process is abored.
#' @param verifyChecksum if the binary file checksum must be verified, it is disabled by default for convenice with really big files.
#' @param convertProcessed2Continuous if true (the default) an imzML file in processed mode will be converted to a continuous mode.
#' @param subImg_rename alternative image name, new rMSI files will be created with the given name.
#' @param subImg_Coords a Complex vector with the motors coordinates to be included in the rMSI data.
#' @param fixBrokenUUID set to FALSE by default to automatically fix an uuid mismatch between the ibd and the imzML files (a warning message will be raised).
#'
#'  Imports an imzML image to an rMSI data object.
#'  It is recomanded to use rMSI2::LoadMsiData directly instead of this function.
#'
#' @return an rMSI data object.
#' @export
#' 
import_imzML <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ),
                         fun_progress = NULL, fun_text = NULL, close_signal = NULL, 
                         verifyChecksum = F, convertProcessed2Continuous = T,
                         subImg_rename = NULL, subImg_Coords = NULL, fixBrokenUUID = F)
{
  setPbarValue<-function(progress)
  {
    setTxtProgressBar(pb, progress)
    return(T)
  }

  #1- Parse XML data
  if(is.null(fun_text))
  {
    cat("Parsing XML data in imzML file...\n")
  }
  else
  {
    fun_text("Parsing XML data in imzML file...")
  }
  xmlRes <- CimzMLParse(path.expand(imzML_File))
  if( !is.null(xmlRes$Error))
  {
    .controlled_loadAbort(paste0(xmlRes$Error, "\n"), close_signal)
  }
  
  #2- Check UUID and fix it if broken
  if( !check_imzMLUUID(imzML_file = imzML_File, fixibd_uuid = fixBrokenUUID, xmlParserResult = xmlRes) )
  {
    .controlled_loadAbort("ERROR: UUID in imzML file does not match UUID in ibd file.\nYou may want to fix this with the argument fixBrokenUUID set to TRUE.\n", close_signal)
  }
  
  if(verifyChecksum)
  {
    if( xmlRes$SHA != "" )
    {
      cat("\nChecking binary data checksum using SHA-1 key... ")
      res <- toupper(digest::digest( ibd_File, algo = "sha1", file = T))
      if( res == xmlRes$SHA )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$SHA, "\nBinary file key:", res,"\n"))
        #.controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
        #Disableing the abort, just showing a warning here as it seams that there is bug in brukers checksum imzml file...
        cat("WARNING: MS data my be corrupt!\n")
      }
    }
    if( xmlRes$MD5 != "")
    {
      cat("Checking binary data checksum using MD5 key... ")
      res <- toupper(digest::digest( ibd_File, algo = "md5", file = T))
      if( res == xmlRes$MD5 )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$MD5, "\nBinary file key:", res,"\n"))
        #.controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
        #Disableing the abort, just showing a warning here as it seams that there is bug in brukers checksum imzml file...
        cat("WARNING: MS data my be corrupt!\n")
      }
    }  
  }
  else
  {
    cat("WARNING: Checksum validation is disabled, data may be corrupt!\n")
  }
  
  #Select specific pixels in the dataset
  if( !is.null(subImg_Coords))
  {
    keepIds <-  which( complex( real = xmlRes$run_data$x, imaginary = xmlRes$run_data$y) %in% subImg_Coords)

    if(length(keepIds) == 0)
    {
      .controlled_loadAbort("ERROR: no subImg_Coords found in current imzML data.\n", close_signal)
    }
    
    xmlRes$run_data <- xmlRes$run_data[keepIds,]
    
    #Provide a new name for the sub image if it has not been provided
    if(is.null(subImg_rename))
    {
      subImg_rename <- paste0(basename(imzML_File), "_subset")
    }
  }
  
  #3- Create a connection to read binary file
  bincon <- file(description = ibd_File, open = "rb")

  #4- Obtain de m/z axis (for continuous mode just read the mass axis, for processed mode calculate a common mass axis)
  pt<-proc.time()
  dataPointEncoding_Mz  <- dataPointBinaryEncoding(xmlRes$mz_dataType)
  dataPointEncoding_Int <- dataPointBinaryEncoding(xmlRes$int_dataType)
  bNoNeed2Resample <- T #Start assuming there is no need to resample the data
  readBin(bincon, integer(), 16, size = 1, signed = F) #Read the first 16 bytes (uuid) to position seek pointer
  if(xmlRes$continuous_mode)
  {
    mzAxis <- readBin(bincon, dataPointEncoding_Mz$dataType, xmlRes$run_data[1, "mzLength"], size = dataPointEncoding_Mz$bytes, signed = T)
    close(bincon)
  }
  else
  {
    close(bincon) #The file open in R must be close before going to c++
    if(convertProcessed2Continuous)
    {
      #Processed mode, so a common mass axis must be calculated and stored in mzAxis var
      rMSIDummyObj <- list(data = list( peaklist = xmlRes ))
      rMSIDummyObj$data$path <- dirname(ibd_File)
      rMSIDummyObj$data$peaklist$path <- dirname(path.expand( ibd_File))
      rMSIDummyObj$data$peaklist$file <- sub("\\.[^.]*$", "", basename(ibd_File))
      img_Dummylst <- list(rMSIDummyObj)
      newCommonMassSingleImzML <- rMSI2:::CcommonMassAxis(img_Dummylst, parallel::detectCores(), 100) 
      mzAxis <-  newCommonMassSingleImzML$mass
      bNoNeed2Resample <- newCommonMassSingleImzML$NoNeed2Resample
      rm(newCommonMassSingleImzML)
      rm(rMSIDummyObj)
      rm(img_Dummylst)
      gc()
      pt<-proc.time() - pt
      cat(paste("\nMass axis calculation time:",round(pt["elapsed"], digits = 1),"seconds\n"))
      cat(paste("The re-sampled mass axis contains", length(mzAxis), "data points\n"))
    }
    else
    {
      mzAxis <- NULL #Setting the common mass axis to NULL since data will be a peak list so it is not possible to directely create an rMSIXBin object
    }
  }
 
  #5- Create the rMSIXBin
  if( is.null(subImg_rename))
  {
    subImg_rename <-  basename(imzML_File)
  }
  img <- CreateEmptyImage(num_of_pixels = nrow(xmlRes$run_data), mass_axis = mzAxis, pixel_resolution = xmlRes$pixel_size_um,
                               img_name = subImg_rename,
                               rMSIXBin_path = path.expand(dirname(imzML_File)),
                               uuid = xmlRes$UUID
                                )

  #Fill missing data
  img$data$rMSIXBin$file <- sub("\\.[^.]*$", "", basename(imzML_File))
  #img$data$rMSIXBin$uuid has been set by CreateEmptyImage
  
  img$data$imzML$file <- sub("\\.[^.]*$", "", basename(imzML_File))
  #img$data$imzML$uuid has been set by CreateEmptyImage
  if( xmlRes$SHA != "" )
  {
    img$data$imzML$SHA <- xmlRes$SHA
  }
  if( xmlRes$MD5 != "" )
  {
    img$data$imzML$MD5 <- xmlRes$MD5
  }
  
  
  if (xmlRes$continuous_mode )
  {
    img$data$imzML$continuous_mode <- T
  }
  else
  {
    if(convertProcessed2Continuous)
    {
      img$data$imzML$continuous_mode <- convertProcessed2Continuous & bNoNeed2Resample
    }
    else
    {
      img$data$imzML$continuous_mode <- F
    }
  }
  
  img$data$imzML$mz_dataType <- xmlRes$mz_dataType
  img$data$imzML$int_dataType <-xmlRes$int_dataType
  img$data$imzML$run <- xmlRes$run_data
  
  #Compute image size and arrange coords to avoid holes
  img$posMotors[, "x"] <- xmlRes$run_data[, "x"]
  img$posMotors[, "y"] <- xmlRes$run_data[, "y"]
  img$pos <- remap2ImageCoords(img$posMotors)
  img$size["x"] <- max(img$pos[,"x"])
  img$size["y"] <- max(img$pos[,"y"])
  
  #6- It is a peak list, so, set references to it
  if(!img$data$imzML$continuous_mode && !convertProcessed2Continuous)
  {
    
    img$data$peaklist <- list( UUID = img$data$imzML$uuid,
                               continuous_mode = F,
                               compression_mz = F,
                               compression_int = F,
                               MD5 = img$data$imzML$MD5,
                               SHA = img$data$imzML$SHA,
                               mz_dataType = img$data$imzML$mz_dataType,
                               int_dataType = img$data$imzML$int_dataType,
                               pixel_size_um = img$pixel_size_um,
                               run_data = img$data$imzML$run,
                               path = img$data$path,
                               file = img$data$imzML$file,
                               rMSIpeakList = xmlRes$rMSIpeakList
                              )
    
    img$data$imzML <- NULL #It is a peak list, so no spectra in profile mode
    img$data$rMSIXBin <- NULL #It is a peak lis, so no pngstream can be created
  }

  #7- And it's done, just return de rMSI object
  gc()
  return(img)
}

#' dataPointBinaryEncoding
#' 
#' Get the number of bytes and data type that must be readed in binary streams for each encoded data point.
#'
#' @param str_dataType a string with the C-style data type: int, long, float and double.
#'
#' @return a list containing an integer with the number of bytes to read for each data point and the R data type.
#'
dataPointBinaryEncoding <- function(str_dataType)
{
  result <- list( bytes = NULL,  dataType = NULL)
  if(str_dataType == "int" || str_dataType == "float")
  {
    result$bytes <- 4
  }
  if(str_dataType == "long" || str_dataType == "double")
  {
    result$bytes <- 8
  }
  
  if(str_dataType == "int" || str_dataType == "long")
  {
    result$dataType <- integer()
  }
  if(str_dataType == "float" || str_dataType == "double")
  {
    result$dataType <- numeric()
  }

  return(result)
}

#' check_imzMLUUID.
#' 
#' Check and fix the UUID of invalid imzML files.
#'
#' @param imzML_file path to an imzML file.
#' @param fixibd_uuid set to TRUE to overwrite the ibd file uuid with the imzML uuid.
#' @param xmlParserResult result of the CimzMLParse() if available
#'
#' @return TRUE if the UUID is OK.
#' @export
#'
check_imzMLUUID <- function(imzML_file, fixibd_uuid = F, xmlParserResult = NULL)
{
  #Set files
  imzML_fname <- path.expand(imzML_file)
  ibd_fname <- basename(imzML_fname)
  ibd_fname <- sub("\\.[^.]*$", "", ibd_fname)
  ibd_fname <- file.path(dirname(imzML_fname), paste0(ibd_fname, ".ibd"))
  
  if(!file.exists(imzML_fname))
  {
    stop("ERROR: imzML file not found")
  }
  
  if(!file.exists(ibd_fname))
  {
    stop("ERROR: ibd file not found")
  }
  
  if(is.null(xmlParserResult))
  {
    cat("Parsing imzML file...\n")
    imzML_XML <- CimzMLParse(imzML_fname) 
    imzMLUUID <- imzML_XML$UUID
    rm(imzML_XML)
    gc()
  }
  else
  {
    imzMLUUID <- xmlParserResult$UUID
  }
  
  #Check imzML UUID
  cat("Checking ibd UUID...")
  bincon <- file(description = ibd_fname, open = "rb")
  binUUID <- paste(sprintf("%.2X", readBin(bincon, integer(), 16, size = 1, signed = F)), collapse = "")
  close(bincon)
  if(binUUID == imzMLUUID)
  {
    cat("OK\n")
    return(TRUE)
  }
  else
  {
    if(fixibd_uuid)
    {
      cat("NOK\n")
      overwriteIbdUUid(ibd_fname, imzMLUUID)
      cat("WARNING: UUID mismatch, overwriting it in the ibd file!\n")
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }
  }
}
