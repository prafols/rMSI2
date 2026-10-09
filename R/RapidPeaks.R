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



#' Fast Peak Matrix Binning Across Multiple Mass Spectrometry Images
#'
#' @description
#' Generates a unified peak intensity matrix from a list of MSI datasets (`rMSI` image objects) 
#' by constructing a common mass axis, detecting reference peaks on a composite overall spectrum, 
#' and extracting peak intensities across all pixel locations. Supports multi-core execution 
#' on Unix platforms via process forking.
#'
#' @param img_lst A \code{list} of \code{rMSI} image objects to merge and bin.
#' @param peak_width_scans \code{numeric} integer scalar. Peak width multiplier in scan units 
#'   used to calculate peak integration boundaries (\code{binSize * peak_width_scans}). Default is \code{6}.
#' @param max_mem_MB \code{numeric} scalar. Maximum memory limit in Megabytes (MB) allowed 
#'   for peak matrix allocation. Defaults to 20\% of available system RAM.
#' @param min_SNR \code{numeric} scalar. Minimum Signal-to-Noise Ratio threshold for 
#'   peak detection on the product spectrum. Default is \code{0.5}.
#' @param n_cores \code{integer} scalar. Number of CPU cores for parallel extraction on 
#'   Unix operating systems. Defaults to available cores minus 1.
#'
#' @details
#' The binning workflow executes in five stages:
#' \enumerate{
#'   \item \strong{Mass Axis Harmonization:} Verifies mass alignment across images or constructs 
#'         a merged common axis using \code{\link{MergeMassAxisAutoBinSize}}.
#'   \item \strong{Spectral Aggregation:} Calculates overall mean and base-peak spectra across 
#'         all images using C++ linear interpolation routines (\code{CaccumulateOverallSpectra}).
#'   \item \strong{Composite Peak Picking:} Performs peak detection on the geometric mean spectrum 
#'         \eqn{\sqrt{\text{mean} \times \text{base}}} via \code{\link{DetectPeaks}}.
#'   \item \strong{Memory Cap Evaluation:} Truncates feature selection to the top \eqn{N} peaks by SNR 
#'         to strictly enforce \code{max_mem_MB}.
#'   \item \strong{Intensity Extraction:} Populates the output matrix by extracting binned peak 
#'         intensities for each pixel coordinate across all images.
#' }
#'
#' @return An \code{rMSI} peak matrix object containing:
#' \item{mass}{\code{numeric} vector of binned peak center \emph{m/z} values.}
#' \item{binSize}{\code{numeric} vector of peak integration widths.}
#' \item{intensity}{\code{matrix} of dimension \eqn{N_{\text{pixels}} \times N_{\text{peaks}}} holding peak intensities.}
#' \item{SNR}{\code{numeric} vector of peak Signal-to-Noise ratios.}
#' \item{area}{\code{numeric} vector of integrated peak areas.}
#' \item{normalizations}{\code{data.frame} or \code{matrix} of aggregated pixel-wise normalization factors across merged images.}
#' \item{pos}{\code{matrix} of combined spatial pixel coordinates (\code{x}, \code{y}).}
#' \item{posMotors}{\code{matrix} of combined stage motor coordinates (\code{x}, \code{y}).}
#' \item{numPixels}{\code{integer} vector containing pixel counts for each merged dataset.}
#' \item{names}{\code{character} vector of merged image dataset names.}
#' \item{uuids}{\code{character} vector of unique identifiers (UUIDs) for merged datasets.}
#'
#' @seealso \code{\link{DetectPeaks}}, \code{\link{MergeMassAxisAutoBinSize}}
#'
#' @keywords internal
#' @export
FastPeakBinning<- function(img_lst, 
                peak_width_scans = 6,
                max_mem_MB = 0.2*ps::ps_system_memory()$avail/(1024^2),
                min_SNR = 0.5,
                n_cores = max(1, parallel::detectCores() - 1, na.rm = T)
                )
{

  # 1. Get Overall average and base by maxing all of them
  common_mass <- img_lst[[1]]$mass
  identicalMassAxis <- TRUE
  if( length(img_lst) > 1)
  {
    for( i in 2:length(img_lst))
    {
      identicalMassAxis <- identicalMassAxis & identical(common_mass, img_lst[[i]]$mass)
    }
  }
  
  # Calculate the new common mass axis 
  if(!identicalMassAxis)
  {
    for( i in 2:length(img_lst))
    {
      massMergeRes <- MergeMassAxisAutoBinSize(common_mass, img_lst[[i]]$mass)
      if(massMergeRes$error)
      {
        stop("ERROR: The mass axis of the images to merge is not compatible because they do not share a common range.\n")
      }
      common_mass <- massMergeRes$mass
    }
  }
  
  # Initialize output accumulation vectors and get common spectra
  n_mass <- length(common_mass)
  overall_mean_spectrum <- numeric(n_mass)
  overall_base_spectrum <- numeric(n_mass)
  
  for (i in seq_along(img_lst))
  {
    img_i <- img_lst[[i]]
    
    if (!identicalMassAxis)
    {
      # C++ handles linear interpolation (mlinterp) and pmax update in-place
      CaccumulateOverallSpectra(
        common_mass,
        img_i$mass,
        img_i$mean,
        img_i$base,
        overall_mean_spectrum,
        overall_base_spectrum )
    } 
    else
    {
      # Direct fast path when mass axes already match perfectly
      overall_mean_spectrum <- pmax(overall_mean_spectrum, img_i$mean, na.rm = TRUE)
      overall_base_spectrum <- pmax(overall_base_spectrum, img_i$base, na.rm = TRUE)
    }
  }
  
  # 2. Peak-pick prod spectrum
  prod <- sqrt(overall_base_spectrum * overall_mean_spectrum)
  ppeaks <- DetectPeaks(common_mass, prod, SNR = min_SNR, OverSampling = 20)  
  
  # Calculate max number of peaks to fit in max_mem_MB
  total_num_of_pixels = sum(unlist(lapply(img_lst, function(x){ 
    return(nrow(x$pos))
    })))
  
  max_cols_in_mem <- floor((max_mem_MB*(1024^2))/(total_num_of_pixels*8))
  max_peaks <- min(length(ppeaks$mass), max_cols_in_mem)
  
  # 3. A dataframe to control de new peak matrix generation
  feat_list <- data.frame( mass = ppeaks$mass, 
                           priority = ppeaks$SNR,  #TODO it will be nice to allow different rules to set peak-priority
                           peakwidth = ppeaks$binSize*peak_width_scans,
                           SNR = ppeaks$SNR,
                           Area = ppeaks$area )
  sub_feat_list <- feat_list[order(feat_list$priority, decreasing = T)[1:max_peaks], ]
  sub_feat_list <- sub_feat_list[order(sub_feat_list$mass, decreasing = F), ] #Sort by mass
  rm(feat_list)
  gc()
  
  # Cache local variables and mem allocation
  masses  <- sub_feat_list$mass
  peakWidths  <- sub_feat_list$peakwidth
  n_peaks <- length(sub_feat_list$mass)
  peakMatrixInt <- matrix(0.0, nrow = total_num_of_pixels, ncol = n_peaks) 
  
  # 4. Execute parallel loop
  pkmat_current_first_row <- 1
  for (i in seq_along(img_lst))
  {
    cat(paste0("Working on image ", i, " of ", length(img_lst), " ...\n"))
    img_i <- img_lst[[i]]
    pos_idx <- img_i$pos
    n_pixels_i <- nrow(pos_idx)
    pkmat_current_last_row <- pkmat_current_first_row + n_pixels_i - 1

    # Call C++ backend
    C_RapidBinning(img_i, masses, peakWidths, peakMatrixInt, pkmat_current_first_row, n_cores)
      
    pkmat_current_first_row <- pkmat_current_last_row + 1
  }
  
  # 5. Combine matrix
  PkMat <- list( mass = sub_feat_list$mass, 
                 binSize = sub_feat_list$peakwidth, 
                 intensity = peakMatrixInt, 
                 SNR = sub_feat_list$SNR,
                 area = sub_feat_list$Area)
  
  gc()

  
  #Append normalizations to the peak matrix
  PkMat$normalizations <- img_lst[[1]]$normalizations
  if(length(img_lst) > 1)
  {
    for( i in 2:length(img_lst))
    {
      PkMat$normalizations <- rbind(PkMat$normalizations, img_lst[[i]]$normalizations)
    }
  }
  
  
  #Add a copy of img$pos to the peakMatrix
  mergedNames <- unlist(lapply(img_lst, function(x){ return(x$name) }))
  mergedNumPixels <- unlist(lapply(img_lst, function(x){ return(nrow(x$pos)) }))
  mergedPos <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
  mergedMotors <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
  mergedUUIDs <- unlist(lapply(img_lst, function(x){ return(x$data$imzML$uuid) }))
  colnames(mergedPos) <- c("x", "y")
  colnames(mergedMotors) <- c("x", "y")
  istart <- 1
  for( i in 1:length(img_lst))
  {
    istop <- istart + nrow(img_lst[[i]]$pos) - 1
    mergedPos[ istart:istop , "x"] <- img_lst[[i]]$pos[, "x"]
    mergedPos[ istart:istop , "y"] <- img_lst[[i]]$pos[, "y"]
    
    mergedMotors[ istart:istop , "x"] <- img_lst[[i]]$posMotors[, "x"]
    mergedMotors[ istart:istop , "y"] <- img_lst[[i]]$posMotors[, "y"]
    
    istart <- istop + 1 
  }
  PkMat <- FormatPeakMatrix(PkMat, mergedPos,  mergedNumPixels, mergedNames, mergedUUIDs, mergedMotors) 
  
  return(PkMat)
}
