#########################################################################
#     
#     Copyright (C) 2022 Lluc Sementé Fernández
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

#' peakAnnotation
#' 
#' Peak annotation tool based on the rMSIannotation algorithm Searches for carbon-based isotopic patterns in the peak matrix and evaluates them using morphology, intensity and mass error criteria. Later,
#' the algorithm searches adduct groups using the isotopic patterns.
#'
#' Complete information on rMSIannotation can be found on the paper:
#' rMSIannotation: A peak annotation tool for mass spectrometry imaging based on the analysis of isotopic intensity ratios
#' Lluc Sementé, Gerard Baquer, María García-Altares, Xavier Correig-Blanchar, Pere Ràfols.
#' 
#' DOI: https://doi.org/10.1016/j.aca.2021.338669
#' 
#' @param rMSIprocPeakMatrix An rMSIprocPeakMatrix object.
#' @param params Processing parameters object of the ProcParams class. 
#' 
#' @export
#' @examples  
#' #First you need to load a peak matrix
#' pks <- rMSI2::LoadPeakMatrix("path/to/matrix.pkmat")
#' 
#' #Generate the results with default parameters
#' params <- rMSI2::ProcessingParameters() # Creates a default rMSI2 ProcParams object.
#' annotations <- rMSIproc::peakAnnotation(pks, params)
#' 
#' #Or with custom parameters
#' params <- rMSI2::rMSI2::LoadProcParams("path/to/parameters.params")
#' annotations <- rMSIproc::peakAnnotation(pks, params)
#' 
peakAnnotation <- function(rMSIprocPeakMatrix, params) 
{
  if (!is(rMSIprocPeakMatrix, "rMSIprocPeakMatrix") )
  {
    stop("The provided files does not contain a valid rMSIprocPeakMatrix object\n")
  }
  
  if (!is(params, "ProcParams") )
  {
    stop("The provided files does not contain a valid ProcParams object\n")
  }
  
  rMSIannotation_C(nrow(rMSIprocPeakMatrix$intensity),
                 ncol(rMSIprocPeakMatrix$intensity),
                 50,
                 params$peakAnnotation$ppmMassTolerance,
                 params$peakAnnotation$isotopeLikelihoodScoreThreshold,
                 FALSE,
                 params$peakAnnotation$adductElementsTable,
                 rMSIprocPeakMatrix$intensity,
                 rMSIprocPeakMatrix$mass,
                 1:50,
                 1:50)

  return("Successful Run!")
}


