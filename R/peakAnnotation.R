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
  if(class(rMSIprocPeakMatrix) != "rMSIprocPeakMatrix")
  {
    stop("The provided files does not contain a valid rMSIprocPeakMatrix object\n")
  }
  
  if(class(params) != "ProcParams")
  {
    stop("The provided files does not contain a valid ProcParams object\n")
  }
  
  isotopicPatterns <- list(NULL)
  adductObj <- list(NULL)
  
  isotopicPatterns <- isotopeAnnotation(rMSIprocPeakMatrix, params)

  if(length(isotopicPatterns$monoisotopicPeaks) == 0)
  {
    stop("Procedure cancelled. No monoisotopic peaks found.", call. = FALSE)
  }
  
  if(length(isotopicPatterns$monoisotopicPeaks) >= 2)
  {
    adductNetworks <- adductAnnotation(rMSIprocPeakMatrix, params, isotopicPatterns)
  }
  
  results <- list()
  results <- annotationOutputFormat(rMSIprocPeakMatrix, params, isotopicPatterns, adductNetworks)
  return(results)
}

#' isotopeAnnotation
#' 
#' Finds and evaluate isotope candidates for each ion mass.
#' 
isotopeAnnotation <- function(rMSIprocPeakMatrix, params)
{
  result <- list()  
  
  if ((params$peakAnnotation$isotopeLikelihoodScoreThreshold > 1) || (params$peakAnnotation$isotopeLikelihoodScoreThreshold < 0))
  {
    stop("Error: isotopeLikelihoodScoreThreshold in the ProcParams object must be between 0 and 1")
  }
  
  if((params$peakAnnotation$ppmMassTolerance>100))
  {
    writeLines("Mass tolerances close or bigger than 100 ppms may produce worng annotations.")
  }

  ##### Isotope annotation C++ method #####
  result <- C_isotopeAnnotator(length(rMSIprocPeakMatrix$mass),                               #number of mass peaks
                               1,                                                             #number of mass channels
                               sum(rMSIprocPeakMatrix$numPixels),                             #number of pixels
                               20,                                                            #number of isotopes
                               rMSIprocPeakMatrix$intensity,                                  #intensity matrix
                               rMSIprocPeakMatrix$mass,                                       #matrix mass axis
                               c(1,2),                                                        #image mass axis 
                               params$peakAnnotation$ppmMassTolerance,                        #tolerance in ppm or scans
                               params$peakAnnotation$isotopeLikelihoodScoreThreshold,         #score threshold
                               FALSE,                                                         #tolerance in scans units
                               1)                                                             #charges of the ions in the isotopic pattern
  
  ##### Output format #####
  result <- isotopeOutputFormat(result, params$peakAnnotation$isotopeLikelihoodScoreThreshold)
  return(result)
}


#' adductAnnotation
#' 
#' Given the monoisotopic ions found by the isotopes test, founds the possible adduct pairs considering the elements in the adductDataFrame.
#'
adductAnnotation <- function(rMSIprocPeakMatrix, params, isotopicPatterns)
{
  adducts <- list()
  adductDataFrame <- params$peakAnnotation$adductElementsTable[order(params$peakAnnotation$adductElementsTable$mass,decreasing = T),]
  #namber of adduct combinations
  combinations <- 0
  for(i in 1:nrow(adductDataFrame)-1)
  {
    combinations <- combinations + i
  }
  
  #adduct labels for the output
  name1 <- 1
  name2 <- 2
  lim <- nrow(adductDataFrame)
  namesVector <- c()
  firstnameVector <- c()
  secondnameVector <- c()
  firstnamePriority <- c()
  secondnamePriority <- c()
  for(i in 1:combinations)
  {
    namesVector[i] <- paste(adductDataFrame$name[name1]," & ",adductDataFrame$name[name2],sep = "")
    firstnameVector[i] <- as.character(adductDataFrame$name[name1])
    secondnameVector[i] <- as.character(adductDataFrame$name[name2])
    firstnamePriority[i] <- adductDataFrame$priority[name1]
    secondnamePriority[i] <- adductDataFrame$priority[name2]
    name2 <- name2 + 1
    if(name2 > lim)
    {
      name1 <- name1 + 1
      name2 <- name1 + 1
    }
  }
  
  #monoisotopic to list order
  ord <- c()
  for (i in 1:length(isotopicPatterns$monoisotopicPeaks)) 
  {
    ord <- c(ord,which.min(abs(rMSIprocPeakMatrix$mass[sort(isotopicPatterns$monoisotopicPeaks)[i]]-as.numeric(names(isotopicPatterns$M1)))))
  }
  ord <- ord-1
  
  #label axis
  labelAxis <- rep(0, times = length(rMSIprocPeakMatrix$mass))
  labelAxis[isotopicPatterns$monoisotopicPeaks] <- 1
  labelAxis[isotopicPatterns$isotopicPeaks] <- 2
  
  M1isotopes <- isotopicPatterns$M1
  
  for(i in 1:length(M1isotopes))
  {
    M1isotopes[[i]] <- M1isotopes[[i]][,which.max(M1isotopes[[i]][5,])]
  }
  
  ##### Calling the C++ method #####
  adducts <- C_adductAnnotation(length(isotopicPatterns$monoisotopicPeaks),
                                nrow(adductDataFrame), 
                                params$peakAnnotation$ppmMassTolerance,
                                length(rMSIprocPeakMatrix$mass),
                                rMSIprocPeakMatrix$mass[sort(isotopicPatterns$monoisotopicPeaks)],
                                adductDataFrame$mass,
                                M1isotopes,
                                ord,
                                rMSIprocPeakMatrix$mass,
                                rMSIprocPeakMatrix$intensity,
                                sum(rMSIprocPeakMatrix$numPixels),
                                labelAxis,
                                sort(isotopicPatterns$monoisotopicPeaks)-1)   
  
  for(i in 1:2)
  {
    adducts[[i]] <- adducts[[i]][-((length(adducts[[i]])-sum(unlist(lapply(adducts[[i]], is.null)))+1):length(adducts[[i]]))]
  }
  names(adducts) <- c("A","B")
  
  ## Data frame for A quality adducts ##
  if(length(adducts$A) > 0)
  {
    adductsA <- data.frame(
      NeutralMass = rep(0, times = length(adducts$A)),
      Adducts = rep(0, times = length(adducts$A)),
      Adduct1Mass = rep(0, times = length(adducts$A)),
      Adduct2Mass = rep(0, times = length(adducts$A)),
      IsotopeIntensityRatioMean = rep(0, times = length(adducts$A)),
      IsotopeIntensityRatioStdError = rep(0, times = length(adducts$A)),
      Correlation = rep(0, times = length(adducts$A)),
      MassError = rep(0, times = length(adducts$A)),
      Adduct1Index = rep(0, times = length(adducts$A)),
      Adduct2Index = rep(0, times = length(adducts$A)),
      AdductPriority = rep(0, times = length(adducts$A)))
    for(i in 1:length(adducts$A))
    {
      name1 <- firstnameVector[(adducts$A[[i]][1]+1)]
      name2 <- secondnameVector[(adducts$A[[i]][1]+1)]
      if(adducts$A[[i]][3] > adducts$A[[i]][5])
      {
        adductsA$Adducts[i] <- paste("[M",name1,"] & [M",name2,"]",sep="")
      }
      else
      {
        adductsA$Adducts[i] <- paste("[M",name2,"] & [M",name1,"]",sep="")
      }
      adductsA$NeutralMass[i]                   <- adducts$A[[i]][2]
      adductsA$Adduct1Mass[i]                   <- adducts$A[[i]][3]
      adductsA$Adduct1Index[i]                  <- adducts$A[[i]][4]
      adductsA$Adduct2Mass[i]                   <- adducts$A[[i]][5]
      adductsA$Adduct2Index[i]                  <- adducts$A[[i]][6]
      adductsA$IsotopeIntensityRatioMean[i]     <- adducts$A[[i]][7]
      adductsA$IsotopeIntensityRatioStdError[i] <- adducts$A[[i]][8]
      adductsA$Correlation[i]                   <- adducts$A[[i]][9]
      adductsA$MassError[i]                     <- adducts$A[[i]][10] 
      adductsA$AdductPriority[i]                <- (firstnamePriority[(adducts$A[[i]][1]+1)]) + (secondnamePriority[(adducts$A[[i]][1]+1)])
    }
    adductsA$Correlation <- signif(adductsA$Correlation, digits = 3)
    adductsA$MassError <- trunc(adductsA$MassError) + signif(adductsA$MassError - trunc(adductsA$MassError), digits = 3)
    adductsA$Adduct1Mass <- trunc(adductsA$Adduct1Mass) + signif(adductsA$Adduct1Mass - trunc(adductsA$Adduct1Mass), digits = 4)
    adductsA$Adduct2Mass <- trunc(adductsA$Adduct2Mass) + signif(adductsA$Adduct2Mass - trunc(adductsA$Adduct2Mass), digits = 4)
    adductsA$NeutralMass <- trunc(adductsA$NeutralMass) + signif(adductsA$NeutralMass - trunc(adductsA$NeutralMass), digits = 4)
    adductsA$IsotopeIntensityRatioMean <- trunc(adductsA$IsotopeIntensityRatioMean) + signif(adductsA$IsotopeIntensityRatioMean - trunc(adductsA$IsotopeIntensityRatioMean), digits = 3)
    adductsA$IsotopeIntensityRatioStdError <- trunc(adductsA$IsotopeIntensityRatioStdError) + signif(adductsA$IsotopeIntensityRatioStdError - trunc(adductsA$IsotopeIntensityRatioStdError), digits = 4)
    
    adductsA <- adductsA[order(adductsA$NeutralMass,decreasing = F),]
    row.names(adductsA) <- NULL
    adducts$A <- adductsA[order(adductsA$NeutralMass,decreasing = F),]
  }
  else {adducts <- adducts[-1]}
  
  ## Data frame for B quality adducts ##
  if(length(adducts$B) > 0)
  {
    adductsB <- data.frame(
      NeutralMass = rep(0, times = length(adducts$B)),
      Adducts = rep(0, times = length(adducts$B)),
      Adduct1Mass = rep(0, times = length(adducts$B)),
      Adduct2Mass = rep(0, times = length(adducts$B)),
      Correlation = rep(0, times = length(adducts$B)),
      MassError = rep(0, times = length(adducts$B)),
      Adduct1Index = rep(0, times = length(adducts$B)),
      Adduct2Index = rep(0, times = length(adducts$B)))
    for(i in 1:length(adducts$B))
    {
      name1 <- firstnameVector[(adducts$B[[i]][1]+1)]
      name2 <- secondnameVector[(adducts$B[[i]][1]+1)]
      if(adducts$B[[i]][3] > adducts$B[[i]][5])
      {
        adductsB$Adducts[i] <- paste("[M",name1,"] & [M",name2,"]",sep="")
      }
      else
      {
        adductsB$Adducts[i] <- paste("[M",name2,"] & [M",name1,"]",sep="")
      }
      adductsB$NeutralMass[i]    <- adducts$B[[i]][2]
      adductsB$Adduct1Mass[i]    <- adducts$B[[i]][3]
      adductsB$Adduct1Index[i]   <- adducts$B[[i]][4]
      adductsB$Adduct2Mass[i]    <- adducts$B[[i]][5]
      adductsB$Adduct2Index[i]   <- adducts$B[[i]][6]
      adductsB$Correlation[i]    <- adducts$B[[i]][7]
      adductsB$MassError[i]      <- adducts$B[[i]][8]
      adductsB$AdductPriority[i] <- firstnamePriority[(adducts$B[[i]][1]+1)] + secondnamePriority[(adducts$B[[i]][1]+1)]
    }
    adductsB$Correlation <- signif(adductsB$Correlation, digits = 3)
    adductsB$MassError <- trunc(adductsB$MassError) + signif(adductsB$MassError - trunc(adductsB$MassError), digits = 3)
    adductsB$Adduct1Mass <- trunc(adductsB$Adduct1Mass) + signif(adductsB$Adduct1Mass - trunc(adductsB$Adduct1Mass), digits = 4)
    adductsB$Adduct2Mass <- trunc(adductsB$Adduct2Mass) + signif(adductsB$Adduct2Mass - trunc(adductsB$Adduct2Mass), digits = 4)
    adductsB$NeutralMass <- trunc(adductsB$NeutralMass) + signif(adductsB$NeutralMass - trunc(adductsB$NeutralMass), digits = 4)
    
    adductsB <- adductsB[order(adductsB$NeutralMass,decreasing = F),]
    row.names(adductsB) <- NULL
    adducts$B <- adductsB[order(adductsB$NeutralMass,decreasing = F),]
  }
  else {adducts <- adducts[-2]}
  
  return(adducts)
}

#' isotopeOutputFormat
#' 
#' Gives format to the isotope module output. 
#'
isotopeOutputFormat <- function(r, ScoreThreshold)
{
  ##### Names the results list #####
  NullVec <- matrix(nrow = length(r)-1, ncol = length(r[[1]]))
  Nullcnt <- 0
  
  for(i in 2:length(r))
  {
    names(r[[i]]) <- signif(r[[1]],digits = 8)
    Nullcnt <- 0
    for(j in 1:length(r[[i]]))
    {
      if(!is.null(r[[i]][[j]]))
      {
        row.names(r[[i]][[j]]) <- c("M+N mass","Abs ppm error","M+0 index", "M+N index","Final score","Morphology score","Intensity score", "Mass error score", "Isotopic ratio", "Number of C atoms")
        colnames(r[[i]][[j]]) <- (1:(dim(r[[i]][[j]])[2]))
      }else
      {
        Nullcnt <- Nullcnt + 1
        NullVec[i-1,Nullcnt] <- j
      }
    }
    r[[i]] <- r[[i]][-(NullVec[i-1,1:Nullcnt])]
  }
  
  
  ##### Removes the mass vector and the empty isotope lists #####
  r <- r[-1]
  flagRemList <- FALSE
  RemList <- c()
  
  for(i in 1:length(r))
  {
    if(length(r[[i]])==0)
    {
      flagRemList <- TRUE
      RemList <- c(RemList,i)
    }
  }
  
  if(flagRemList)
  {
    r <- r[-RemList]
  }
  
  if(length(r) == 0)
  {  
    cat(paste("Number of monoisotopic ions found: 0\n"))
    return()
  }
  
  ##### Creats the isotopic & monoisotopic ions index vectors #####
  MonoVec <- c()
  CompVec <- c()
  for(i in (1:(length(r))))
  {
    for(j in (1:length(r[[i]])))
    {
      if(i == 1)
      {
        if((!is.null(r[[i]][[j]])))
        {
          if((as.numeric(r[[i]][[j]][5,which.max(r[[i]][[j]][5,])]) >= ScoreThreshold))
          {
            MonoVec <- c(MonoVec, r[[i]][[j]][3,which.max(r[[i]][[j]][5,])])
          }
        }
      }
      
      bestIsotope <- which.max(r[[i]][[j]][5,])
      if(i == 1)
      {
        if(as.numeric(r[[i]][[j]][5, bestIsotope]) >= ScoreThreshold)
        {
          CompVec <- c(CompVec, r[[i]][[j]][4, bestIsotope])
        }
      }
      else
      {
        if((as.numeric(r[[i]][[j]][6, bestIsotope]) >= ScoreThreshold) &
           (as.numeric(r[[i]][[j]][7, bestIsotope]) >= ScoreThreshold) &
           (as.numeric(r[[i]][[j]][8, bestIsotope]) >= ScoreThreshold))
        {
          CompVec <- c(CompVec, r[[i]][[j]][4, bestIsotope])
        }
      }
    }
  }
  
  names(r) <- paste("M",1:length(r),sep = "")
  r$isotopicPeaks <- (unique(CompVec))
  r$monoisotopicPeaks <- (unique(MonoVec))
  r$monoisotopicPeaks <- setdiff(r$monoisotopicPeaks, r$isotopicPeaks)
  MonoVec <- r$monoisotopicPeaks
  cat(paste("Number of monoisotopic ions found:", length(MonoVec),"\n"))
  
  return(r)
}


#' annotationOutputFormat
#' 
#' This function gives format to the output
#'
annotationOutputFormat <- function(rMSIprocPeakMatrix, params, isotopeObj, adductObj)
{
  #Isotope results
  C <- data.frame(MonoisotopicMass = rep(0, times = length(isotopeObj$monoisotopicPeaks)),
                  ILS = rep(0, times = length(isotopeObj$monoisotopicPeaks)),
                  IsotopicIntensityRatio = rep(0, times = length(isotopeObj$monoisotopicPeaks)),
                  EstimatedCarbonAtoms = rep(0, times = length(isotopeObj$monoisotopicPeaks)),
                  MonoisotopicIndex = rep(0, times = length(isotopeObj$monoisotopicPeaks)))
  ord <- c()
  for (i in 1:length(isotopeObj$monoisotopicPeaks)) 
  {
    ord <- c(ord,which.min(abs(rMSIprocPeakMatrix$mass[sort(isotopeObj$monoisotopicPeaks)[i]]-as.numeric(names(isotopeObj$M1)))))
  }
  for(i in 1:length(ord))
  {
    C$MonoisotopicMass[i]       <- rMSIprocPeakMatrix$mass[isotopeObj$M1[[ord[i]]][3]]
    C$MonoisotopicIndex[i]      <- isotopeObj$M1[[ord[i]]][3,which.max(isotopeObj$M1[[ord[i]]][5,])]
    C$ILS[i]                    <- isotopeObj$M1[[ord[i]]][5,which.max(isotopeObj$M1[[ord[i]]][5,])]
    C$IsotopicIntensityRatio[i] <- isotopeObj$M1[[ord[i]]][9,which.max(isotopeObj$M1[[ord[i]]][5,])]
    C$EstimatedCarbonAtoms[i]   <- isotopeObj$M1[[ord[i]]][10,which.max(isotopeObj$M1[[ord[i]]][5,])]
  }
  C$MonoisotopicMass       <- trunc(C$MonoisotopicMass) + signif(C$MonoisotopicMass-trunc(C$MonoisotopicMass), digits = 4) 
  C$MonoisotopicIndex      <- trunc(C$MonoisotopicIndex) + signif(C$MonoisotopicIndex-trunc(C$MonoisotopicIndex), digits = 4)
  C$ILS                    <- trunc(C$ILS) + signif(C$ILS-trunc(C$ILS), digits = 4)
  C$IsotopicIntensityRatio <- trunc(C$IsotopicIntensityRatio) + signif(C$IsotopicIntensityRatio-trunc(C$IsotopicIntensityRatio), digits = 4)
  C$EstimatedCarbonAtoms   <- round(C$EstimatedCarbonAtoms)
  IsotopicTestData <- isotopeObj[-((length(isotopeObj)-1):length(isotopeObj))]
  
  #Isotopic patterns
  Patterns <- list()
  for(i in 1:length(isotopeObj$monoisotopicPeaks))
  {
    monoIsoMass <- rMSIprocPeakMatrix$mass[isotopeObj$monoisotopicPeaks[i]]
    tmp <- list()
    cntPattern <- 0
    for(j in 1:length(IsotopicTestData))
    {
      if(any(names(IsotopicTestData[[j]]) == signif(monoIsoMass,digits = 8)))
      {
        instance <- which(names(IsotopicTestData[[j]]) == signif(monoIsoMass,digits = 8))
        
        if((max(IsotopicTestData[[j]][[instance]][5,]) >= params$peakAnnotation$isotopeLikelihoodScoreThreshold) && 
           (j==1))
        {
          tmp[[j]] <- IsotopicTestData[[j]][[instance]][,which.max(IsotopicTestData[[j]][[instance]][5,])]
          cntPattern <- cntPattern + 1
        }
        
        if((IsotopicTestData[[j]][[instance]][6,which.max(IsotopicTestData[[j]][[instance]][5,])] >= params$peakAnnotation$isotopeLikelihoodScoreThreshold) && 
           (IsotopicTestData[[j]][[instance]][7,which.max(IsotopicTestData[[j]][[instance]][5,])] >= params$peakAnnotation$isotopeLikelihoodScoreThreshold) && 
           (IsotopicTestData[[j]][[instance]][8,which.max(IsotopicTestData[[j]][[instance]][5,])] >= params$peakAnnotation$isotopeLikelihoodScoreThreshold) && 
           (j!=1))
        {
          tmp[[j]] <- IsotopicTestData[[j]][[instance]][,which.max(IsotopicTestData[[j]][[instance]][5,])]
          cntPattern <- cntPattern + 1
        }
        
      }
      else  
      {
        break
      }
    }
    rnames <- paste0("M+",0:cntPattern)
    cnames <- c("M+N m/z", "M+N index", "Mass error (ppm)", "ILS", "Morphology Score", "Intensity Score", "Mass Error Score", "Estimated C atoms", "Intensity ratio")
    if(length(tmp) > 1)
    {
      Patterns[[i]] <- t(as.data.frame(tmp))
      Patterns[[i]] <- rbind(Patterns[[i]][1,], Patterns[[i]])
      Patterns[[i]] <- Patterns[[i]][,c(1,4,2,5,6,7,8,10,9)]
      rownames(Patterns[[i]]) <- rnames
      colnames(Patterns[[i]]) <- cnames
      Patterns[[i]][1,] <- c(monoIsoMass,
                             isotopeObj$monoisotopicPeaks[i], 
                             weighted.mean(Patterns[[i]][-1,3], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,4], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,5], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,6], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,7], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,8], (10^-((1:cntPattern)-1))),
                             weighted.mean(Patterns[[i]][-1,9], (10^-((1:cntPattern)-1))))
      
    }
    else
    {
      Patterns[[i]] <- as.vector(t(as.data.frame(tmp,optional = TRUE)))
      Patterns[[i]] <- rbind(Patterns[[i]][1], Patterns[[i]])
      Patterns[[i]] <- Patterns[[i]][,c(1,4,2,5,6,7,8,10,9)]
      Patterns[[i]][1,] <- c(monoIsoMass,
                             isotopeObj$monoisotopicPeaks[i], 
                             weighted.mean(Patterns[[i]][-1,3], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,4], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,5], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,6], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,7], 1/(1:cntPattern)),
                             weighted.mean(Patterns[[i]][-1,8], (10^-((1:cntPattern)-1))),
                             weighted.mean(Patterns[[i]][-1,9], (10^-((1:cntPattern)-1))))
      rownames(Patterns[[i]]) <- rnames
      colnames(Patterns[[i]]) <- cnames
    }
  }
  names(Patterns) <- signif(rMSIprocPeakMatrix$mass[isotopeObj$monoisotopicPeaks],digits = 8)
  Patterns <- Patterns[order(as.numeric(names(Patterns)))]
  
  #Adducts table
  flagA <- FALSE
  flagB <- FALSE
  if(!is.null(adductObj[[1]]))
  {
    if(!is.null(adductObj$A))
    {
      A <- adductObj$A
      flagA <- TRUE
    }
    if(!is.null(adductObj$B))
    {
      B <- adductObj$B
      flagB <- TRUE
    }
  }
  
  
  if(flagA | flagB)
  {
    if(flagA & flagB)
    {
      adductTables <- rbind(adductObj$A[,c(2,3,4,9,10)],adductObj$B[,c(2,3,4,7,8)])
      adductTables$group <- rep(c("A","B"), times = c(nrow(adductObj$A), nrow(adductObj$B)))
    } else if(flagA)
    {
      adductTables <- adductObj$A[,c(2,3,4,9,10)]
      adductTables$group <- rep("A", times = c(nrow(adductObj$A)))
    } else if(flagB)
    {
      adductTables <- adductObj$B[,c(2,3,4,7,8)]
      adductTables$group <- rep("B", times = c(nrow(adductObj$B)))
    }
    adducts12 <- strsplit(adductTables$Adducts,split = " & ")
    adducts12 <- t(as.data.frame(adducts12))
    adductsAndIndex <- data.frame(adduct = c(adducts12[,1],adducts12[,2]),
                                  index = c(adductTables$Adduct1Index,adductTables$Adduct2Index),
                                  group = rep(adductTables$group, times = 2))
    rm(adductTables, adducts12)
    
    adductNames <- paste0("[M",params$peakAnnotation$adductElementsTable$name,"]")
    
    aT <- data.frame()
    for(i in 1:length(isotopeObj$monoisotopicPeaks))
    {
      aT[1:nrow(params$peakAnnotation$adductElementsTable), i] <- rep("", times =nrow(params$peakAnnotation$adductElementsTable))  
    }
    
    for(i in 1:length(isotopeObj$monoisotopicPeaks))
    {
      peakIndex <- isotopeObj$monoisotopicPeaks[i] #search by index, not mass
      matches <- unique(adductsAndIndex$adduct[which(adductsAndIndex$index == peakIndex)])
      for(j in 1:length(matches))
      {
        aT[which(adductNames==matches[j]), i] <- adductsAndIndex$group[which(adductsAndIndex$index == peakIndex)[j]]
      }
    }
    
    colnames(aT) <- signif(rMSIprocPeakMatrix$mass[isotopeObj$monoisotopicPeaks],digits = 8)
    rownames(aT) <- adductNames
    aT <- aT[,order(rMSIprocPeakMatrix$mass[isotopeObj$monoisotopicPeaks],decreasing = F)]
  }
  
  
  #Final results list
  if(!is.null(adductObj[[1]]))
  {
    if(is.null(adductObj$A))
    {
      results <- list(B = B,C = C, isotopicPatterns = Patterns, isotopicTestData = IsotopicTestData, adductAssignation = aT, monoisotopicPeaks = isotopeObj$monoisotopicPeaks, isotopicPeaks = isotopeObj$isotopicPeaks)
    }
    if(is.null(adductObj$B))
    {
      results <- list(A = A,C = C, isotopicPatterns = Patterns, isotopicTestData = IsotopicTestData, adductAssignation = aT, monoisotopicPeaks = isotopeObj$monoisotopicPeaks, isotopicPeaks = isotopeObj$isotopicPeaks)
    }
    if(!is.null(adductObj$B) & !is.null(adductObj$A))
    {
      results <- list(A = A,B = B,C = C, isotopicPatterns = Patterns, isotopicTestData = IsotopicTestData, adductAssignation = aT, monoisotopicPeaks = isotopeObj$monoisotopicPeaks, isotopicPeaks = isotopeObj$isotopicPeaks)
    }
  }
  else
  {
    results <- list(C = C, isotopicPatterns = Patterns, isotopicTestData = IsotopicTestData, monoisotopicPeaks = isotopeObj$monoisotopicPeaks, isotopicPeaks = isotopeObj$isotopicPeaks)
  }
  return(results)
}

