#########################################################################
#     rMSI2 - R package for MSI data handling and visualization
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

#' plotIsotopicPattern
#' 
#' Plots the mean spectra of the region close to an isotopic pattern produced by rMSIannotation
#'
#' @param peakMatrix An rMSIprocPeakMatrix.
#' @param rMSIannotationObj The complete rMSIannotation object returned by rMSIannotation.
#' @param isotopicPattern An element of the list returned by rMSIannotation in the slot isotopicPatterns.
#' @param onlyInPattern If TRUE only the peaks in the pattern will be plotted. Otherwise, all the peaks in the mass range will be plotted
#'
#'
plotIsotopicPattern <- function(peakMatrix, rMSIannotationObj, isotopicPattern, onlyInPattern = TRUE)
{
  
  if(!onlyInPattern){
    plotData <- data.frame(mass = peakMatrix$mass[isotopicPattern[1,2]:isotopicPattern[nrow(isotopicPattern),2]],
                           index = isotopicPattern[1,2]:isotopicPattern[nrow(isotopicPattern),2],
                           intensity = apply(peakMatrix$intensity[,isotopicPattern[1,2]:isotopicPattern[nrow(isotopicPattern),2]], 2, mean),
                           color = "unknown",
                           size = "other")
    plotData$size[match(isotopicPattern[,2],isotopicPattern[1,2]:isotopicPattern[nrow(isotopicPattern),2])] <- "in pattern"
    plotData$color[which(!is.na(match(plotData$index,rMSIannotationObj$monoisotopicPeaks)))] <- "M+0"
    plotData$color[which(!is.na(match(plotData$index,rMSIannotationObj$isotopicPeaks)))] <- "M+N"
    g <- ggplot2::ggplot(data = plotData) + ggplot2::geom_segment(mapping = ggplot2::aes(x = mass, y = intensity,
                                                                                         xend = mass, yend = 0, col = color, linetype = size), size = 1) + 
      ggplot2::theme_classic() + ggplot2::scale_y_continuous(expand = c(0, 0)) + ggplot2::labs(x = "m/z", y = "mean intensity", color = "ion", linetype = "") 
  }
  else {
    plotData <- data.frame(mass = peakMatrix$mass[isotopicPattern[,2]],
                           index = isotopicPattern[,2],
                           intensity = apply(peakMatrix$intensity[,isotopicPattern[,2]], 2, mean),
                           color = "")
    plotData$color[which(!is.na(match(plotData$index, rMSIannotationObj$monoisotopicPeaks)))] <- "M+0"
    plotData$color[which(!is.na(match(plotData$index, rMSIannotationObj$isotopicPeaks)))] <- "M+N"
    g <- ggplot2::ggplot(data = plotData) + ggplot2::geom_segment(mapping = ggplot2::aes(x = mass, y = intensity, xend = mass, yend = 0, col = color), size = 1) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme_classic() + ggplot2::labs(x = "m/z", y = "mean intensity", color = "ion") 
  }
  print(g)
}



#' plotAnnotatedSpectra
#' 
#'Interactive visualization of the mean spectra of the peak matrix with color labaels for the monoisotopic and isotopic peaks. 
#'
#' @param peakMatrix  An rMSIprocPeakMatrix.
#' @param rMSIannotationObj Output of rMSIannotation
#'
plotAnnotatedSpectra <- function(peakMatrix, rMSIannotationObj)
{
  if (!requireNamespace("plotly", quietly = TRUE)) 
  {
    stop("Package plotly needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  candidates <- c()
  for(j in 1:length(rMSIannotationObj$isotopicTestData[[1]]))
  {
    candidates<-c(candidates,rMSIannotationObj$isotopicTestData[[1]][[j]][3])
  }
  
  colors<-rep("unknown",times=length(peakMatrix$mass))
  colors[candidates]<-"M+0 candidate"
  colors[rMSIannotationObj$monoisotopicPeaks]<-"M+0"
  colors[rMSIannotationObj$isotopicPeaks]<-"M+N"
  
  texts<-rep(" ",times=length(peakMatrix$mass))
  for(i in 1:length(rMSIannotationObj$isotopicTestData))
  {
    for(j in 1:length(rMSIannotationObj$isotopicTestData[[i]]))
    {
      candidate<-which.max(rMSIannotationObj$isotopicTestData[[i]][[j]][5,])
      texts[rMSIannotationObj$isotopicTestData[[i]][[j]][4,candidate]]<-format(rMSIannotationObj$isotopicTestData[[i]][[j]][5,candidate],nsmall = 3, digits = 3)
    }
  }
  
  
  fig <- plotly::plot_ly(x=peakMatrix$mass,
                         y=(apply(peakMatrix$intensity,2,mean)),
                         type="bar",
                         color=colors,
                         colors = c("#0CFF00","#FF2D00","#00F0FF","#000000"))
  
  #fig <- fig %>% add_annotations(text=texts[-which(texts==" ")],
  # x=pks$mass[-which(texts==" ")],
  # y=(apply(pks$intensity,2,mean)/max(apply(pks$intensity,2,mean)))[-which(texts==" ")],
  # showarrow=FALSE,
  # font=list(color="black",size=6)) %>% layout(scene = list(yaxis = list(range = c(0,1),rangemode = "tozero")))
  return(fig)
}



#' adductNetworkAnalysis
#' 
#' This function reduces the adducts in tables A and/or B into clusters using network analysis.
#'
#' @param peakMatrix An rMSIprocPeakMatrix.
#' @param rMSIannotationObj Output of rMSIannotation
#' @param adductGroup Flag to indicate from which group adduct should be used. Allowed values are "A", "B" and "AB". Refer to rMSIannotation documentation for details on each group.
#' @param correlationThreshold Only edges with correlation higher than the threshold will be considered for the analysis.
#' 
#' @return List containing all the adduct networks found, solved or not.
#'
adductNetworkAnalysis <- function(peakMatrix, rMSIannotationObj, adductGroup = "A", correlationThreshold = 0.75)
{
  if (!requireNamespace("igraph", quietly = TRUE)) 
  {
    stop("Package igraph needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("tidygraph", quietly = TRUE)) 
  {
    stop("Package tidygraph needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  targetAdducts <- rMSIannotationObj$parameters$adduct_ions
  targetAdducts$name <- paste0("[M",targetAdducts$name,"]")
  # Preparing data structures for the network analysis
  if(adductGroup == "A")
  {
    if(!any(names(rMSIannotationObj)=="A"))
    {
      stop("Group A annotations not found in the rMSIannotation object.",
           call. = FALSE)
    }
    rMSIannotationObj$A$Correlation[rMSIannotationObj$A$Correlation<0] <- 0.01
    relations <- data.frame(from = rMSIannotationObj$A$Adduct1Mass, 
                            to = rMSIannotationObj$A$Adduct2Mass,
                            index1 = rMSIannotationObj$A$Adduct1Index,
                            index2 = rMSIannotationObj$A$Adduct2Index,
                            adductLink = as.factor(rMSIannotationObj$A$Adducts),
                            adduct1 = unlist(lapply(strsplit(rMSIannotationObj$A$Adducts,split = " & "), function(x) return(x[1]))), 
                            adduct2 = unlist(lapply(strsplit(rMSIannotationObj$A$Adducts,split = " & "), function(x) return(x[2]))),
                            correlation = rMSIannotationObj$A$Correlation,
                            priority = rMSIannotationObj$A$AdductPriority,
                            neutralMass = rMSIannotationObj$A$NeutralMass)
    vertices <- unique(c(rMSIannotationObj$A$Adduct1Mass,rMSIannotationObj$A$Adduct2Mass))
  }
  
  if(adductGroup == "B")
  {
    if(!any(names(rMSIannotationObj)=="B"))
    {
      stop("Group B annotations not found in the rMSIannotation object.",
           call. = FALSE)
    }
    cor <- rMSIannotationObj$B$Correlation
    cor[cor<0] <- 0.01
    relations <- data.frame(from = rMSIannotationObj$B$Adduct1Mass, 
                            to = rMSIannotationObj$B$Adduct2Mass,
                            index1 = rMSIannotationObj$B$Adduct1Index,
                            index2 = rMSIannotationObj$B$Adduct2Index,
                            adductLink = as.factor(rMSIannotationObj$B$Adducts),
                            adduct1 = unlist(lapply(strsplit(rMSIannotationObj$B$Adducts,split = " & "), function(x) return(x[1]))), 
                            adduct2 = unlist(lapply(strsplit(rMSIannotationObj$B$Adducts,split = " & "), function(x) return(x[2]))),
                            correlation = cor,
                            priority = rMSIannotationObj$B$AdductPriority,
                            neutralMass = rMSIannotationObj$B$NeutralMass)
    vertices <- unique(c(rMSIannotationObj$B$Adduct1Mass,rMSIannotationObj$B$Adduct2Mass))
  }
  
  if(adductGroup == "AB")
  {
    if((!any(names(rMSIannotationObj)=="A")) | (!any(names(rMSIannotationObj)=="B")))
    {
      stop("Group A or B annotations not found in the rMSIannotation object.",
           call. = FALSE)
    }
    cor <- c(rMSIannotationObj$A$Correlation, rMSIannotationObj$B$Correlation)
    cor[cor<0] <- 0.01
    relations <- data.frame(from = c(rMSIannotationObj$A$Adduct1Mass, rMSIannotationObj$B$Adduct1Mass), 
                            to = c(rMSIannotationObj$A$Adduct2Mass, rMSIannotationObj$B$Adduct2Mass),
                            index1 = c(rMSIannotationObj$A$Adduct1Index,rMSIannotationObj$B$Adduct1Index),
                            index2 = c(rMSIannotationObj$A$Adduct2Index,rMSIannotationObj$B$Adduct1Index),
                            adductLink = as.factor(c(rMSIannotationObj$A$Adducts, rMSIannotationObj$B$Adducts)),
                            adduct1 = unlist(lapply(strsplit(c(rMSIannotationObj$A$Adducts, rMSIannotationObj$B$Adducts),split = " & "), function(x) return(x[1]))), 
                            adduct2 = unlist(lapply(strsplit(c(rMSIannotationObj$A$Adducts, rMSIannotationObj$B$Adducts),split = " & "), function(x) return(x[2]))),
                            correlation = cor,
                            priority = c(rMSIannotationObj$A$AdductPriority, priority = rMSIannotationObj$B$AdductPriority),
                            neutralMass = c(rMSIannotationObj$A$NeutralMass, rMSIannotationObj$B$NeutralMass))
    vertices <- unique(c(rMSIannotationObj$A$Adduct1Mass, rMSIannotationObj$B$Adduct1Mass,
                         rMSIannotationObj$A$Adduct2Mass, rMSIannotationObj$B$Adduct2Mass))
  }
  
  #Network analysis using igraph
  N <- igraph::decompose(igraph::graph_from_data_frame(d = relations, directed=T, vertices=vertices))
  adductNetworks <- list()
  
  #Processing each subgraph
  for(i in 1:length(N))
  {
    g <- igraph::as_data_frame(N[[i]])
    g <- g[g$correlation >= correlationThreshold,]
    if(nrow(g) == 0)
    {
      next
    }
    #Vertices Data
    ionsLabels <- unique(c(g$from,g$to))
    g$v1 <- 0
    g$v2 <- 0
    for(j in 1:length(ionsLabels))
    {
      g$v1[which(g$from == ionsLabels[j])] <- j
      g$v2[which(g$to == ionsLabels[j])] <- j
    }
    
    #Edges Data
    edgeLabels <- unique(paste0(g$v1,g$v2))
    g$edge <- 0
    for(j in 1:length(edgeLabels))
    {
      g$edge[which(paste0(g$v1,g$v2) == edgeLabels[j])] <- j
    }
    
    #Candidates Network
    candidatesNetwork <- igraph::graph_from_data_frame(g,directed = F)
    numberOfEdges <- igraph::gsize(candidatesNetwork)
    
    #Simplified network
    simpNet <- igraph::simplify(candidatesNetwork)
    numberOfEdgesSimp <- igraph::gsize(simpNet)
    
    #Matrix with all the adduct combinations and the priority
    ionsInNetwork <- length(igraph::V(candidatesNetwork))
    edgesNumberVector <- as.numeric(igraph::get.adjacency(candidatesNetwork,type = "upper"))
    numbAddCombinations <- prod(edgesNumberVector[edgesNumberVector!=0])
    adductCombination <- edgesToCombinatios(g, candidatesNetwork, numbAddCombinations)
    adductCombination$possible <- 1
    adductCombination$correctEdges <- 0
    adductCombinationEdgeVerification <- matrix(data = 0, ncol = numberOfEdgesSimp,nrow = numbAddCombinations)
    
    #Verifying which combinations are possible (it will require the capacity of breaking edges...)
    combinationTest <- list()
    for(j in 1:nrow(adductCombination))
    {
      combinationTest[[j]] <- rep(NA, times = ionsInNetwork)
      for(k in 1:numberOfEdgesSimp)
      {
        adductsInEdge <- unlist(strsplit(adductCombination[j,k],split = " & "))
        v1 <- g$v1[which(g$edge==k)][1]
        v2 <- g$v2[which(g$edge==k)][1]
        
        if(is.na(combinationTest[[j]][v1]) & is.na(combinationTest[[j]][v2])) {
          combinationTest[[j]][v1] <- adductsInEdge[1]
          combinationTest[[j]][v2] <- adductsInEdge[2]
        } else if(is.na(combinationTest[[j]][v1])) {
          combinationTest[[j]][v1] <- adductsInEdge[1]
        } else if(is.na(combinationTest[[j]][v2])) {
          combinationTest[[j]][v2] <- adductsInEdge[2]
        } 
        
        if((combinationTest[[j]][v1] == adductsInEdge[1]) & (combinationTest[[j]][v2] == adductsInEdge[2])) {
          adductCombinationEdgeVerification[j, k] <- 1
          adductCombination$correctEdges[j] <- adductCombination$correctEdges[j] + 1
        } else {
          adductCombinationEdgeVerification[j, k] <- 0
          adductCombination$possible[j] <- 0
        }
      }
    }
    
    if(any(adductCombination$possible==1))
    {
      topCorrectEdge <- which(adductCombination$correctEdges==max(adductCombination$correctEdges))
      if(length(topCorrectEdge) > 1)
      {
        topPriorityEdge <- which(adductCombination$priority[topCorrectEdge]==min(adductCombination$priority[topCorrectEdge]))
        optimalSolution <- topPriorityEdge
        ionsTable <- data.frame(mz = unique(c(g$from,g$to)), 
                                adduct = combinationTest[[optimalSolution]])
      }else
      {
        optimalSolution <- topCorrectEdge
        ionsTable <- data.frame(mz = unique(c(g$from,g$to)), 
                                adduct = combinationTest[[optimalSolution]])
      }
      
      #Generating the best graph
      numberOfMultipleEdges<- as.vector(table(g$edge))
      edgeCounter <- rep(1, times = length(numberOfMultipleEdges))
      if(optimalSolution != 1)
      {
        for(j in 1:(optimalSolution-1))
        {
          edgeCounter <- nextCombination(edgeCounter, numberOfMultipleEdges)
        }
      }
      edgeRows <- list()
      for(j in 1:numberOfEdgesSimp)
      {
        edgeRows[[j]] <- which(g$edge == j)
      }
      solvedG <- data.frame(from = 0,to = 0,index1 = 0,index2 = 0,adductLink = 0,
                            adduct1 = 0,adduct2 = 0,correlation = 0,priority = 0,
                            neutralMass = 0,v1 = 0,v2 = 0,edge = 0)
      for(j in 1:numberOfEdgesSimp)
      {
        solvedG[j,] <- g[edgeRows[[j]][edgeCounter[j]],]
      }
      
      adductNetworks[[i]] <- list()
      class(adductNetworks[[i]]) <- "AdductNetwork"
      adductNetworks[[i]]$completNetwork <- candidatesNetwork
      solvedNetwork <- igraph::graph_from_data_frame(solvedG,directed = F)
      adductNetworks[[i]]$edgeLists <- list(complete = g, solved = solvedG)
      adductNetworks[[i]]$ionsTable <- ionsTable
      adductNetworks[[i]]$solvedNetwork <- solvedNetwork
      adductNetworks[[i]]$adductCombination <- adductCombination
    } else
    {
      adductNetworks[[i]] <- list()
      class(adductNetworks[[i]]) <- "AdductNetwork"
      adductNetworks[[i]]$completNetwork <- candidatesNetwork
      adductNetworks[[i]]$edgeLists <- list(complete = g)
      adductNetworks[[i]]$adductCombination <- adductCombination
    }
    #TODO Make the weights in the adducts data frame be consecutive integer numbers (do not allow 0,10,123 but 0,1,2)
  }
  nullLists <- which(unlist(lapply(adductNetworks, is.null)))
  if(length(nullLists) >= 1)
  {
    adductNetworks <- adductNetworks[-nullLists]
  }
  return(adductNetworks)
}

#' plotAdductNetwork
#' 
#' Visualization of an adduct network
#'
#' @param adductNetwork One of the elements returned by the adduct network analysis.
#' @param network  String indicating which network to visualize: "complete", "solved"  or "reduced". 
#'
#' 
plotAdductNetwork <- function(adductNetwork, network = "solved")
{
  if(class(adductNetwork) == "AdductNetwork")
  {
    if(network == "complete")
    {
      igraph::E(adductNetwork$completNetwork)$color <- as.factor(adductNetwork$edgeLists$complete$priority)
      igraph::E(adductNetwork$completNetwork)$label <- paste0("cor:",as.character(adductNetwork$edgeLists$complete$correlation)," \n",
                                                              "V",as.character(adductNetwork$edgeLists$complete$v1),":",as.character(adductNetwork$edgeLists$complete$adduct1),"\n",
                                                              "V",as.character(adductNetwork$edgeLists$complete$v2),":",as.character(adductNetwork$edgeLists$complete$adduct2))
      igraph::V(adductNetwork$completNetwork)$label <- paste0(unique(c(adductNetwork$edgeLists$complete$from,adductNetwork$edgeLists$complete$to)),"\n",
                                                              "V",unique(c(adductNetwork$edgeLists$complete$v1,adductNetwork$edgeLists$complete$v2)))
      igraph::E(adductNetwork$completNetwork)$label.cex <- .8
      return(plot(adductNetwork$completNetwork, layout= igraph::layout.circle, main = "complete network"))
    } 
    
    
    if(network == "solved")
    {
      if(any(names(adductNetwork)=="solvedNetwork")) {
        igraph::E(adductNetwork$solvedNetwork)$color <- as.factor(adductNetwork$edgeLists$solved$priority)
        igraph::E(adductNetwork$solvedNetwork)$label <- paste0("cor:",as.character(adductNetwork$edgeLists$solved$correlation))
        igraph::V(adductNetwork$solvedNetwork)$label <- paste0(adductNetwork$ionsTable$mz, "\n",adductNetwork$ionsTable$adduct)
        return(plot(adductNetwork$solvedNetwork, layout= igraph::layout.circle, main = "solved network"))
      } else  {
        warning("No solved network found for this instance. Plotting the complete graph instead...")
        igraph::E(adductNetwork$completNetwork)$color <- as.factor(adductNetwork$edgeLists$complete$priority)
        igraph::E(adductNetwork$completNetwork)$label <- paste0("cor:",as.character(adductNetwork$edgeLists$complete$correlation)," \n",
                                                                "V",as.character(adductNetwork$edgeLists$complete$v1),":",as.character(adductNetwork$edgeLists$complete$adduct1),"\n",
                                                                "V",as.character(adductNetwork$edgeLists$complete$v2),":",as.character(adductNetwork$edgeLists$complete$adduct2))
        igraph::V(adductNetwork$completNetwork)$label <- paste0(unique(c(adductNetwork$edgeLists$complete$from,adductNetwork$edgeLists$complete$to)),"\n",
                                                                "V",unique(c(adductNetwork$edgeLists$complete$v1,adductNetwork$edgeLists$complete$v2)))
        igraph::E(adductNetwork$completNetwork)$label.cex <- .8
        return(plot(adductNetwork$completNetwork, layout= igraph::layout.circle, main = "complete network"))
      }
    }
    
    
    if(network == "reduced") {
      if(!any(names(adductNetwork)=="reducedNetwork")) {
        warning("No reduced network found for this instance. Plotting the complete graph instead...")
        igraph::E(adductNetwork$completNetwork)$color <- as.factor(adductNetwork$edgeLists$complete$priority)
        igraph::E(adductNetwork$completNetwork)$label <- paste0("cor:",as.character(adductNetwork$edgeLists$complete$correlation)," \n",
                                                                "V",as.character(adductNetwork$edgeLists$complete$v1),":",as.character(adductNetwork$edgeLists$complete$adduct1),"\n",
                                                                "V",as.character(adductNetwork$edgeLists$complete$v2),":",as.character(adductNetwork$edgeLists$complete$adduct2))
        igraph::V(adductNetwork$completNetwork)$label <- paste0(unique(c(adductNetwork$edgeLists$complete$from,adductNetwork$edgeLists$complete$to)),"\n",
                                                                "V",unique(c(adductNetwork$edgeLists$complete$v1,adductNetwork$edgeLists$complete$v2)))
        igraph::E(adductNetwork$completNetwork)$label.cex <- .8
        return(plot(adductNetwork$completNetwork, layout= igraph::layout.circle, main = "complete network"))
      } else {
        igraph::E(adductNetwork$reducedNetwork)$color <- as.factor(adductNetwork$edgeLists$reduced$priority)
        igraph::E(adductNetwork$reducedNetwork)$label <- paste0("cor:",as.character(adductNetwork$edgeLists$reduced$correlation)," \n",
                                                                "V",as.character(adductNetwork$edgeLists$reduced$v1),":",as.character(adductNetwork$edgeLists$reduced$adduct1),"\n",
                                                                "V",as.character(adductNetwork$edgeLists$reduced$v2),":",as.character(adductNetwork$edgeLists$reduced$adduct2))
        igraph::V(adductNetwork$reducedNetwork)$label <- paste0(unique(c(adductNetwork$edgeLists$reduced$from,adductNetwork$edgeLists$reduced$to)),"\n",
                                                                "V",unique(c(adductNetwork$edgeLists$reduced$v1,adductNetwork$edgeLists$reduced$v2)))
        igraph::E(adductNetwork$reducedNetwork)$label.cex <- .8
        return(plot(adductNetwork$reducedNetwork, layout= igraph::layout.circle, main = "reduced network"))
      }
    }
    
    
  } else {
    stop("Input is not an 'AdductNetwork' class instance. This function only accepts elements in the list returned by the function 'adductNetworkAnalysis")
  }
}

#' edgesToCombinatios
#' 
#' 
#'
edgesToCombinatios <- function(g, candidatesNetwork, numbAddCombinations)
{
  numIons <- length(igraph::V(candidatesNetwork)) 
  numEdges <- length(unique(g$edge))
  numberOfMultipleEdges<- as.vector(table(g$edge)) #Vector with the number of repeating edges per edge
  edgeCounter <- rep(1, times = length(numberOfMultipleEdges)) #Vector to know the actual combination. When this vector is equal to numberOfMultipleEdges then stop
  
  #g <- g[unlist(edgeRows[order(numberOfMultipleEdges)]),]
  #numberOfMultipleEdges <- numberOfMultipleEdges[order(numberOfMultipleEdges)]
  edgeRows <- list()
  for(j in 1:numEdges)
  {
    edgeRows[[j]] <- which(g$edge == j)
  }
  
  adductCombination <- matrix(data = 0, ncol = numEdges+1, nrow = numbAddCombinations)
  adductCombination <- as.data.frame(adductCombination)
  
  for(i in 1:numbAddCombinations) 
  {
    for(j in 1:numEdges)
    {
      adductCombination[i,j] <- g$adductLink[edgeRows[[j]][edgeCounter[j]]]
      adductCombination[i,numEdges+1] <- adductCombination[i,numEdges+1] + g$priority[edgeRows[[j]][edgeCounter[j]]]
    }
    edgeCounter <- nextCombination(edgeCounter, numberOfMultipleEdges)
  }
  colnames(adductCombination) <- c(paste0("edge",1:numEdges), "priority")
  return(adductCombination)
}

#' nextCombination
#' This functions computes the following permutation of edges to obtain all the possible adduct pairs combinations in the graph
#'
nextCombination <- function(edgeCounter, numberOfMultipleEdges)
{
  for(n in 1:length(numberOfMultipleEdges))
  {
    if(!((numberOfMultipleEdges[n]-edgeCounter[n]) == 0))
    {
      edgeCounter[n] <- edgeCounter[n] + 1
      return(edgeCounter)
    }
    else if(n != length(numberOfMultipleEdges))
    {
      edgeCounter[1:n] <- 1
    }
  }
  return(numberOfMultipleEdges)
}

#' truncateAndSolveAdductNetwork
#'
#' Manually truncates vertices in adduct networks to solve the network. This is a manual approach, you should look for vertices that are not likly to be part of the network.
#'
#' @param adductNetwork An adductNetwork 
#' @param vertices One or more integers refering to the edge to remove. The integer values of a network can be seen plotting the adductNetwork under analysis.
#'
#' @return An adductNetwork without the vertices removed.
#'
truncateAndSolveAdductNetwork <- function(adductNetwork, vertices)
{
  verticeRows <- which(adductNetwork$edgeLists$complete$v1 %in% vertices | (adductNetwork$edgeLists$complete$v2 %in% vertices))
  adductNetwork$edgeLists$reduced <- adductNetwork$edgeLists$complete[-verticeRows,]
  
  if(nrow(adductNetwork$edgeLists$reduced) != 0)
  {
    cnt <- 1
    for(i in unique(adductNetwork$edgeLists$reduced$edge)){
      adductNetwork$edgeLists$reduced$edge[which(adductNetwork$edgeLists$reduced$edge == i)] <- cnt
      cnt <- cnt + 1
    }
    
    adductNetwork$reducedNetwork <- igraph::graph_from_data_frame(adductNetwork$edgeLists$reduced,directed = T)
    plotAdductNetwork(adductNetwork,network = "reduced")
    {
      #Simplified network
      simpNet <- igraph::simplify(adductNetwork$reducedNetwork)
      numberOfEdgesSimp <- igraph::gsize(simpNet)
      
      #Matrix with all the adduct combinations and the priority
      ionsInNetwork <- length(igraph::V(adductNetwork$reducedNetwork))
      edgesNumberVector <- as.numeric(igraph::get.adjacency(adductNetwork$reducedNetwork,type = "upper"))
      numbAddCombinations <- prod(edgesNumberVector[edgesNumberVector!=0])
      adductCombination <- edgesToCombinatios(adductNetwork$edgeLists$reduced, adductNetwork$reducedNetwork, numbAddCombinations)
      adductCombination$possible <- 1
      adductCombination$correctEdges <- 0
      adductCombinationEdgeVerification <- matrix(data = 0, ncol = numberOfEdgesSimp,nrow = numbAddCombinations)
      
      #Verifying which combinations are possible (it will require the capacity of breaking edges...)
      combinationTest <- list()
      for(j in 1:nrow(adductCombination))
      {
        combinationTest[[j]] <- rep(NA, times = ionsInNetwork)
        for(k in 1:numberOfEdgesSimp)
        {
          adductsInEdge <- unlist(strsplit(adductCombination[j,k],split = " & "))
          v1 <- adductNetwork$edgeLists$reduced$v1[which(adductNetwork$edgeLists$reduced$edge==k)][1]
          v2 <- adductNetwork$edgeLists$reduced$v2[which(adductNetwork$edgeLists$reduced$edge==k)][1]
          
          if(is.na(combinationTest[[j]][v1]) & is.na(combinationTest[[j]][v2])) {
            combinationTest[[j]][v1] <- adductsInEdge[1]
            combinationTest[[j]][v2] <- adductsInEdge[2]
          } else if(is.na(combinationTest[[j]][v1])) {
            combinationTest[[j]][v1] <- adductsInEdge[1]
          } else if(is.na(combinationTest[[j]][v2])) {
            combinationTest[[j]][v2] <- adductsInEdge[2]
          } 
          
          if((combinationTest[[j]][v1] == adductsInEdge[1]) & (combinationTest[[j]][v2] == adductsInEdge[2])) {
            adductCombinationEdgeVerification[j, k] <- 1
            adductCombination$correctEdges[j] <- adductCombination$correctEdges[j] + 1
          } else {
            adductCombinationEdgeVerification[j, k] <- 0
            adductCombination$possible[j] <- 0
          }
        }
      }
      
      if(any(adductCombination$possible==1))
      {
        topCorrectEdge <- which(adductCombination$correctEdges==max(adductCombination$correctEdges))
        if(length(topCorrectEdge) > 1)
        {
          topPriorityEdge <- which(adductCombination$priority[topCorrectEdge]==min(adductCombination$priority[topCorrectEdge]))
          optimalSolution <- topPriorityEdge
          ionsTable <- data.frame(mz = unique(c(adductNetwork$edgeLists$reduced$from,adductNetwork$edgeLists$reduced$to)), 
                                  adduct = combinationTest[[optimalSolution]])
        }else
        {
          optimalSolution <- topCorrectEdge
          ionsTable <- data.frame(mz = unique(c(adductNetwork$edgeLists$reduced$from,adductNetwork$edgeLists$reduced$to)), 
                                  adduct = combinationTest[[optimalSolution]])
        }
        
        #Generating the best graph
        numberOfMultipleEdges<- as.vector(table(adductNetwork$edgeLists$reduced$edge))
        edgeCounter <- rep(1, times = length(numberOfMultipleEdges))
        if(optimalSolution != 1)
        {
          for(j in 1:(optimalSolution-1))
          {
            edgeCounter <- nextCombination(edgeCounter, numberOfMultipleEdges)
          }
        }
        edgeRows <- list()
        for(j in 1:numberOfEdgesSimp)
        {
          edgeRows[[j]] <- which(adductNetwork$edgeLists$reduced$edge == j)
        }
        solvedG <- data.frame(from = 0,to = 0,index1 = 0,index2 = 0,adductLink = 0,
                              adduct1 = 0,adduct2 = 0,correlation = 0,priority = 0,
                              neutralMass = 0,v1 = 0,v2 = 0,edge = 0)
        for(j in 1:numberOfEdgesSimp)
        {
          solvedG[j,] <- adductNetwork$edgeLists$reduced[edgeRows[[j]][edgeCounter[j]],]
        }
        
        adductNetwork$edgeLists$solved <- solvedG
        adductNetwork$ionsTable <- ionsTable
        adductNetwork$solvedNetwork <- igraph::graph_from_data_frame(solvedG,directed = F)
        adductNetwork$adductCombinationReduced <- adductCombination
      } else
      {
        stop("No solution found for this network operation")
      }
    }
  }
  return(adductNetwork)
}
