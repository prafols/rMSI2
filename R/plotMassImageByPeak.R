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

#Idem k anterior xo sense:
# - selectedPixels: Ho solucionare amb terra::crop terra::zoom i limatant el valor maxim de Z en plot del raster per cada layer RGB
# - rotate: Ho solucionare directament rotant objecte de raster en plotar, no em cal re-generar imatge per plotar
# - Retorna: una llista on
#         img.data correspondra al raster de la imatge,
#         mass.peak
#         tolerance
# El resultat d'aquest metode es guardara en la classe MSImagePlotWidget en una llista de capes R,G,B
.buildImageByPeak<-function(img, mass.peak, tolerance=0.25, NormCoefs = NULL)
{
  #Grab only MassPeaks Lists
  if(typeof(img) != "list")
  {
    stop("Error in rMSI2::.buildImageByPeak(), img is not a list")
  }

  #Get the image slice
  img_slice<-builRasterImageFromMass( img, mass.peak, tolerance, NormCoefs)

  #Create the raster
  my_raster <- terra::rast( nrow = ncol(img_slice$pixels), ncol = nrow(img_slice$pixels), xmin= 0, xmax= nrow(img_slice$pixels), ymin= 0, ymax= ncol(img_slice$pixels))
  terra::values(my_raster) <- as.vector(img_slice$pixels)

  #Return zplots matrix and some metadata in a list
  list(raster = my_raster, mass = img_slice$Mass, tolerance = img_slice$Tolerance, cal_resolution = img$pixel_size_um)
}


#Create an empty raster object
.InitRGBEmptyRaster<-function( width, height)
{
  new_raster<-terra::rast( nrow = height, ncol = width, xmin= 0, xmax= width, ymin= 0, ymax= height)
  terra::values(new_raster) <- rep(0, width*height)
  return( list( raster = new_raster, mass = NULL, tolerance = NULL, cal_resolution = NULL))
}


#Normalize to 255 (8 bits per color channel)
.NormalizeTo255 <-function( m, scaling )
{
  maxN<- max(terra::values(m))
  if(maxN > 0)
  {
    terra::values(m) <- scaling*255 * terra::values(m) / (maxN*3)
    terra::values(m)[terra::values(m)>255]<- 255 #Clipping
  }
  return(m)
}

#Normalize by 0 to 1
.NormalizeFrom0to1 <-function( values )
{
  maxN <- max(values)
  minN <- min(values)
  if((maxN - minN) != 0)
  {
    m <- 1/(maxN - minN)
  }
  else
  {
    m <- 0
  }
  return( values * m - m*minN)
}

#Idem k anterior plotMassImageRGB
#Com a param img_RGB es una llista amb atributs R,G,B on cada objecte (layer) es un objecte retornat per builImageByPeak
# XResLevel es ara un integer ja que delego la interpolacio a raster
# All RGB layer must have the same size
.BuildRGBImage <- function( imgR, imgG, imgB, XResLevel = 3 , light = 3)
{
  #Normalize from 0 to 1
  terra::values(imgR$raster)<-.NormalizeFrom0to1(terra::values(imgR$raster))
  terra::values(imgG$raster)<-.NormalizeFrom0to1(terra::values(imgG$raster))
  terra::values(imgB$raster)<-.NormalizeFrom0to1(terra::values(imgB$raster))

  #Normalize to 255 using light
  imgR <- .NormalizeTo255(imgR$raster, light)
  imgG <- .NormalizeTo255(imgG$raster, light)
  imgB <- .NormalizeTo255(imgB$raster, light)

  #Create an RGB image space
  RGB_raster <- imgR
  terra::add(RGB_raster) <- imgG
  terra::add(RGB_raster) <- imgB
  
  interpolated_raster <- terra::rast( nrow= XResLevel*terra::nrow(RGB_raster), ncol= XResLevel*terra::ncol(RGB_raster), xmin= 0, xmax= terra::ncol(RGB_raster), ymin= 0, ymax= terra::nrow(RGB_raster))
  RGB_raster<-terra::resample(RGB_raster, interpolated_raster)
  terra::values(RGB_raster)[ terra::values(RGB_raster) < 0  ] <- 0 #Values below zero are dube interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

#Remap a intensity vector to HSV coloring in a rainbow function
.ReMappingIntensity2HSV<-function(single_channel_raster, maxN = max(terra::values(single_channel_raster)), value_multiplier = 3)
{
  auxVals <- terra::values(single_channel_raster) #Keep for vMApped

  #Normalize to a range 0 -> 1
  minRas <- min(terra::values(single_channel_raster))
  if(( maxN - minRas ) == 0)
  {
    m <- 0
    m_V <- 0
  }
  else
  {
    m <- 1/( maxN - minRas )
    m_V <- 1/( max(c( abs(minRas), abs(maxN) ) ) )
  }
  terra::values(single_channel_raster) <- terra::values(single_channel_raster) * m - m*minRas

  #Remapping hue space
  hue_top <- 0.7
  hue_bottom <- 0.85
  hMapped <- (1 + hue_top - hue_bottom) * (-1*terra::values(single_channel_raster)  + 1) + hue_bottom
  over_one <- which(hMapped > 1)
  hMapped[ over_one ] <- hMapped[over_one] -1

  #Remapping value space
  vMapped <- abs(auxVals*m_V)
  vMapped <- vMapped  * value_multiplier
  vMapped[ vMapped > 1] <- 1

  #Creating HSV image
  rgbSpace <- col2rgb(hsv( h =  hMapped, s = rep(1, length(terra::values(single_channel_raster) )), v = vMapped))

  #Layering RGB raster
  R_raster<- single_channel_raster
  G_raster<- single_channel_raster
  B_raster<- single_channel_raster
  terra::values(R_raster) <-  rgbSpace[1, ] #Extract R
  terra::values(G_raster) <-  rgbSpace[2, ] #Extract G
  terra::values(B_raster) <-  rgbSpace[3, ] #Extract B

  #Create an RGB image space
  RGB_raster <- R_raster
  terra::add(RGB_raster) <- G_raster
  terra::add(RGB_raster) <- B_raster
  return(RGB_raster)
}

#Build a RGB images raster using rainbow colors from only one raster layer
.BuildSingleIonRGBImage<-function( img,   XResLevel = 3 , global_intensity_scaling_factor = NULL, light = 3)
{
  #Create an RGB image space
  if(is.null(global_intensity_scaling_factor))
  {
    RGB_raster <- .ReMappingIntensity2HSV(img$raster, value_multiplier = light)
  }
  else
  {
    RGB_raster <- .ReMappingIntensity2HSV(img$raster, maxN = global_intensity_scaling_factor, value_multiplier = light)
  }
  
  interpolated_raster <- terra::rast( nrow= XResLevel*terra::nrow(RGB_raster), ncol= XResLevel*terra::ncol(RGB_raster), xmin= 0, xmax= terra::ncol(RGB_raster), ymin= 0, ymax= terra::nrow(RGB_raster))
  RGB_raster<-terra::resample(RGB_raster, interpolated_raster)
  terra::values(RGB_raster)[ terra::values(RGB_raster) < 0  ] <- 0 #Values below zero are due interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

#roi_rectangle  is c(left, rigth, bottom, top)
.plotMassImageRGB <- function(rasterRGB, cal_um2pixels = 1,  rotation=0, flipV=F, flipH=F, display_axes=T, display_syscoords=T, roi_rectangle = NULL, zoom_rectangle = NULL, border = 10)
{
  img_Xmax <- terra::ext(rasterRGB)$xmax
  img_Ymax <- terra::ext(rasterRGB)$ymax

  #Crop raster according the zoom window
  if( !is.null(zoom_rectangle) )
  {
    rasterRGB<- terra::crop( rasterRGB, terra::ext( c(zoom_rectangle[1] -1, zoom_rectangle[2],  terra::ext(rasterRGB)$ymax - zoom_rectangle[4],  terra::ext(rasterRGB)$ymax - zoom_rectangle[3] + 1)))
  }

  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.7, mar = c(1,1,1,1), mgp = c(2, 0.5, 0.5))

  #Apply flip
  if((flipV && rotation == 0) || (flipV && rotation == 180) || (flipH && rotation == 90) || (flipH && rotation == 270))
  {
    aux <- img_Ymax - roi_rectangle[4]
    roi_rectangle[4] <- img_Ymax - roi_rectangle[3]
    roi_rectangle[3] <- aux
    rasterRGB <- terra::flip(rasterRGB, direction = "vertical")
  }
  if((flipH && rotation == 0) || (flipH && rotation == 180) || (flipV && rotation == 90) || (flipV && rotation == 270))
  {
    aux <- img_Xmax - roi_rectangle[2]
    roi_rectangle[2] <- img_Xmax - roi_rectangle[1]
    roi_rectangle[1] <- aux
    rasterRGB <- terra::flip(rasterRGB, direction = "horizontal")
  }

  #Apply rotation
  roi_rectangle_pre <- roi_rectangle
  if( rotation == 0 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- roi_rectangle_pre[1] - 1 #Left maps Left
      roi_rectangle[2] <- roi_rectangle_pre[2] #Rigth maps Right
      roi_rectangle[3] <- img_Ymax - roi_rectangle_pre[4] #Bottom  maps Top
      roi_rectangle[4] <- img_Ymax - roi_rectangle_pre[3] + 1 #Top maps Bottom
    }
  }

  if( rotation == 90 )
  {
    if(!is.null(roi_rectangle))
    {
      if(is.null(zoom_rectangle))
      {
        roi_rectangle[1] <-  roi_rectangle_pre[3] - 1 #Left maps Bottom
        roi_rectangle[2] <- roi_rectangle_pre[4] #Rigth maps Top
      }
      else
      {
        roi_rectangle[1] <- img_Ymax - zoom_rectangle[4] - zoom_rectangle[3] + roi_rectangle_pre[3] #Left maps Bottom
        roi_rectangle[2] <- img_Ymax - zoom_rectangle[4] - zoom_rectangle[3] + roi_rectangle_pre[4] +1 #Rigth maps Top
      }
      roi_rectangle[3] <- roi_rectangle_pre[1] - 1 #Bottom  maps Left
      roi_rectangle[4] <- roi_rectangle_pre[2] #Top maps Right
    }


    rasterRGB <- terra::flip(terra::t(rasterRGB), direction = "vertical") #90degree rotation
  }
  if( rotation == 270 )
  {
    if(!is.null(roi_rectangle))
    {
      roi_rectangle[1] <- img_Ymax - roi_rectangle_pre[4] #Left maps top
      roi_rectangle[2] <- img_Ymax - roi_rectangle_pre[3] + 1 #Rigth maps bottom
      if(is.null(zoom_rectangle))
      {
        roi_rectangle[3] <- img_Xmax - roi_rectangle_pre[2] #Bottom  maps Right
        roi_rectangle[4] <- img_Xmax - roi_rectangle_pre[1] + 1 #Top maps Left
      }
      else
      {
        roi_rectangle[3] <- zoom_rectangle[2] + zoom_rectangle[1] - roi_rectangle_pre[2] - 1 #Bottom  maps Right
        roi_rectangle[4] <- zoom_rectangle[2] + zoom_rectangle[1] - roi_rectangle_pre[1] #Top maps Left
      }
    }

    rasterRGB <- terra::flip(terra::t(rasterRGB), direction = "horizontal") #270degree rotation
  }
  if( rotation == 180 )
  {
    if(!is.null(roi_rectangle))
    {
      if(is.null(zoom_rectangle))
      {
        roi_rectangle[1] <- img_Xmax - roi_rectangle_pre[2] #Left maps Right
        roi_rectangle[2] <- img_Xmax - roi_rectangle_pre[1] + 1 #Rigth maps Left
        roi_rectangle[3] <- roi_rectangle_pre[3] -1 #Bottom  maps Top
        roi_rectangle[4] <- roi_rectangle_pre[4] #Top maps Bottom
      }
      else
      {
        roi_rectangle[1] <- zoom_rectangle[2] + zoom_rectangle[1]  - roi_rectangle_pre[2] - 1 #Left maps Right
        roi_rectangle[2] <- zoom_rectangle[2] + zoom_rectangle[1]  - roi_rectangle_pre[1] #Rigth maps Left
        roi_rectangle[3] <- img_Ymax - zoom_rectangle[3] - zoom_rectangle[4]  + roi_rectangle_pre[3] #Bottom  maps Top
        roi_rectangle[4] <- img_Ymax - zoom_rectangle[3] - zoom_rectangle[4]  +roi_rectangle_pre[4] + 1 #Top maps Bottom
      }
    }

    rasterRGB <- terra::flip(terra::flip(rasterRGB, direction = "vertical"), direction = "horizontal")#180degree rotation
  }

  #Add a border pixels and plot image
  terra::ext(rasterRGB) <- border + c( terra::ext(rasterRGB)$xmin,  terra::ext(rasterRGB)$xmax,  terra::ext(rasterRGB)$ymin,  terra::ext(rasterRGB)$ymax)
  rasterRGB <- terra::extend(rasterRGB, terra::ext( terra::ext(rasterRGB)$xmin - border,  terra::ext(rasterRGB)$xmax + border,  terra::ext(rasterRGB)$ymin - border,  terra::ext(rasterRGB)$ymax + border), fill = 0)
  terra::plotRGB(rasterRGB, axes = T, smooth = F )

  if(display_axes)
  {
    #Add calibrated axes
    xAxis<- seq(0,  terra::ext(rasterRGB)$xmax, by = ( terra::ext(rasterRGB)$xmax/10))
    xLabels <- sprintf( "%0.1f", xAxis * cal_um2pixels)
    yAxis<- seq(0,  terra::ext(rasterRGB)$ymax, by = ( terra::ext(rasterRGB)$ymax/10))
    yLabels <- sprintf( "%0.1f", yAxis * cal_um2pixels)
    par(xaxt = "l", yaxt = "l")
    axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = yLabels, las = 1) #Y left axes
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos =  terra::ext(rasterRGB)$xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
    axis(side=1, tck = -0.015, cex.axis = 0.7, pos = 0, at = xAxis, labels = xLabels ) #X below axes
    axis(side=3, tck = -0.015, cex.axis = 0.7, pos =  terra::ext(rasterRGB)$ymax, at = xAxis, labels = xLabels ) #X avobe axes
  }
  else
  {
    #Add calibrated scale bar
    Lp <- 0.001
    Hp <- 0.015
    WpTarget <- 0.15

    #Cal the most elegant nearest value
    legend_possible_values <- as.vector(sapply(10^(1:4), function(x){ x*(1:9) }))
    cal_length <- legend_possible_values[which.min(abs(legend_possible_values - WpTarget*terra::ext(rasterRGB)$xmax*cal_um2pixels))]
    Wp <- cal_length/(cal_um2pixels* terra::ext(rasterRGB)$xmax)
    yB <- Lp* terra::ext(rasterRGB)$ymax
    yT <- (Lp + Hp)* terra::ext(rasterRGB)$ymax

    if( rotation == 180  && display_syscoords)
    {
      #Avoid overlapping scale and axes
      xL <-  terra::ext(rasterRGB)$xmin + (Lp + Wp)*( terra::ext(rasterRGB)$xmax -  terra::ext(rasterRGB)$xmin)
      xR <-  terra::ext(rasterRGB)$xmin + Lp * ( terra::ext(rasterRGB)$xmax -  terra::ext(rasterRGB)$xmin)
    }
    else
    {
      xL <-  terra::ext(rasterRGB)$xmax - (Lp + Wp)*( terra::ext(rasterRGB)$xmax -  terra::ext(rasterRGB)$xmin)
      xR <-  terra::ext(rasterRGB)$xmax - Lp * ( terra::ext(rasterRGB)$xmax -  terra::ext(rasterRGB)$xmin)
    }
    lines( c( xL, xL, xR, xR), c( yB, yT, yT, yB ), col = "white", lwd = 2 )
    text( x = (xL + 0.5*(xR - xL)), y = 1.4*yT, labels = sprintf("%0.0f um", cal_length), col = "white", cex = 0.8, adj = c(0.5,0))

    #Add coors system arrows
    if(display_syscoords)
    {
      raster_size <- c( terra::ext(rasterRGB)$xmax- terra::ext(rasterRGB)$xmin,  terra::ext(rasterRGB)$ymax- terra::ext(rasterRGB)$ymin)
      arrow_length <- 0.25*min(raster_size)
      if( rotation == 0  )
      {
        P_0 <- c( terra::ext(rasterRGB)$xmin + 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymax - 0.01*(raster_size[2]))
        P_X <- c( terra::ext(rasterRGB)$xmin + arrow_length,  terra::ext(rasterRGB)$ymax - 0.01*(raster_size[2]))
        P_Y <- c( terra::ext(rasterRGB)$xmin + 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymax - arrow_length)
        Txt_Adj <- c(0, 1)
      }
      if( rotation == 90 )
      {
        P_0 <- c( terra::ext(rasterRGB)$xmin + 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymin + 0.01*(raster_size[2]))
        P_X <- c( terra::ext(rasterRGB)$xmin + 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymin + arrow_length)
        P_Y <- c( terra::ext(rasterRGB)$xmin + arrow_length,  terra::ext(rasterRGB)$ymin + 0.01*(raster_size[2]))
        Txt_Adj <- c(0, 0)
      }
      if( rotation == 270 )
      {
        P_0 <- c( terra::ext(rasterRGB)$xmax - 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymax - 0.01*(raster_size[2]))
        P_X <- c( terra::ext(rasterRGB)$xmax - 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymax - arrow_length)
        P_Y <- c( terra::ext(rasterRGB)$xmax - arrow_length,  terra::ext(rasterRGB)$ymax - 0.01*(raster_size[2]))
        Txt_Adj <- c(1, 1)
      }
      if( rotation == 180 )
      {
        P_0 <- c( terra::ext(rasterRGB)$xmax - 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymin + 0.01*(raster_size[2]))
        P_X <- c( terra::ext(rasterRGB)$xmax - arrow_length,  terra::ext(rasterRGB)$ymin + 0.01*(raster_size[2]))
        P_Y <- c( terra::ext(rasterRGB)$xmax - 0.01*(raster_size[1]),  terra::ext(rasterRGB)$ymin + arrow_length)
        Txt_Adj <- c(1, 0)
      }

      arrows(x0 = P_0[1], y0 = P_0[2], x1 = P_X[1], y1 = P_X[2], code = 2, lwd = 2, col = "white", length = 0.1)
      text ( x = P_X[1], y = P_X[2], labels = " X", col = "white", adj = Txt_Adj, cex = 0.8)
      arrows(x0 = P_0[1], y0 = P_0[2], x1 = P_Y[1], y1 = P_Y[2], code = 2, lwd = 2, col = "white", length = 0.1)
      text ( x = P_Y[1], y = P_Y[2], labels = " Y", col = "white", adj = Txt_Adj, cex = 0.8)
    }
  }

  if(!is.null(roi_rectangle))
  {
    #roi_rectangle  is c(left, rigth, top, bottom)
    rect( xleft = roi_rectangle[1] + border, xright = roi_rectangle[2]  + border, ytop = roi_rectangle[4] + border, ybottom = roi_rectangle[3] + border, border = "red", lwd = 2 )
  }
}

.plotIntensityScale<-function(img, color = NULL, light=3, intensity_limit = NULL, fixGtkMargin = -1) #fixGtkMargin is used as param to set the line of mtext in gtk
{
  max_int<-max(terra::values(img$raster))
  min_int<-min(terra::values(img$raster))

  #Create the raster
  ncols_scale <- 10
  scale_raster <- terra::rast( nrow = 255, ncol = ncols_scale, xmin= 0, xmax= ncols_scale, ymin= 0, ymax= 255)
  terra::values(scale_raster) <- as.vector(matrix(seq(from=max_int, to=min_int, length.out = 255), nrow = ncols_scale, ncol = 255, byrow = T))

  #Check RGB channels
  if( is.null(color))
  {
    #Remap Color 2 rainbow Space  (24bits color space, 8 bits per channel 255 steps)
    if(is.null(intensity_limit))
    {
      RGB_raster<-.ReMappingIntensity2HSV(scale_raster, value_multiplier = light )
    }
    else
    {
      RGB_raster<-.ReMappingIntensity2HSV(scale_raster, value_multiplier = light, maxN = intensity_limit)
    }
  }
  else
  {
    img_zero<-.InitRGBEmptyRaster( terra::ncol(scale_raster), terra::nrow(scale_raster) )
    terra::values(scale_raster) <- .NormalizeFrom0to1(terra::values(scale_raster))
    img_255 <-  .NormalizeTo255(scale_raster, light)
    if(color == "R")
    {
      RGB_raster <- img_255
      terra::add(RGB_raster) <- img_zero$raster
      terra::add(RGB_raster) <- img_zero$raster
    }
    if( color == "G" )
    {
      RGB_raster <- img_zero$raster
      terra::add(RGB_raster) <- img_255
      terra::add(RGB_raster) <- img_zero$raster
    }
    if( color == "B")
    {
      RGB_raster <- img_zero$raster
      terra::add(RGB_raster) <- img_zero$raster
      terra::add(RGB_raster) <- img_255
    }
  }


  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.7, mar = c(3,2,3,2),  mgp = c(3, 0.5, 0.5), oma = c(2,2,2,2))

  terra::plotRGB(RGB_raster, axes = T, smooth = T  )

  #Add axes
  yAxis<- seq(0, terra::ext(RGB_raster)$ymax, length.out = 11)
  yLabels <- sprintf( "%0.1e", seq(min_int, max_int, length.out = 11))
  par(xaxt = "l", yaxt = "l")
  axis(side=2, tck = -0.015, cex.axis = 0.6, pos = 0, at = yAxis, labels = F, las = 1) #Y left axes
  if( max_int == 0 )
  {
    axis(side=4, tck = -0.015, cex.axis = 0.6, pos = terra::ext(RGB_raster)$xmax, at = yAxis, labels = F) #Y right axes
  }
  else
  {
    axis(side=4, tck = -0.015, cex.axis = 0.6, pos = terra::ext(RGB_raster)$xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
  }

  axis(side = 1, tck = -0.015, cex.axis = 0.6, labels = F, pos = 0, at = c(0,ncols_scale))
  axis(side = 3, tck = -0.015, cex.axis = 0.6, labels = F, pos = terra::ext(RGB_raster)$ymax, at = c(0,ncols_scale))

  #Add the main title
  if(max_int == 0)
  {
    mtext("Channel Disabled", side = 2, line = fixGtkMargin, cex = 0.8, adj = 0.5  )
  }
  else
  {
    if(is.character(img$mass))
    {
      mtext(img$mass, side = 2, line = fixGtkMargin, cex = 0.8, adj = 0.5  )
    }
    else
    {
      mtext(sprintf("m/z: %0.4f+/-%0.3f Da", img$mass, img$tolerance), side = 2, line = fixGtkMargin, cex = 0.8, adj = 0.5)
    }
  }
}

#Combinar funcions anteior per fer aixo facil
#La idea es que si mass.peak i tolerance es passen com a vectors de 1 a 3 valors es fan imatges RGB
#' Plot a mass image of a up to 3 selected ions.
#'
#' @param img an rMSI data object.
#' @param mass.peak a vector of up to 3 elements containing the mass of ions to plot.
#' @param tolerance a vector of up to 3 elements containing the mass range to plot for each ion.
#' @param XResLevel the interpolation factor (default is 3).
#' @param NormalizationCoefs optional parameter with a matrix containing the scaling factors for each pixel.
#' @param rotation the rotation of image expressed in deg. Valid values are: 0, 90, 180 and 270.
#' @param show_axes if true axis will be plotted. Otherwise a um scale and xy arrows  are drawn.
#' @param scale_to_global_intensity scale the image intensity to fit into the global intensity (only for RGB).
#' @param vlight the lighting of the plotted image.
#' @param crop_area the region of image to show in a vector formated as c(left, right, bottom, top) if null the whole image is plot.
#' @param intensity_limit limit the all pixels intensities to a given value. If null no limiting is used.
#'
#' Plots a mass image. If only one mass ion are used a rainbow color code image is generated. If more ions are used each mass will
#' be encoded in an RGB color channel.
#'
#' @export
#'
plotMassImageByPeak<-function(img, mass.peak, tolerance=0.25, XResLevel = 3, NormalizationCoefs = NULL, rotation = 0, show_axes = F, scale_to_global_intensity = F, vlight = 3, crop_area = NULL, intensity_limit = NULL)
{
  numberOfChannels <- 1

  if(length(mass.peak) == 1 )
  {
    #Single Ion image
    im_sgn<-.buildImageByPeak(img, mass.peak, tolerance, NormalizationCoefs)

    #Apply limit intensity to raster object
    if(!is.null(intensity_limit) && !is.na(intensity_limit[1]) )
    {
      terra::values(im_sgn$raster)[ terra::values(im_sgn$raster) > intensity_limit ] <- intensity_limit[1]
    }

    if(scale_to_global_intensity)
    {
      raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel, global_intensity_scaling_factor = max(img$mean), light = vlight)
    }
    else
    {
      raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel, global_intensity_scaling_factor = intensity_limit[1], light = vlight)
    }

  }
  else
  {
    #Multiple Ions image
    if(length(tolerance) == 1)
    {
      tolerance[2] <- tolerance[1]
    }

    im_R<-.buildImageByPeak(img, mass.peak[1], tolerance[1], NormalizationCoefs)
    im_G<-.buildImageByPeak(img, mass.peak[2], tolerance[2], NormalizationCoefs)

    #Apply limit intensity to raster objects R and G
    if(!is.null(intensity_limit))
    {
      if(!is.na(intensity_limit[1]))
      {
        terra::values(im_R$raster)[ terra::values(im_R$raster) > intensity_limit ] <- intensity_limit[1]
      }
      if(!is.na(intensity_limit[2]))
      {
        terra::values(im_G$raster)[ terra::values(im_G$raster) > intensity_limit ] <- intensity_limit[2]
      }
    }

    numberOfChannels <- 2

    if( length(mass.peak) == 2 )
    {
      #Use Red and Green only
      im_B<-.InitRGBEmptyRaster( img$size["x"], img$size["y"] )
    }
    else
    {
      #Use RGB
      if(length(tolerance) == 2)
      {
        tolerance[3] <- tolerance[1]
      }
      im_B<-.buildImageByPeak(img, mass.peak[3], tolerance[3], NormalizationCoefs)

      #Apply limit intensity to raster object B
      if(!is.null(intensity_limit))
      {
        if(!is.na(intensity_limit[3]))
        {
          terra::values(im_B$raster)[ terra::values(im_B$raster) > intensity_limit ] <- intensity_limit[3]
        }
      }

      numberOfChannels <- 3
    }

    raster_RGB <-.BuildRGBImage( im_R, im_G, im_B,  XResLevel = XResLevel)
  }

  oldPar <- par(no.readonly = T)
  layout( matrix( (numberOfChannels+1):1, ncol = (1+numberOfChannels), nrow = 1, byrow = TRUE ), widths = c(7, rep(1, numberOfChannels)) )

  if(numberOfChannels == 1 )
  {
    .plotIntensityScale(im_sgn, light =  vlight, intensity_limit = intensity_limit[1])
  }
  else
  {
    if( numberOfChannels == 3 )
    {
      .plotIntensityScale(im_B, "B" )
    }
    .plotIntensityScale(im_G, "G" )
    .plotIntensityScale(im_R, "R" )
  }

  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotation, display_axes = show_axes, roi_rectangle = crop_area )
  par(oldPar)
}

.FillSimpleRaster <- function(img, values, text)
{
  #Fill data matrix
  layer<-matrix(0, nrow=img$size["x"], ncol=img$size["y"]) #Now I'm using a zero instead of NA to display a completely black background
  for( i in 1:nrow(img$pos))
  {
    layer[img$pos[ i , "x" ], img$pos[ i , "y" ]] <- values[i]
  }

  #Create the raster objects
  Ras <- terra::rast( nrow = ncol(layer), ncol = nrow(layer), xmin= 0, xmax= nrow(layer), ymin= 0, ymax= ncol(layer))
  terra::values(Ras) <- as.vector(layer)

  return(list(raster = Ras, mass = text, tolerance = "", cal_resolution = img$pixel_size_um))
}


#' plotRGBDataOnImg: Function to plot generic data on a RGB raster using rMSI pixel locations.
#'
#' @param img an rMSI object.
#' @param Rvalues a vector of values to be ploted as red layer sorted according IDs of rMSI object.
#' @param Gvalues a vector of values to be ploted as green layer sorted according IDs of rMSI object.
#' @param Bvalues a vector of values to be ploted as blue layer sorted according IDs of rMSI object.
#' @param RText Text to plot as Red scale label.
#' @param GText Text to plot as Green scale label.
#' @param BText Text to plot as Blue scale label.
#' @param XResLevel the interpolation to use in plot.
#' @param light the light to aplied to raster.
#' @param rotation image rotation in degree.
#'
#' @export
#'
plotRGBDataOnImg <- function(img, Rvalues, Gvalues, Bvalues, RText, GText, BText, XResLevel = 3 , light = 3, rotation = 0)
{
  Rraster<-.FillSimpleRaster(img, Rvalues, RText)
  Graster<-.FillSimpleRaster(img, Gvalues, GText)
  Braster<-.FillSimpleRaster(img, Bvalues, BText)

  raster_RGB <-.BuildRGBImage( Rraster, Graster, Braster,  XResLevel = XResLevel, light = light)

  oldPar <- par(no.readonly = T)
  layout( matrix( (4:1), ncol = 4, nrow = 1, byrow = TRUE ), widths = c(7, rep(1, 3)) )

  .plotIntensityScale(Braster, "B" )
  .plotIntensityScale(Graster, "G" )
  .plotIntensityScale(Rraster, "R" )
  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotation, display_axes = F)
  
  par(oldPar)
}

#' Plots a image using various MS image objects and the same intensity scale for every MS image.
#'
#' @param ... MS images in rMSI format as various arguments to be plotted together.
#' @param mass.peak m/z selected to plot MS image
#' @param tolerance mass tolerance to represent the MS image.
#' @param NormalizationName a string representing the normalization to use in the plot.
#' @param rotation a vector of rotation in degree to perform on plot for each image.
#'        If only one value is supplied, then all image will be rotated the same.
#' @param flipV vector of booleans indicating if images must be flipped vertically.
#' @param flipH vector of booleans indicating if images must be flipped horizontally.
#' @param light the lighting factor used in the plot.
#'
#' @export
plotVariousMassImagesByPeak <- function(..., mass.peak, tolerance = 0.25 , NormalizationName = NULL, rotation = 0, flipV = F, flipH = F, light = 5)
{
  MAIN_BORDER_SIZE <- 20
  SUB_BORDER_SIZE <- 15
  imgList <- list(...)

  #Check resultions
  pxRes <- unlist(lapply(imgList, function(x){ return(x$pixel_size_um) }))
  if(any(pxRes[1] != pxRes) )
  {
    cat("WARNING: At least one image have a different pixel resolution\n The used size scale bar may be inacurate\n")
  }

  #Compute rasters
  if(is.null(NormalizationName))
  {
    rasterImgs <- lapply(imgList, function(x){ .buildImageByPeak(x, mass.peak, tolerance)})
  }
  else
  {
    notAvailNorms <-unlist(lapply(imgList, function(x) {
      nullNorm <- is.null(x$normalizations[[NormalizationName]]);
      if(nullNorm)
      {
        cat(paste("Error: Image", x$name, "does not have normalization with name:", NormalizationName, "available\n"))
      }
      return(nullNorm);
    }))
    if(any(notAvailNorms))
    {
      stop("Error, not available normalization method for at least one image.\n")
    }
    rasterImgs <- lapply(imgList, function(x){ .buildImageByPeak(x, mass.peak, tolerance, x$normalizations[[NormalizationName]])})
  }

  multiImgWin <- gWidgets2::gwindow(title = "Multiple MS images view", visible = F)
  grpTop <- gWidgets2::ggroup(horizontal = F, container = multiImgWin)
  imaging_dev <- gWidgets2::ggraphics(spacing = 5 )
  gWidgets2::size( imaging_dev )<- c(800, 600)
  gWidgets2::add(obj = grpTop, child = imaging_dev,  fill = T, expand = T)
  grpCtl <- gWidgets2::ggroup(horizontal = T, container = grpTop)
  lblMz <- gWidgets2::glabel(text = paste("m/z", round(mass.peak, digits = 4), "+/-", round(tolerance, digits = 2)), container = grpCtl)
  lblNorm <- gWidgets2::glabel(text = paste("Normalization:", NormalizationName), container = grpCtl)
  gWidgets2::addSpring(grpCtl)
  btnSave <- gWidgets2::gbutton(text = "Save to png", container = grpCtl)
  gWidgets2::addHandlerClicked(btnSave, .saveGgraphicsPlot2PngFile, action = imaging_dev)
  gWidgets2::visible(multiImgWin) <- T
  gWidgets2::visible(imaging_dev) <- T
  plot.new()
  title(main = "Generating multi MS image, please wait...")


  #Method to apply rotations to a raster object
  rotateRaster <- function(rasterImg, rotation)
  {
    if(rotation == 0)
    {
      return(rasterImg)
    }
    if(rotation == 90)
    {
      return(terra::flip(terra::t(rasterImg), direction = "vertical"))
    }
    if(rotation == 180)
    {
      return(terra::flip(terra::flip(rasterImg, direction = "vertical"), direction = "horizontal"))
    }
    if(rotation == 270)
    {
      return(terra::flip(terra::t(rasterImg), direction = "horizontal"))
    }
    stop("ERROR, invalid rotation parameter. Valid parameters are: 0, 90, 180 and 270\n")
  }

  #Compute max intenssity for all images
  idMaxInt <- 1
  maxIntensity <- 0
  rotation <- c(rotation, rep(rotation[1], length(rasterImgs) - length(rotation)))
  flipH <- c(flipH, rep(F, length(rasterImgs) - length(flipH)))
  flipV <- c(flipV, rep(F, length(rasterImgs) - length(flipV)))
  for( i in 1:length(rasterImgs))
  {
    currMax <- max(terra::values(rasterImgs[[i]]$raster))
    if(currMax > maxIntensity)
    {
      maxIntensity <- currMax
      idMaxInt <- i
    }
    if(flipH[i])
    {
      rasterImgs[[i]]$raster <- terra::flip(rasterImgs[[i]]$raster, direction = "horizontal")
    }
    if(flipV[i])
    {
      rasterImgs[[i]]$raster <- terra::flip(rasterImgs[[i]]$raster, direction = "vertical")
    }
    rasterImgs[[i]]$raster <- rotateRaster(rasterImgs[[i]]$raster, rotation[i])
  }

  #Layout various rasters horizontally
  xrangeH <- sum(unlist(lapply(rasterImgs, function(x){ terra::ext(x$raster)$xmax })))
  yrangeH <- max(unlist(lapply(rasterImgs, function(x){ terra::ext(x$raster)$ymax })))

  #Layout various rasters vertically
  xrangeV <- max(unlist(lapply(rasterImgs, function(x){ terra::ext(x$raster)$xmax })))
  yrangeV <- sum(unlist(lapply(rasterImgs, function(x){ terra::ext(x$raster)$ymax })))

  #Choose betweeb horizontal vs vertical layout
  aspectRatio <- c( xrangeH/yrangeH, xrangeV/yrangeV )
  names(aspectRatio) <- c("H", "V")
  selLayout <- names(which.min(abs(1 - aspectRatio)))
  if( selLayout == "H")
  {
    xrange <- xrangeH + (SUB_BORDER_SIZE * length(rasterImgs))
    yrange <- yrangeH
    vMode <- 0
    hMode <- 1
    imgLabelRotation <- 90
  }
  else
  {
    xrange <- xrangeV
    yrange <- yrangeV + (SUB_BORDER_SIZE * length(rasterImgs))
    vMode <- 1
    hMode <- 0
    imgLabelRotation <- 0
  }

  #Prepara main raster
  globCanvas <- terra::ext(0, xrange, 0, yrange)
  xoffset <- 0
  yoffset <- yrange + SUB_BORDER_SIZE
  rasterList <- lapply(rasterImgs, function(x){ return(x$raster) })
  txtXPos <- c()
  txtYPos <- c()
  for( i in 1:length(rasterList))
  {
    yoffset <- vMode * (yoffset - terra::ext(rasterList[[i]])$ymax - SUB_BORDER_SIZE) #This is only true for vertical layout
    txtXPos[i] <- xoffset
    txtYPos[i] <- yoffset
    terra::ext(rasterList[[i]]) <- c( xoffset , terra::ext(rasterList[[i]])$xmax + xoffset, yoffset, terra::ext(rasterList[[i]])$ymax + yoffset)
    xoffset <- hMode * (terra::ext(rasterList[[i]])$xmax + SUB_BORDER_SIZE)#This is only true for horitzontal layout
    rasterList[[i]] <- terra::ext(rasterList[[i]], globCanvas, value = 0)
  }

  #Combain all rasters in a mosaic
  rasterList$fun <- max
  rasterList$na.rm <- TRUE
  RGBMosaic <- do.call(terra::mosaic, rasterList)
  rm(rasterList)
  RGBMosaic <- list(raster = RGBMosaic) #Workaround to be able to use .BuildSingleIonRGBImage function
  RGBMosaic <- .BuildSingleIonRGBImage(RGBMosaic, XResLevel = 3, global_intensity_scaling_factor = maxIntensity, light = light)

  #Plot data
  layout( matrix( 1:2, ncol = 2, nrow = 1, byrow = TRUE ), widths = c(7 , 1))
  .plotMassImageRGB(RGBMosaic, cal_um2pixels = imgList[[i]]$pixel_size_um, display_axes = F, display_syscoords = F, border = MAIN_BORDER_SIZE)
  txtYPos <- txtYPos + MAIN_BORDER_SIZE
  txtXPos <- txtXPos + MAIN_BORDER_SIZE
  rect( xleft = txtXPos,
        xright = txtXPos +  unlist(lapply(rasterImgs, function(x){ return( terra::ext(x$raster)$xmax - terra::ext(x$raster)$xmin )  }) ),
        ybottom = txtYPos,
        ytop = txtYPos +  unlist(lapply(rasterImgs, function(x){ return( terra::ext(x$raster)$ymax - terra::ext(x$raster)$ymin )  }) ),
        border = "white"
  )

  text(x = txtXPos - hMode, y = vMode + txtYPos + vMode*unlist(lapply(rasterImgs, function(x){ return( terra::ext(x$raster)$ymax - terra::ext(x$raster)$ymin ) }) ),
       labels = unlist(lapply(imgList, function(x){return(x$name)})), adj = c(0, 0), cex = 0.8, col ="white", srt = imgLabelRotation)
  .plotIntensityScale(rasterImgs[[idMaxInt]], intensity_limit = maxIntensity)
}

#' Save a gWidgets2::ggraphic widget to a png file.
#' A file dialog is presented to allow the user selecting where png file will be placed.
#'
#' @param evt the evt$action must be a pointer to gWidgets2::ggraphic to save.
#' @param ...
#'
.saveGgraphicsPlot2PngFile <- function(evt, ...)
{
  png_width <- evt$action$get_size()["width"]
  png_height <- evt$action$get_size()["height"]
  fname <- gWidgets2::gfile(text = "Choose location to save plot as png file", type = "save", initial.filename = "~/", multi = F)
  if( is.null(fname))
  {
    return()
  }
  if( length(fname) == 0)
  {
    return()
  }
  if( tolower(unlist(strsplit(fname, split = "\\."))[length( unlist(strsplit(fname, split = "\\.")))]) != "png" )
  {
    fname <- paste(fname, ".png", sep = "")
  }
  gWidgets2::visible(evt$action) <- T
  dev.print(png, fname, width = png_width, height = png_height)
}

