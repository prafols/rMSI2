#NEW API functions definitions

#Idem k anterior xo sense:
# - selectedPixels: Ho solucionare amb raster::crop raster::zoom i limatant el valor maxim de Z en plot del raster per cada layer RGB
# - rotate: Ho solucionare directament rotant objecte de raster en plotar, no em cal re-generar imatge per plotar
# - Retorna: una llista on
#         img.data correspondra al raster de la imatge,
#         mass.peak
#         tolerance
# El resultat d'aquest metode es guardara en la classe MSImagePlotWidget en una llista de capes R,G,B
.buildImageByPeak<-function(img, mass.peak, tolerance=0.25, NormCoefs = NULL )
{
  #Grab only MassPeaks Lists
  if(typeof(img) != "list")
  {
    stop("Error in rMSI::.buildImageByPeak(), img is not a list")
  }

  img_width<-img$size["x"]
  img_height<-img$size["y"]

  #Normalization coeficients, 1 if null
  if(is.null(NormCoefs))
  {
    NormCoefs<-rep(1, nrow(img$pos))
  }

  l1 <-which.min( abs( img$mass - ( mass.peak - tolerance ) ) )
  l2 <-which.min( abs( img$mass - ( mass.peak + tolerance ) ) )

  l1 <- l1[1] #Store only the first element
  l2 <- l2[length(l2)] #store only the last element
  mass.peak <- round(mean(c(img$mass[l2] , img$mass[l1])), digits = 4)
  tolerance <- round(0.5*(img$mass[l2] - img$mass[l1]), digits = 4)

  z_indexes <- l1:l2
  if(length(z_indexes) == 1)
  {
    z_indexes <- c(z_indexes,z_indexes)
  }

  #Fill raster matrix
  zplots<-matrix(0, nrow=img_width, ncol=img_height) #Now I'm using a zero instead of NA to display a completely black background
  curr_id <- 0
  for( iCube in 1:length(img$data))
  {
    if( nrow(img$data[[iCube]]) > 1 )
    {
      valpixels<-apply(img$data[[iCube]][ , z_indexes], 1, mean)
      ##valpixels<-apply(img$data[[iCube]][ , z_indexes], 1, max)
    }
    else
    {
      valpixels<-mean(img$data[[iCube]][, z_indexes])
      #valpixels<-max(img$data[[iCube]][, z_indexes])
    }

    for( i in 1:nrow(img$data[[iCube]]))
    {
      curr_id <- curr_id + 1
      x<-img$pos[ curr_id , "x" ]
      y<- 1 + img_height - img$pos[ curr_id, "y" ] #Flip image verticaly to match ploting #TODO oju k si la nova rutina xuta be aixo canviara!
      zplots[x,y]<- valpixels[i] * NormCoefs[curr_id] #Draw pixel apply normalization
    }
  }

  #Create the raster
  my_raster <- raster::raster( nrow = nrow(zplots), ncol = ncol(zplots), xmn= 0, xmx= ncol(zplots), ymn= 0, ymx= nrow(zplots))
  raster::values(my_raster) <- as.vector(t(zplots))

  #Return zplots matrix and some metadata in a list
  list(raster = my_raster, mass = mass.peak, tolerance = tolerance, cal_resolution = img$pixel_size_um)
}


#Create an empty raster object
.InitRGBEmptyRaster<-function( width, height)
{
  new_raster<-raster::raster( nrow = width, ncol = height, xmn= 0, xmx= height, ymn= 0, ymx= width)
  raster::values(new_raster) <- rep(0, width*height)
  return( list( raster = new_raster, mass = NULL, tolerance = NULL, cal_resolution = NULL))
}

#Idem k anterior plotMassImageRGB
#Com a param img_RGB es una llista amb atributs R,G,B on cada objecte (layer) es un objecte retornat per builImageByPeak
# XResLevel es ara un integer ja que delego la interpolacio a raster
# All RGB layer must have the same size
.BuildRGBImage <- function( imgR, imgG, imgB, XResLevel = 3 )
{
  #Normalize to 255 (8 bits per color channel)
  NormalizeTo255 <-function( m )
  {
    maxN<- max(raster::values(m))
    if(maxN > 0)
    {
      raster::values(m) <- 255 * raster::values(m) / maxN
    }
    return(m)
  }
  imgR <- NormalizeTo255(imgR$raster)
  imgG <- NormalizeTo255(imgG$raster)
  imgB <- NormalizeTo255(imgB$raster)

  #Create an RGB image space
  RGB_raster <- raster::addLayer(imgR,imgG,imgB )
  interpolated_raster <- raster::raster( nrow= XResLevel*RGB_raster@nrows, ncol= XResLevel*RGB_raster@ncols, xmn= 0, xmx= RGB_raster@ncols, ymn= 0, ymx= RGB_raster@nrows)
  RGB_raster<-raster::resample(RGB_raster, interpolated_raster)
  raster::values(RGB_raster)[ raster::values(RGB_raster) < 0  ] <- 0 #Values below zero are dube interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

#Remap a intensity vector to HSV coloring in a rainbow function
.ReMappingIntensity2HSV<-function(single_channel_raster)
{
  #Normalize to 1
  maxN<- max(raster::values(single_channel_raster))
  if(maxN > 0)
  {
    raster::values(single_channel_raster) <- raster::values(single_channel_raster) / maxN
  }

  #Remapping hue space
  hue_top <- 0.7
  hue_bottom <- 0.85
  hMapped <- (1 + hue_top - hue_bottom) * (-1*raster::values(single_channel_raster)  + 1) + hue_bottom
  over_one <- which(hMapped > 1)
  hMapped[ over_one ] <- hMapped[over_one] -1

  #Remapping value space
  vMapped <- raster::values(single_channel_raster) *3
  vMapped[ vMapped > 1] <- 1

  #Creating HSV image
  rgbSpace <- col2rgb(hsv( h =  hMapped, s = rep(1, length(raster::values(single_channel_raster) )), v = vMapped))

  #Layering RGB raster
  R_raster<- single_channel_raster
  G_raster<- single_channel_raster
  B_raster<- single_channel_raster
  raster::values(R_raster) <-  rgbSpace[1, ] #Extract R
  raster::values(G_raster) <-  rgbSpace[2, ] #Extract G
  raster::values(B_raster) <-  rgbSpace[3, ] #Extract B

  #Create an RGB image space
  return (raster::addLayer(R_raster,G_raster,B_raster ))
}

#Build a RGB images raster using rainbow colors from only one raster layer
.BuildSingleIonRGBImage<-function( img,   XResLevel = 3 )
{
  #Create an RGB image space
  RGB_raster <- .ReMappingIntensity2HSV(img$raster)
  interpolated_raster <- raster::raster( nrow= XResLevel*RGB_raster@nrows, ncol= XResLevel*RGB_raster@ncols, xmn= 0, xmx= RGB_raster@ncols, ymn= 0, ymx= RGB_raster@nrows)
  RGB_raster<-raster::resample(RGB_raster, interpolated_raster)
  raster::values(RGB_raster)[ raster::values(RGB_raster) < 0  ] <- 0 #Values below zero are dube interpolation artifacts, clip it to zero.

  return(RGB_raster)
}

.plotMassImageRGB <- function(rasterRGB, cal_um2pixels = 1,  rotation=0, flipV=F, flipH=F, display_axes=T)
{
  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.6, mar = c(1,1,1,1), mgp = c(2, 0.5, 0.5))



  ###TODO: amb el param zlim = c(min, max) pot fer que tots els valors per sobre de max agafin color de max i per sota min de min
  ###      seria un bon metode per implementar un limitador a nivell de raster (molt eficient)

  #Apply rotation
  if( rotation == 90 )
  {
    rasterRGB <- raster::flip(raster::t(rasterRGB), direction = "y") #90º rotation
  }
  if( rotation == 270 )
  {
    rasterRGB <- raster::flip(raster::t(rasterRGB), direction = "x") #270º rotation
  }
  if( rotation == 180 )
  {
    rasterRGB <- raster::flip(raster::flip(rasterRGB, direction = "y"), direction = "x")#180º rotation
  }

  raster::plotRGB(rasterRGB, axes = T, asp = 1, interpolate = T )

  if(display_axes)
  {
    #Add calibrated axes
    xAxis<- seq(0, rasterRGB@extent@xmax, by = (rasterRGB@extent@xmax/10))
    xLabels <- sprintf( "%0.1f", xAxis * cal_um2pixels)
    yAxis<- seq(0, rasterRGB@extent@ymax, by = (rasterRGB@extent@ymax/10))
    yLabels <- sprintf( "%0.1f", yAxis * cal_um2pixels)
    par(xaxt = "l", yaxt = "l")
    axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = yLabels, las = 1) #Y left axes
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = rasterRGB@extent@xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
    axis(side=1, tck = -0.015, cex.axis = 0.7, pos = 0, at = xAxis, labels = xLabels ) #X below axes
    axis(side=3, tck = -0.015, cex.axis = 0.7, pos = rasterRGB@extent@ymax, at = xAxis, labels = xLabels ) #X avobe axes
  }
  else
  {
    #Add calibrated scale bar
    Lp <- 0.05
    Hp <- 0.015
    WpTarget <- 0.15

    #Cal the most elegant nearest value
    legend_possible_values <- as.vector(sapply(10^(1:4), function(x){ x*(1:9) }))
    cal_length <- legend_possible_values[which.min(abs(legend_possible_values - WpTarget*rasterRGB@extent@xmax*cal_um2pixels))]
    Wp <- cal_length/(cal_um2pixels*rasterRGB@extent@xmax)

    xL <- Lp*rasterRGB@extent@xmax
    xR <- (Lp + Wp)*rasterRGB@extent@xmax
    yB <- Lp*rasterRGB@extent@ymax
    yT <- (Lp + Hp)*rasterRGB@extent@ymax
    lines( c( xL, xL, xR, xR  ), c( yB, yT, yT, yB ), col = "white", lwd = 2 )
    text( x = (( Lp + 0.5*Wp )*rasterRGB@extent@xmax), y = 1.3*yT, labels = sprintf("%0.0f um", cal_length), col = "white", cex = 0.8, adj = c(0.5,0))
  }
}

.plotIntensityScale<-function(img, color = NULL)
{
  max_int<-max(raster::values(img$raster))

  #Create the raster
  ncols_scale <- 10
  scale_raster <- raster::raster( nrow = 255, ncol = ncols_scale, xmn= 0, xmx= ncols_scale, ymn= 0, ymx= 255)
  raster::values(scale_raster) <- as.vector(matrix(seq(from=max_int, to=0, length.out = 255), nrow = ncols_scale, ncol = 255, byrow = T))

  #Normalize to 255 (8 bits per color channel)
  NormalizeTo255 <-function( m )
  {
    maxN<- max(raster::values(m))
    if(maxN > 0)
    {
      raster::values(m) <- 255 * raster::values(m) / maxN
    }
    return(m)
  }

  #Check RGB channels
  if( is.null(color))
  {
    #Remap Color 2 rainbow Space  (24bits color space, 8 bits per channel 255 steps)
    RGB_raster<-.ReMappingIntensity2HSV(scale_raster)
  }
  else
  {
    img_zero<-.InitRGBEmptyRaster( scale_raster@nrows, scale_raster@ncols )
    img_255 <-  NormalizeTo255(scale_raster)
    if(color == "R")
    {
      RGB_raster <- raster::addLayer( img_255, img_zero$raster, img_zero$raster )
    }
    if( color == "G" )
    {
      RGB_raster <- raster::addLayer( img_zero$raster, img_255, img_zero$raster )
    }
    if( color == "B")
    {
      RGB_raster <- raster::addLayer( img_zero$raster, img_zero$raster, img_255 )
    }
  }


  #Setting my tricky par values...
  par( bg = "black", fg =  "white", col.lab="white", xaxt="n", yaxt="n", col.axis = "white", col.main = "white", col.sub = "white",
       cex.axis = 0.6, mar = c(1,0,1,2),  mgp = c(3, 0.5, 0.5))

  raster::plotRGB(RGB_raster, axes = T, asp = 1, interpolate = T  )

  #Add axes
  yAxis<- seq(0, RGB_raster@extent@ymax, length.out = 10)
  yLabels <- sprintf( "%0.1e", seq(0, max_int, length.out = 10))
  par(xaxt = "l", yaxt = "l")
  axis(side=2, tck = -0.015, cex.axis = 0.7, pos = 0, at = yAxis, labels = F, las = 1) #Y left axes
  if( max_int == 0 )
  {
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = RGB_raster@extent@xmax, at = yAxis, labels = F) #Y right axes
  }
  else
  {
    axis(side=4, tck = -0.015, cex.axis = 0.7, pos = RGB_raster@extent@xmax, at = yAxis, labels = yLabels, las = 1) #Y right axes
  }

  axis(side = 1, tck = -0.015, cex.axis = 0.7, labels = F, pos = 0, at = c(0,ncols_scale))
  axis(side = 3, tck = -0.015, cex.axis = 0.7, labels = F, pos = RGB_raster@extent@ymax, at = c(0,ncols_scale))

  #Add the main title
  if(max_int == 0)
  {
    mtext("Channel Disabled", side = 2, line = -1, cex = 0.8, adj = 0.5  )
  }
  else
  {
    mtext(sprintf("m/z: %0.3f+/-%0.2f Da", img$mass, img$tolerance), side = 2, line = -1, cex = 0.8, adj = 0.5  )
  }
}

#Combinar funcions anteior per fer aixo facil
#La idea es que si mass.peak i tolerance es passen com a vectors de 1 a 3 valors es fan imatges RGB
plotMassImageByPeak<-function(img, mass.peak, tolerance=0.25, XResLevel = 3, NormalizationCoefs = NULL, rotation = 0)
{
  numberOfChannels <- 1

  if(length(mass.peak) == 1 )
  {
    #Single Ion image
    im_sgn<-.buildImageByPeak(img, mass.peak, tolerance, NormalizationCoefs)
    raster_RGB<-.BuildSingleIonRGBImage(im_sgn, XResLevel = XResLevel)
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
      numberOfChannels <- 3
    }

    raster_RGB <-.BuildRGBImage( im_R, im_G, im_B,  XResLevel = XResLevel)
  }


  layout( matrix( (numberOfChannels+1):1, ncol = (1+numberOfChannels), nrow = 1, byrow = TRUE ), widths = c(7, rep(1, numberOfChannels)) )


  if(numberOfChannels == 1 )
  {
    .plotIntensityScale(im_sgn )
  }
  else
  {
    .plotIntensityScale(im_R, "R" )
    .plotIntensityScale(im_G, "G" )
    if( numberOfChannels == 3 )
    {
      .plotIntensityScale(im_B, "B" )
    }
  }

  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotation, display_axes = F)
}