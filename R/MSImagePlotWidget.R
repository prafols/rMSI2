#########################################################################
#     rMSI2 - R package for MSI data handling and visualization
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

###A GUI to display MS images from some ions/features

.MSImagePlotWidget <- function( in_img, parent_widget=gwindow ( "Default MSImagePlotWidget" , visible = FALSE ), AddSpectra_function = NULL, GetSpectraInfo_function = NULL, ClearSpectraPlot_function = NULL, meanSpectrumColor = "red", widget_name = "")
{
  options(guiToolkit="tcltk") 
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  img <- in_img
  rm(in_img)
  parent <- parent_widget
  rm(parent_widget)
  mz_tolerance <- c()
  mz_selected <- c()
  GUI_RASTER_BORDER <- 10 #Pixels of border around ms image to allow better placement of axes and scales
  imaging_dev <- 0
  Spin_Smooth <- 0
  Spin_Xmin <- 0
  Spin_Xmax <- 0
  Spin_Ymin <- 0
  Spin_Ymax <- 0
  Rotation <- 0
  flipV <- F
  flipH <- F
  Tbl_spotList <- 0
  Scale_light <- 0
  plotting_raster <- .BuildSingleIonRGBImage( .InitRGBEmptyRaster( img$size["x"], img$size["y"] ),   XResLevel = 1, light =  1 ) #Current plotted MS image object
  ROI <- NULL #Current used ROI on image, NULL means no ROI
  ZOOM_win <- NULL #Current zoom windows
  IntLimit_ROI <- NULL #Intensity limiting ROI
  IntLimits <- NULL #A vector of intensity limits for each channel
  SidePanel_position <- 0 #Storing side panel position of spectra list in order to resotre it when it is hide/show
  AddSpectra_ptr <- AddSpectra_function #Pointer to a AddSpectra method of a spectraWidget to be able of ploting directly
  rm(AddSpectra_function)
  GetSpectraInfo_ptr <- GetSpectraInfo_function #Pointer to a GetSpectra method of a spectraWidget to be able to get current ploted spectra directly
  rm(GetSpectraInfo_function)
  ClearSpectraPlot_ptr <- ClearSpectraPlot_function #Pointer to a RmSpectra method of a spectraWidget to be able to clear ploted spectra directly
  rm(ClearSpectraPlot_function)
  myName <- widget_name
  NormalizationCoefs <- rep(1, nrow(img$pos))
  previous_spectralist_width<-180 #width in pixels

  #Current image RGB layers
  Rlayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Glayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Blayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  
  #Current plotting RGB layers
  red_layer <- Rlayer_raster
  green_layer <- Glayer_raster
  blue_layer <- Blayer_raster
  
  #Keep track of the number of channels to plot (RGB vs. single ion mode)
  ch_count <- 0

  ImgBuildFun <- function(channel, mass, tol )
  {
    this$mz_tolerance[channel] <- tol
    this$mz_selected[channel] <- mass
    img_new<-.buildImageByPeak(this$img, mass.peak = mass, tolerance = tol, NormCoefs = this$NormalizationCoefs)

    #Apply intensity limitation directly to the raster object
    if( !is.null(this$IntLimit_ROI))
    {
      roi_c <- c(this$IntLimit_ROI[1] -1, this$IntLimit_ROI[2],  terra::ext(img_new$raster)$ymax - this$IntLimit_ROI[4],  terra::ext(img_new$raster)$ymax - this$IntLimit_ROI[3] + 1)
      if(roi_c[1] == roi_c[2])
      {
        roi_c[2] <- roi_c[2] + 1
      }
      if(roi_c[1] > roi_c[2])
      {
        raux <- roi_c[1]
        roi_c[1] <- roi_c[2]
        roi_c[2] <- raux
      }
      
      if(roi_c[3] == roi_c[4])
      {
        roi_c[4] <- roi_c[4] + 1
      }
      if(roi_c[3] > roi_c[4])
      {
        raux <- roi_c[3]
        roi_c[3] <- roi_c[4]
        roi_c[4] <- raux
      }
      
      crop_ext <- terra::ext(roi_c)
       
      if(terra::ext(crop_ext)$xmin == terra::ext(crop_ext)$xmax && terra::ext(crop_ext)$ymin == terra::ext(crop_ext)$ymax)
      {
         #Single pixel selected
         this$IntLimits[channel] <- as.numeric(terra::extract(img_new$raster,cbind(terra::ext(crop_ext)$xmin, terra::ext(crop_ext)$ymin) ))
      }
      else
      {
         this$IntLimits[channel] <- max(terra::values( terra::crop( img_new$raster, crop_ext)))
      }
      
      if(this$IntLimits[channel] <= 0.0)
      {
        this$IntLimits <- NULL
      }
      else
      {
        terra::values(img_new$raster)[ terra::values(img_new$raster) > this$IntLimits[channel] ] <- this$IntLimits[channel]
      }
    }
    else
    {
      this$IntLimits <- NULL
    }

    mz_str <- this$mz_selected[channel]
    if(mz_str < 1000)
    {
      mz_str <- paste( substr( sprintf("%.10f", mz_str), 1, 6), "Da" )
    }
    else
    {
      mz_str <- paste( substr( sprintf("%.10f", mz_str/1000), 1, 5), "kDa" )
    }

    if( channel == 1)
    {
      this$Rlayer_raster<-img_new
      tcltk::tkconfigure(this$Btn_RedEnable, text = mz_str)
    }
    else if (channel == 2)
    {
      this$Glayer_raster<-img_new
      tcltk::tkconfigure(this$Btn_GreenEnable, text = mz_str)
    }
    else if (channel == 3)
    {
      this$Blayer_raster<-img_new
      tcltk::tkconfigure(this$Btn_BlueEnable, text = mz_str)
    }

    #Return del buildImage
    return(list(selMz = img_new$mass, selTol = img_new$tolerance))
  }

  RedrawMSImage <-function()
  {
    if(this$ch_count >0)
    {
      .plotMassImageRGB (this$plotting_raster, cal_um2pixels = this$img$pixel_size_um,  rotation = this$Rotation, flipV = this$flipV, flipH = this$flipH,
                           display_axes = F, roi_rectangle =  this$ROI, zoom = this$ZOOM_win, border = this$GUI_RASTER_BORDER)
    }
    else
    {
      par(bg = "black", fg = "white")
      plot.new( )
      plot.window( xlim=c(-50,50), ylim=c(-5,5) )
      text(0,0,"Selected an m/z channel to display its image", cex = 1.5 )
    }
  }

  PlotMassImageRGB <- function()
  {
    this$ch_count <- 0
    if( getValue_coloredCheckBox(this$Btn_RedEnable))
    {
      this$red_layer <- this$Rlayer_raster
      unique_layer <- this$red_layer
      this$ch_count <- this$ch_count + 1
    }
    else
    {
      this$red_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( getValue_coloredCheckBox(this$Btn_GreenEnable))
    {
      this$green_layer <- this$Glayer_raster
      unique_layer <- this$green_layer
      this$ch_count <- this$ch_count + 1
    }
    else
    {
      this$green_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( getValue_coloredCheckBox(this$Btn_BlueEnable))
    {
      this$blue_layer <- this$Blayer_raster
      unique_layer <- this$blue_layer
      this$ch_count <- this$ch_count + 1
    }
    else
    {
      this$blue_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }

    inter_level<-switch(svalue(this$Combo_Xres), x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5)
    if(this$ch_count == 1)
    {
      this$plotting_raster<-.BuildSingleIonRGBImage( unique_layer,   XResLevel = inter_level, light =  svalue(this$Scale_light) )
    }
    else
    {
      this$plotting_raster<-.BuildRGBImage( imgR = this$red_layer, imgG = this$green_layer, imgB = this$blue_layer, XResLevel = inter_level, light =  svalue(this$Scale_light) )
    }

    redrawPlotWidget(this$imaging_dev)
    redrawPlotWidget(this$scaleRed_dev)
    redrawPlotWidget(this$scaleGreen_dev)
    redrawPlotWidget(this$scaleBlue_dev)
  }
  
  ReDrawRedScale <- function()
  {
    if(this$ch_count == 1)
    {
      .plotIntensityScale(this$red_layer, light = svalue(this$Scale_light), fixGtkMargin = -1) 
    }
    else
    {
      .plotIntensityScale(this$red_layer, "R", light = svalue(this$Scale_light), fixGtkMargin = -1 ) 
    }
  }
  
  ReDrawGreenScale <- function()
  {
    if(this$ch_count == 1)
    {
      .plotIntensityScale(this$green_layer, light = svalue(this$Scale_light), fixGtkMargin = -1)
    }
    else
    {
      .plotIntensityScale(this$green_layer, "G", light = svalue(this$Scale_light), fixGtkMargin = -1)
    }
  }
  
  ReDrawBlueScale <- function()
  {
    if(this$ch_count == 1)
    {
      .plotIntensityScale(this$blue_layer, light = svalue(this$Scale_light), fixGtkMargin = -1)
    }
    else
    {
      .plotIntensityScale(this$blue_layer, "B", light = svalue(this$Scale_light), fixGtkMargin = -1)
    }
  }

  BtnClearSpotList <- function( mass, tol, ... )
  {
    this$Tbl_spotList$set_items(data.frame(this$Tbl_spotList$get_items())[1,])
    tweaksSpectraTable(this$Tbl_spotList)
  }

  SpectraListSelChange <- function( ... )
  {
    selected <- gWidgets2::svalue(this$Tbl_spotList)
    df<-data.frame(this$Tbl_spotList$get_items())  #get data frame...
    color_list<-c()
    id_list<-c()

    for( i in selected)
    {
      selDf<-df[ df$ID == i,] #Get the correct data row
      color_list <- c(color_list, as.character(selDf$Colour))
      id_list <- c( id_list, i)
    }

    #Add spectra to plot
    plotSpectra( id_list, color_list)
  }

  #Plot spectra using parent plot function including normalization
  plotSpectra <- function( id, colors, mz_min = NULL, mz_max = NULL)
  {
    intensity_list <- list()
    color_list <- list()
    id_list <- list()
    norm_vect <- c()

    id0_pos <- which(id == 0)
    if(length(id0_pos > 0))
    {
      #There is a mean spectrum
      if(class(this$img$mean) == "MassSpectrum")
      {
        #Addap to old data mean spectrum using MALDIquant object
        intensity_list[[length(intensity_list) + 1]] <- this$img$mean@intensity
      }
      else
      {
        intensity_list[[length(intensity_list) + 1]] <- this$img$mean
      }
      color_list[[length(color_list) + 1 ]] <- colors[id0_pos]

      inf_index <- which(is.infinite(this$NormalizationCoefs))
      if(length(inf_index) > 0)
      {
        norm_vect <- c(norm_vect, mean(this$NormalizationCoefs[-inf_index]))
      }
      else
      {
        norm_vect <- c(norm_vect, mean(this$NormalizationCoefs))
      }
      id_list[[length(id_list) + 1]] <-  0
    }

    #Data
    id_fil_pos <- which(id > 0)
    if(length(id_fil_pos) > 0)
    {
      intensity_data <- loadImgChunkFromIds(this$img, id[id_fil_pos])
      for( i in 1:length(id_fil_pos))
      {
        intensity_list[[length(intensity_list) + 1]] <- intensity_data[i, ]
        color_list[[length(color_list) + 1 ]] <- colors[id_fil_pos[i]]
        norm_vect <- c(norm_vect, this$NormalizationCoefs[id[id_fil_pos[i]]])
        id_list[[length(id_list) + 1]] <-  id[id_fil_pos[i]]
      }
    }

    if(!is.null(this$AddSpectra_ptr))
    {
      this$AddSpectra_ptr(this$img$mass, intensity_list, color_list, id_list, this$myName, norm_vect, mz_min, mz_max)
    }
  }

  BtnExportSpotList <- function( ... )
  {
    indexes<-as.vector(data.frame(this$Tbl_spotList$get_items())$ID)
    indexes<-indexes[indexes > 0] #Remove zero!
    store_paths<-gWidgets2::gfile("Save current spots to txt files", type="selectdir", multi = F, initial.dir = path.expand("~/"))
    if(length(store_paths) == 0)
    {
      return ()
    }

    #Export indexes to ID.txt file
    if(length(indexes) > 0)
    {
      #Save the id file
      id_fname <- file.path( store_paths, "ID.txt")
      write(indexes, file = id_fname, ncolumns = 1)
    }

    #Continue exporting TXT spectra
    bExportData <- gWidgets2::gconfirm(paste("ID's have been exported to ID.txt file\nDo you also want to store spectra as TXT?"), title = "Export spectra as TXT", icon = "question")
    if(!bExportData)
    {
      return()
    }

    #display a BIG WaRNIng IF to much ID's are seected!
    if( length(indexes) > 25)
    {
      bExportData<-gWidgets2::gconfirm(paste("You are exporting a lot of data. (", length(indexes) ,"mass spectrums )\nThis may take a long time and expend a lot of memory.\nDo you want to store spectra as TXT?"), title = "Warning: Large data export!", icon = "warning")
    }

    mPbar<-.ProgressBarDialog("Exporting spectra to txt...")

    #Create a dir to store all data inside
    store_paths <- file.path(store_paths, paste("Export_", tools::file_path_sans_ext(this$img$name),"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), sep = "" ))
    dir.create( store_paths )

    #Save mean spectra
    if(class(this$img$mean) == "MassSpectrum")
    {
      #Handling old mean based on MALDIquant
      spc <- matrix(data= c(this$img$mean@mass, this$img$mean@intensity), ncol = 2, byrow = F)
    }
    else
    {
      spc <- matrix(data= c(this$img$mass, this$img$mean), ncol = 2, byrow = F)
    }
    write( x = t(spc), file = file.path( store_paths, "average.txt" ), ncolumns = 2  )


    if(length(indexes) > 0 && bExportData)
    {
      #Save each ID in the lit
      dataChunck <- loadImgChunkFromIds(this$img, indexes)

      for( i in 1:length(indexes))
      {
        spc <- matrix(data= c(this$img$mass, dataChunck[i, ]), ncol = 2, byrow = F)
        spc_fname<-  file.path( store_paths, paste("ID_", indexes[i] ,".txt", sep =""))
        write( x = t(spc), file = spc_fname, ncolumns = 2  )
        if( !mPbar$setValue( 100* i/length(indexes) ))
        {
          #Aborted by user
          rm(dataChunck)
          unlink( store_paths, recursive = T )
          gc()
          return()
        }
      }
    }

    mPbar$close()

    rm(dataChunck)
    gc()
    #Show a message of export completed
    smsg<- paste("Export complete!\nYour data is stored at:\n\t", store_paths , sep = "")
    if(bExportData)
    {
      if(length(indexes) > 0)
      {
        smsg <- paste(smsg,"\nAverage spectrum and", length(indexes), "spectrums exported.")
      }
      else
      {
        smsg <- paste(smsg,"\nOnly average spetrum exported.")
      }

    }
    else
    {
      smsg <- paste(smsg, "\nLarge data export disabled. Ony average spectrum and ID list exported.")
    }
    gWidgets2::gmessage(smsg, title = "Export complete!", icon = "info")

  }

  OnPixelSelection <- function( x, y )
  {
    X_left<-round(min(x)) - this$GUI_RASTER_BORDER
    X_right<-round(max(x)) - this$GUI_RASTER_BORDER
    Y_bottom<-round(min(y)) - this$GUI_RASTER_BORDER #Transform raster coords to image coords (only Y axis is affected)
    Y_top<-round(max(y)) - this$GUI_RASTER_BORDER #Transform raster coords to image coords (only Y axis is affected)

    #Apply rotation!
    if(this$Rotation == 0)
    {
      this$ROI <- c(X_left + 1, X_right, this$img$size["y"] - Y_top + 1, this$img$size["y"] - Y_bottom )
    }
    if(this$Rotation == 90)
    {
      if(is.null(this$ZOOM_win))
      {
        this$ROI <- c(Y_bottom + 1, Y_top, X_left + 1, X_right)
      }
      else
      {
        this$ROI <- c(Y_bottom + 1, Y_top, this$ZOOM_win[3] + this$ZOOM_win[4] - this$img$size["y"] + X_left, this$ZOOM_win[3] + this$ZOOM_win[4] - this$img$size["y"] + X_right - 1)
      }
    }
    if(this$Rotation == 180)
    {
      if(is.null(this$ZOOM_win))
      {
        this$ROI <- c( this$img$size["x"] - X_right + 1, this$img$size["x"] - X_left, Y_bottom + 1, Y_top)
      }
      else
      {
        this$ROI <- c( this$ZOOM_win[2] + this$ZOOM_win[1] - X_right,
                       this$ZOOM_win[2] + this$ZOOM_win[1] - X_left - 1,
                       this$ZOOM_win[3] + this$ZOOM_win[4] + Y_bottom - this$img$size["y"],
                       this$ZOOM_win[3] + this$ZOOM_win[4] + Y_top - this$img$size["y"] - 1)
      }
    }
    if(this$Rotation == 270)
    {
      if(is.null(this$ZOOM_win))
      {
        this$ROI <- c( this$img$size["x"] - Y_top + 1, this$img$size["x"] - Y_bottom, this$img$size["y"] - X_right + 1, this$img$size["y"] - X_left)
      }
      else
      {
        this$ROI <- c( this$ZOOM_win[2] + this$ZOOM_win[1] - Y_top, this$ZOOM_win[2] + this$ZOOM_win[1] - Y_bottom - 1, this$img$size["y"] - X_right + 1, this$img$size["y"] - X_left)
      }
    }

    #Apply flip
    if((this$flipV && this$Rotation == 0) || (this$flipV && this$Rotation == 180) || (this$flipH && this$Rotation == 90) || (this$flipH && this$Rotation == 270))
    {
      aux <- this$img$size["y"]  - this$ROI[4]
      this$ROI[4] <- this$img$size["y"]  - this$ROI[3]
      this$ROI[3] <- aux
    }
    if((this$flipH && this$Rotation == 0) || (this$flipH && this$Rotation == 180) || (this$flipV && this$Rotation == 90) || (this$flipV && this$Rotation == 270))
    {
      aux <- this$img$size["x"] - this$ROI[2]
      this$ROI[2] <- this$img$size["x"] - this$ROI[1]
      this$ROI[1] <- aux
    }

    #Set it to ROI spinbuttons
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()
    svalue(this$Spin_Xmin) <- this$ROI[1]
    svalue(this$Spin_Xmax) <- this$ROI[2]
    svalue(this$Spin_Ymin) <- this$ROI[3]
    svalue(this$Spin_Ymax) <- this$ROI[4]
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    redrawPlotWidget(this$imaging_dev)
    gWidgets2::enabled(this$Btn_RoiZoom) <- T
    gWidgets2::enabled(this$Frame_RoiCtl) <- T
  }
  
  OnMouseMotion <- function( x, y )
  {
    #Empty method but must be delcared to allow the real-time roi drawing on the MS image
  }

  SaveImg2Png <- function( ... )
  {
    mass_sel <- c()
    tol_sel <-c ()
    if( getValue_coloredCheckBox(this$Btn_RedEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[1])
      tol_sel<-c(tol_sel, mz_tolerance[1])
    }
    if( getValue_coloredCheckBox(this$Btn_GreenEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[2])
      tol_sel<-c(tol_sel, mz_tolerance[2])
    }
    if( getValue_coloredCheckBox(this$Btn_BlueEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[3])
      tol_sel<-c(tol_sel, mz_tolerance[3])
    }
    if(length(mass_sel) == 0)
    {
      gWidgets2::gmessage("No channel enabled, nothing is exported.", icon = "info")
      return()
    }


    fname<-gWidgets2::gfile("Save current MSI plot to png file", type="save", multi = F, filter =  c("tiff"="tiff", "svg"="svg"), initial.dir = path.expand(getwd()))
    if(length(fname) == 0)
    {
      return ()
    }

    #Auto append the image file extension
    if(grepl(".svg", basename(fname)))
    {
      fname<-paste(fname, ".svg", sep = "")
      svg( filename = fname , width = 1200, height = 500)
    }
    else
    {
      if(!grepl(".tiff", basename(fname)))
      {
        fname<-paste(fname, ".tiff", sep = "")
      }
      tiff( filename = fname , width = 1200, height = 500, compression = "none", res = 160)
    }

    plotMassImageByPeak(this$img,  mass.peak = mass_sel, tolerance = tol_sel,
                        XResLevel = switch(gWidgets2::svalue(this$Combo_Xres), x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5),
                        rotation = this$Rotation, vlight= gWidgets2::svalue(this$Scale_light),
                        crop_area = ZOOM_win, intensity_limit = this$IntLimits, NormalizationCoefs = this$NormalizationCoefs)
    dev.off()

    gWidgets2::gmessage(paste("Image saved at:", fname), icon = "info")
  }

  SpinImageRangeChanged <- function( ... )
  {
    #Set ROI from spinbuttons
    this$ROI[1]<- svalue(this$Spin_Xmin)
    this$ROI[2] <- svalue(this$Spin_Xmax)
    this$ROI[3] <- svalue(this$Spin_Ymin)
    this$ROI[4] <- svalue(this$Spin_Ymax)
    redrawPlotWidget(this$imaging_dev)
  }

  SliderLightChanged<- function( ... )
  {
    #Set it to ROI spinbuttons
    this$PlotMassImageRGB()
  }

  BtnRotateCCW <- function ( ... )
  {
    this$Rotation <- this$Rotation + 90
    if(this$Rotation == 360)
    {
      this$Rotation <- 0
    }
    this$RotateImage(this$Rotation)
  }

  BtnRotateCW <- function ( ... )
  {
    this$Rotation <- this$Rotation - 90
    if(this$Rotation == -90)
    {
      this$Rotation <- 270
    }
    this$RotateImage(this$Rotation)
  }

  RotateImage <- function( angle )
  {
    rotateLabel <- this$Rotation
    if( rotateLabel == 90 ) {rotateLabel<-270}
    else if(rotateLabel == 270){ rotateLabel<-90}
    Lbl_Rotation$set_value(paste("Rotation:", rotateLabel))

    #Plot rotated image
    redrawPlotWidget(this$imaging_dev)
  }

  ComboBox_XRes_Changed <- function( ... )
  {
    this$PlotMassImageRGB()
  }

  ComboBox_Norm_Changed <- function( ... )
  {
    if(gWidgets2::svalue(this$Combo_Norm) == "RAW")
    {
      this$NormalizationCoefs <- rep(1, nrow(this$img$pos))
    }
    else
    {
      this$NormalizationCoefs <- this$img$normalization[[gWidgets2::svalue(this$Combo_Norm)]]
    }

    #Re-Build the raster with the new Normalization factor
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

    #Re-Plot spectra using norm factors
    if(!is.null(this$GetSpectraInfo_ptr))
    {
      plotted_data <- this$GetSpectraInfo_ptr(this$myName)

      if(!is.null(this$ClearSpectraPlot_ptr))
      {
        this$ClearSpectraPlot_ptr(plotted_data$ID, this$myName)
        this$plotSpectra(plotted_data$ID, plotted_data$color, mz_min = plotted_data$mz_min, mz_max = plotted_data$mz_max)
      }
    }
  }

  ROI_Deleted <-function (...)
  {
    #Set it to ROI spinbuttons
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()
    svalue(this$Spin_Xmin) <- 1
    svalue(this$Spin_Xmax) <- img$size["x"]
    svalue(this$Spin_Ymin) <- 1
    svalue(this$Spin_Ymax) <- img$size["y"]
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    this$ROI <- NULL
    redrawPlotWidget(this$imaging_dev)
    gWidgets2::enabled(this$Frame_RoiCtl) <- F
    if(is.null(this$ZOOM_win))
    {
      gWidgets2::enabled(this$Btn_RoiZoom) <- F
    }
  }

  ROI_Zoom <- function( ... )
  {
    this$ZOOM_win <- switch(svalue(this$Btn_RoiZoom) , this$ROI)
    redrawPlotWidget(this$imaging_dev)

    if( is.null(this$ZOOM_win) )
    {
      this$Btn_RoiZoom[] <-"Zoom in ROI"
    }
    else
    {
      this$Btn_RoiZoom[] <-"Zoom Out"
    }

    #No roi defined and zoom out
    if(is.null(this$ZOOM_win) && is.null(this$ROI))
    {
      gWidgets2::enabled(this$Btn_RoiZoom) <- F
    }
  }

  ROI_GetSpectra <- function( ... )
  {
    if(!is.null(this$ROI))
    {
      Left <- this$ROI[1]
      Right <- this$ROI[2]
      Bottom <- this$ROI[3]
      Top <- this$ROI[4]

      #Limits to image size
      Left<-max(1, Left)
      Right<-max(1, Right)
      Bottom<-max(1, Bottom)
      Top<-max(1, Top)

      Left<-min(this$img$size["x"], Left)
      Right<-min(this$img$size["x"], Right)
      Bottom<-min(this$img$size["y"], Bottom)
      Top<-min(this$img$size["y"], Top)

      ID<-c()
      X<-c()
      Y<-c()
      Colour<-c()
      Zpos <-complex(real = this$img$pos[,"x"], imaginary = this$img$pos[,"y"]) #Convert positions to a complex numbers vector to find pointed coords fast and easy
      currID <- data.frame(this$Tbl_spotList$get_items())$ID
      for(xi in Left:Right)
      {
        for(yi in Bottom:Top)
        {
          preID <- which( Zpos == complex(real = xi, imaginary = yi) )

          if(length(preID) > 0)
          {
            if(!(preID %in% currID))
            {
              X<-c(X, xi)
              Y<-c(Y, yi)
              Colour<- c(Colour ,hsv( h = preID/length(Zpos), s = 0.7, v = 1))
              ID<-c(ID, preID)
            }
          }
        }
      }
      this$Tbl_spotList$set_items(rbind(data.frame(this$Tbl_spotList$get_items()), data.frame(ID, X, Y, Colour)))
      
      tweaksSpectraTable(this$Tbl_spotList)
    }
  }
  
  tweaksSpectraTable <- function(tbl)
  {
    tktable <- gWidgets2::getToolkitWidget(tbl)
    all_child <- as.character(tcltk::tcl(tktable, "children", ""))
    
    my_color_tag_list <-c()
    for(curr_child in all_child)
    {
      curr_color <- as.character( tcltk::tcl(tktable, "item", curr_child, "-values"))[4]
      tcltk::tcl(tktable, "item", curr_child, "-tags", curr_color)
      my_color_tag_list <-unique( c(my_color_tag_list, curr_color)  )
    }
    
    for(curr_tag_color in my_color_tag_list)
    {
      tcltk::tcl(tktable, "tag", "configure", curr_tag_color, "-background", curr_tag_color)
    }
    
    #Hide Color columns & teaks
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "configure","-displaycolumns", 1:3)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 4, "-minwidth", 0)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 4, "-width", 0)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 4, "-stretch", 0)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 1, "-stretch", 1)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 2, "-minwidth", 40)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 3, "-minwidth", 40)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 2, "-width", 40)
    tcltk::tcl(gWidgets2::getToolkitWidget(tbl), "column", 3, "-width", 40)
  }

  ROI_IntensityLimit <- function( ... )
  {
    if( is.null(this$ROI))
    {
      gWidgets2::gmessage("To apply intensity limitation you must define a ROI ", icon = "info")
      return()
    }
    else
    {
      this$IntLimit_ROI <- this$ROI
      #Limit Roi to image range
      this$IntLimit_ROI[1] <- max( c(this$IntLimit_ROI[1], 1 ) )
      this$IntLimit_ROI[2] <- min( c(this$IntLimit_ROI[2], this$img$size["x"] ) )
      this$IntLimit_ROI[3] <- max( c(this$IntLimit_ROI[3], 1 ) )
      this$IntLimit_ROI[4] <- min( c(this$IntLimit_ROI[4], this$img$size["y"] ) )
    }

    #Re-Build the raster with the new intensity limit
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

    gWidgets2::enabled(this$Btn_RoiIntUnLimit) <- T
    gWidgets2::svalue(this$Btn_RoiIntUnLimit)<- "Remove Intensity Limit"
    this$Btn_RoiIntUnLimit$add_handler_changed(this$ROI_IntensityUnLimit)
  }

  ROI_IntensityUnLimit <- function( ... )
  {
    #Disable intensity limitation
    this$IntLimit_ROI <- NULL

    #Re-Build the raster with the new intensity limit
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

    this$Btn_RoiIntUnLimit$remove_handlers()
    gWidgets2::enabled(this$Btn_RoiIntUnLimit) <- F
    gWidgets2::svalue(this$Btn_RoiIntUnLimit)<- "No Intensity Limit"
  }

  HideShowSpectraList <- function ( ... )
  {
    if(gWidgets2::svalue(this$Btn_ShowSpectraList))
    {
      #show spectra list
      setSpectraListWidth(this$previous_spectralist_width)
    }
    else
    {
      #hide spectra list
      curr_widget_width <- as.numeric(tcltk::tkwinfo("width", gWidgets2::getToolkitWidget(this$Top_frm)))
      this$previous_spectralist_width <- gWidgets2::svalue(this$Panel_Img) * curr_widget_width
      gWidgets2::svalue(this$Panel_Img) <- 0
    }
  }

  BtnFlipH <- function ( ... )
  {
    this$flipH <-  !(this$flipH)
    redrawPlotWidget(this$imaging_dev)
  }

  BtnFlipV <- function ( ... )
  {
    this$flipV <- !(this$flipV)
    redrawPlotWidget(this$imaging_dev)
  }
  
  onExpose <- function()
  {
    if(gWidgets2::svalue(this$Btn_ShowSpectraList))
    {
      setSpectraListWidth(this$previous_spectralist_width)
    }
  }
  
  setSpectraListWidth <- function(desired_width_pixels)
  {
    curr_widget_width <- as.numeric(tcltk::tkwinfo("width", gWidgets2::getToolkitWidget(this$Top_frm)))
    gWidgets2::svalue(this$Panel_Img) <- this$previous_spectralist_width /  curr_widget_width# % of space allocated to spectra list  
  }

  #Build the GUI
  Top_frm <- gWidgets2::gframe( text =  img$name, container = parent, fill = T, expand = T, spacing = 2 ) 
  tcltk::tkbind( gWidgets2::getToolkitWidget(Top_frm), "<Expose>", this$onExpose )
  
  Panel_Img<- gWidgets2::gpanedgroup(horizontal = T, container = Top_frm,  fill = T, expand = T )
  spectraListFrame<-gWidgets2::gframe("Spectra List", container = Panel_Img,  fill = T, spacing = 5, expand = T )
  Grp_Tbl <- gWidgets2::ggroup(horizontal = F, container = spectraListFrame,  expand=TRUE, fill = TRUE)
  
  ID<-0
  X<-0
  Y<-0
  Tbl_spotList<-gWidgets2::gtable( data.frame(ID,X,Y, Colour = meanSpectrumColor), container = Grp_Tbl, multiple = T, chosen.col = 1,  fill = T, expand = T )
  tweaksSpectraTable(Tbl_spotList)
  
  Btn_PlotSelSpotList<-gWidgets2::gbutton("Plot", container= Grp_Tbl,  handler = this$SpectraListSelChange)
  Btn_ClearSpotList<-gWidgets2::gbutton("Clear", container= Grp_Tbl,  handler = this$BtnClearSpotList)
  Btn_ExportSpotList<-gWidgets2::gbutton("Export", container= Grp_Tbl,  handler = this$BtnExportSpotList)

  Grp_TopImg <- gWidgets2::ggroup(horizontal = F, container = Panel_Img, expand = T, fill = T)
  Grp_Buttons <- gWidgets2::ggroup(horizontal = T, container = Grp_TopImg, expand = F, fill = F)

  Btn_ShowSpectraList <- gWidgets2::gcheckbox("Spectra list", container = Grp_Buttons, handler = this$HideShowSpectraList, use.togglebutton = T, checked = T )
  Lbl_Rotation<- gWidgets2::glabel(text = "Rotation: 0", container = Grp_Buttons)
  Btn_rotate_CCW <- gbutton_icon(image_file = file.path(system.file(package = "rMSI2", "icons"),"Rotate_CCW.png"), container = Grp_Buttons, handler = this$BtnRotateCCW)
  Btn_rotate_CW <- gbutton_icon(image_file = file.path(system.file(package = "rMSI2", "icons"),"Rotate_CW.png"), container = Grp_Buttons, handler = this$BtnRotateCW)
  Btn_flipV <-  gbutton_icon(image_file = file.path(system.file(package = "rMSI2", "icons"),"FlipV.png"), container = Grp_Buttons, handler = this$BtnFlipV)
  Btn_flipH <- gbutton_icon(image_file = file.path(system.file(package = "rMSI2", "icons"),"FlipH.png"), container = Grp_Buttons, handler = this$BtnFlipH)
  
  Lbl_Xres<- gWidgets2::glabel(text = "Interpolation:", container = Grp_Buttons)
  Combo_Xres <- gcombobox_rMSI( items = c("x1","x2","x3","x4","x5"), selected = 2, container = Grp_Buttons, handler = this$ComboBox_XRes_Changed)
  Lbl_Normalitzation <- gWidgets2::glabel(text = "Normalization:", container = Grp_Buttons)
  Combo_Norm <- gcombobox_rMSI( items = c("RAW", names(img$normalizations)), selected = 1, container = Grp_Buttons, handler = this$ComboBox_Norm_Changed)
  gWidgets2::glabel("Light:", container = Grp_Buttons)
  Scale_light <- gWidgets2::gslider( from = 0.6, to = 10, by = 0.2, value = 3, horizontal = T, handler = this$SliderLightChanged, container =  Grp_Buttons)
  gWidgets2::addSpring(Grp_Buttons)
  Btn_plot2file<- gWidgets2::gbutton("Save in image file", container = Grp_Buttons, handler = this$SaveImg2Png)

  Grp_ImgTop<-gWidgets2::ggroup( horizontal = T, container =  Grp_TopImg,  fill = T, expand = T)
  Grp_ImgRoi<-gWidgets2::ggroup( horizontal = F, container =  Grp_ImgTop,  fill = T, expand = T)
  
  imaging_dev <-createPlotWidget(parent = Grp_ImgRoi, 
                     redraw_function = this$RedrawMSImage,
                     MouseSelection_callback = this$OnPixelSelection,
                     MouseHover_callback = this$OnMouseMotion,
                     initial_width = 200, 
                     initial_height = 200 ) 
  
  Grp_ScalesV <- gWidgets2::ggroup( horizontal = F, container =  Grp_ImgTop,  fill = T, expand = F) 
  Grp_ScalesH <- gWidgets2::ggroup( horizontal = T, container =  Grp_ScalesV,  fill = T, expand = T, spacing = 0)
  Grp_ZoomLimit<-gWidgets2::ggroup(horizontal = T, container = Grp_ScalesV)
  Btn_RoiZoom<-gWidgets2::gcheckbox("Zoom in ROI", checked = F, use.togglebutton = T, container = Grp_ZoomLimit, handler = this$ROI_Zoom)
  Btn_RoiIntUnLimit<-gWidgets2::gbutton("No Intensity Limit", container = Grp_ZoomLimit)

  #Red Color Scale
  Grp_RedScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleRed_dev <-createPlotWidget(parent = Grp_RedScale, 
                                 redraw_function = this$ReDrawRedScale,
                                 initial_width = 100, 
                                 initial_height = 200 ) 
 
  Btn_RedEnable <- coloredCheckBox(text  = "m/z", checked = T, handler =  this$PlotMassImageRGB, container = Grp_RedScale, foreground = "red", bold = T)
 
  #Green Color scale
  Grp_GreenScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleGreen_dev <-createPlotWidget(parent = Grp_GreenScale, 
                                  redraw_function = this$ReDrawGreenScale,
                                  initial_width = 100, 
                                  initial_height = 200 ) 
  
  Btn_GreenEnable <- coloredCheckBox(text  = "m/z", checked = F, handler =  this$PlotMassImageRGB, container = Grp_GreenScale, foreground = "darkgreen", bold = T)
  
  #Blue Color scale
  Grp_BlueScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleBlue_dev <-createPlotWidget(parent = Grp_BlueScale, 
                                  redraw_function = this$ReDrawBlueScale,
                                  initial_width = 100, 
                                  initial_height = 200 ) 
  
  Btn_BlueEnable <- coloredCheckBox(text  = "m/z", checked = F, handler =  this$PlotMassImageRGB, container = Grp_BlueScale, foreground = "darkblue", bold = T)  
  
  #ROI CTL
  Frame_RoiCtl<-gWidgets2::gframe("ROI Controls", container = Grp_ImgRoi, expand = F, fill = F, spacing = 1 )
  Grp_RoiCtl<-gWidgets2::ggroup(horizontal = F, container = Frame_RoiCtl, expand = F, fill = F, spacing = 1)
  Grp_RoiCtl_1stRow<-gWidgets2::ggroup(horizontal = T, container = Grp_RoiCtl)
  Grp_RoiCtl_2ndRow<-gWidgets2::ggroup(horizontal = T, container = Grp_RoiCtl)
  Btn_RoiDelete<-gWidgets2::gbutton("Delete", container = Grp_RoiCtl_1stRow, handler = this$ROI_Deleted)
  Btn_RoiGetSpectra<-gWidgets2::gbutton("Get Spectra", container = Grp_RoiCtl_1stRow, handler = this$ROI_GetSpectra)
  gWidgets2::addSpring(Grp_RoiCtl_1stRow)
  Lbl_XImgRange<- gWidgets2::glabel(text = "X range:", container = Grp_RoiCtl_1stRow)
  Spin_Xmin<- gWidgets2::gspinbutton(from = 1, to =  img$size["x"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_1stRow)
  Spin_Xmax<- gWidgets2::gspinbutton(from = 1, to = img$size["x"], digest = 0, by = 1 , value = img$size["x"], handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_1stRow)
  Btn_RoiIntLimit<-gWidgets2::gbutton("Apply Intensity Limit", container = Grp_RoiCtl_2ndRow, handler = this$ROI_IntensityLimit)
  gWidgets2::addSpring(Grp_RoiCtl_2ndRow)
  Lbl_YImgRange<- gWidgets2::glabel(text = "Y range:", container = Grp_RoiCtl_2ndRow)
  Spin_Ymin<- gWidgets2::gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_2ndRow)
  Spin_Ymax<- gWidgets2::gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = img$size["y"], handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_2ndRow)

  gWidgets2::enabled(Btn_RoiZoom) <- F
  gWidgets2::enabled(Frame_RoiCtl) <- F
  gWidgets2::enabled(Btn_RoiIntUnLimit) <- F
  
  # Set the name for the class
  class(this) <- append(class(this),"MSImagePlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
