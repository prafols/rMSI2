###A GUI to display MS images from some ions/features

.MSImagePlotWidget <- function( in_img, parent_widget=gwindow ( "Default MSImagePlotWidget" , visible = FALSE ), AddSpectra_function = NULL)
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
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
  imaging_dev <- 0
  Spin_Smooth <- 0
  Spin_Xmin <- 0
  Spin_Xmax <- 0
  Spin_Ymin <- 0
  Spin_Ymax <- 0
  Rotation <- 0
  Tbl_spotList <- 0
  iSel <- 1:nrow(img$pos) #Initially grab the whole image
  image_range <- c( 1, 1, img$size["x"], img$size["y"] )
  AddSpectra_ptr <- AddSpectra_function #Pointer to a AddSpectra method of a spectraWidget to be able of ploting directly

  #Current image RGB layers
  Rlayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Glayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Blayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )


  #==================================================================================================
  ImgBuildFun <- function(channel, mass, tol )
  {
    this$mz_tolerance[channel] <- tol
    this$mz_selected[channel] <- mass
    img_new<-.buildImageByPeak(this$img, mass.peak = mass, tolerance = tol, NormCoefs = NULL ) #TODO some day I can use this to include various normalizations

    if( channel == 1)
    {
      this$Rlayer_raster<-img_new
    }
    else if (channel == 2)
    {
      this$Glayer_raster<-img_new
    }
    else if (channel == 3)
    {
      this$Blayer_raster<-img_new
    }

    #Return del buildImage
    return(list(selMz = img_new$mass, selTol = img_new$tolerance))
  }

  #==================================================================================================
  PlotMassImageRGB <- function()
  {
    ch_count <- 0
    if( svalue(this$Btn_RedEnable))
    {
      red_layer <- this$Rlayer_raster
      unique_layer <- red_layer
      ch_count <- ch_count + 1
    }
    else
    {
      red_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( svalue(this$Btn_GreenEnable))
    {
      green_layer <- this$Glayer_raster
      unique_layer <- green_layer
      ch_count <- ch_count + 1
    }
    else
    {
      green_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( svalue(this$Btn_BlueEnable))
    {
      blue_layer <- this$Blayer_raster
      unique_layer <- blue_layer
      ch_count <- ch_count + 1
    }
    else
    {
      blue_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }

    if(ch_count < 1)
    {
      print("No selected data to plot image")
      return()
    }
    if(ch_count == 1)
    {
      plotting_raster<-.BuildSingleIonRGBImage( unique_layer,   XResLevel = 1 )##TODO now the XResLvel is set with an integer addapt it
    }
    else
    {
      plotting_raster<-.BuildRGBImage( imgR = red_layer, imgG = green_layer, imgB = blue_layer, XResLevel = 1 )##TODO now the XResLvel is set with an integer addapt it
    }

    visible(this$imaging_dev)<-TRUE
    .plotMassImageRGB (plotting_raster, cal_um2pixels = this$img$pixel_size_um,  rotation = this$Rotation, display_axes = F)


    visible(this$scaleRed_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(red_layer)
    }
    else
    {
      .plotIntensityScale(red_layer, "R")
    }

    visible(this$scaleGreen_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(green_layer)
    }
    else
    {
      .plotIntensityScale(green_layer, "G")
    }


    visible(this$scaleBlue_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(blue_layer)
    }
    else
    {
      .plotIntensityScale(blue_layer, "B")
    }


    ###TODO this is the old implementation to b remove when the new method works
    #plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$Combo_Xres), selectedPixels = this$iSel, rotation = this$Rotation)
  }

  #==================================================================================================
  BtnClearSpotList <- function( mass, tol, ... )
  {
    this$spectraWidget$OnLostFocus() #This is just a test

    this$Tbl_spotList$set_items(data.frame(this$Tbl_spotList$get_items())[1,])
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0),
                               gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]],
                               background = 3 )

    render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]]
    render$set( font = "bold")
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0), render)
    gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 3)$set(visible = F)
  }

  #==================================================================================================
  SpectraListSelChange <- function( ... )
  {
    selected <- svalue(this$Tbl_spotList)
    this$spectraWidget$ClearSpectra()
    df<-data.frame(this$Tbl_spotList$get_items())  #get data frame...
    max_nrow<-nrow(this$img$data[[1]])

    for( i in selected)
    {
      selDf<-df[ df$ID == i,] #Get the correct data row
      if(i == 0)
      {
        this$spectraWidget$AddSpectra(this$spectra_mass, this$spectra_intensity, col = as.character(selDf$Colour))
      }
      else
      {
        icube<-(1+((i-1) %/% max_nrow))
        irow<- (i - (icube -1) * max_nrow)
        this$spectraWidget$AddSpectra(this$img$mass, this$img$data[[icube]][ irow ,] , col =  as.character(selDf$Colour))
      }
    }
  }

  #==================================================================================================
  BtnExportSpotList <- function( ... )
  {
    indexes<-as.vector(data.frame(this$Tbl_spotList$get_items())$ID)
    indexes<-indexes[indexes > 0] #Remove zero!
    store_paths<-gfile("Save current spots to txt files", type="selectdir", multi = F, initial.dir = path.expand("~/"))
    if(length(store_paths) == 0)
    {
      return ()
    }

    #display a BIG WaRNIng IF to much ID's are seected!
    bExportData <- T
    if( length(indexes) > 25)
    {
      bExportData<-gconfirm(paste("You are exporting a lot of data. (", length(indexes) ,"mass spectrums )\nThis may take a long time and expend a lot of memory.\nDo you want to store spectra as TXT?"), title = "Warning: Large data export!", icon = "warning")
    }

    mPbar<-.ProgressBarDialog("Exporting spectra to txt...")

    #Create a dir to store all data inside
    store_paths <- file.path(store_paths, paste("Export_", tools::file_path_sans_ext(this$img$name),"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), sep = "" ))
    dir.create( store_paths )


    #Save mean spectra
    spc <- matrix(data= c(this$img$mean@mass, this$img$mean@intensity), ncol = 2, byrow = F)
    write( x = t(spc), file = file.path( store_paths, "average.txt" ), ncolumns = 2  )


    if(length(indexes) > 0)
    {
      #Save the id file
      id_fname <- file.path( store_paths, "ID.txt")
      write(indexes, file = id_fname, ncolumns = 1)

      if(bExportData)
      {
        #Save each ID in the lit
        dataChunck <- loadImgCunckFromIds(this$img, indexes)

        for( i in 1:length(indexes))
        {
          spc <- matrix(data= c(this$img$mass, dataChunck[i, ]), ncol = 2, byrow = F)
          spc_fname<-  file.path( store_paths, paste("ID_", indexes[i] ,".txt", sep =""))
          write( x = t(spc), file = spc_fname, ncolumns = 2  )
          if( !mPbar$setValue( 100* i/length(indexes) ))
          {
            #Aborted by user
            rm(dataChunk)
            unlink( store_paths, recursive = T )
            gc()
            return()
          }
        }
      }
    }

    mPbar$close()

    rm(dataChunk)
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
    gmessage(smsg, title = "Export complete!", icon = "info")

  }

  #==================================================================================================
  OnPixelSelection <- function( evt, ...)
  {
    xResDiv<-switch(svalue(this$Combo_Xres), x1 = 1, x2 = 2, x4 = 4)

    X_left<-round(min(evt$x)/xResDiv)
    X_right<-round(max(evt$x)/xResDiv)
    Y_bottom<-round(min(evt$y)/xResDiv)
    Y_top<-round(max(evt$y)/xResDiv)

    #Apply rotation
    if(this$Rotation >= 0 && this$Rotation < 90)
    {
      Left <- X_left
      Right <- X_right
      Bottom <- Y_bottom
      Top <- Y_top
    }
    if(this$Rotation >= 90 && this$Rotation < 180)
    {
      Left <- Y_bottom
      Right <- Y_top
      Bottom <- this$img$size["y"] - X_right
      Top <- this$img$size["y"] - X_left
    }
    if(this$Rotation >= 180 && this$Rotation < 270)
    {
      Left <- this$img$size["x"] - X_right
      Right <- this$img$size["x"] - X_left
      Bottom <- this$img$size["y"] - Y_top
      Top <- this$img$size["y"] - Y_bottom
    }
    if(this$Rotation >= 270 && this$Rotation < 360)
    {
      Left <- this$img$size["x"] - Y_top
      Right <- this$img$size["x"] - Y_bottom
      Bottom <- X_left
      Top <- X_right
    }

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
        yimg <- 1 + this$img$size["y"] - yi
        preID <- which( Zpos == complex(real = xi, imaginary = yimg) )
        if(!(preID %in% currID))
        {
          X<-c(X, xi)
          Y<-c(Y, yi)
          Colour<- c(Colour ,hsv( h = preID/length(Zpos), s = 0.7, v = 1))
          ID<-c(ID, preID)
        }
      }
    }
    this$Tbl_spotList$set_items(rbind(data.frame(this$Tbl_spotList$get_items()), data.frame(ID, X, Y, Colour)))
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0),
                               gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]],
                               background = 3
    )

    render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]]
    render$set( font = "bold")
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0), render)

    gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 3)$set(visible = F)
  }

  #==================================================================================================
  SaveImg2Png <- function( ... )
  {
    mass_sel <- c()
    tol_sel <-c ()
    if( svalue(this$Btn_RedEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[1])
      tol_sel<-c(tol_sel, mz_tolerance[1])
    }
    if( svalue(this$Btn_GreenEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[2])
      tol_sel<-c(tol_sel, mz_tolerance[2])
    }
    if( svalue(this$Btn_BlueEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[3])
      tol_sel<-c(tol_sel, mz_tolerance[3])
    }
    if(length(mass_sel) == 0)
    {
      gWidgets2::gmessage("No channel enabled, nothing is exported.", icon = "info")
      return()
    }


    fname<-gfile("Save current MSI plot to png file", type="save", multi = F, filter =  c("tiff"="tiff"), initial.dir = path.expand(getwd()))
    if(length(fname) == 0)
    {
      return ()
    }

    #Auto append the image file extension
    if(!grepl(".tiff", basename(fname)))
    {
      fname<-paste(fname, ".tiff", sep = "")
    }

     #This was the old implementation
#     visible(this$imaging_dev)<-TRUE
#     dev.print(jpeg, filename = fname, quality = "99", width = size(this$imaging_dev)["width"], height = size(this$imaging_dev)["height"])

    #New implementation
    tiff( filename = fname , width = 1200, height = 500, compression = "none", res = 160)
    plotMassImageByPeak(this$img,  mass.peak = mass_sel, tolerance = tol_sel, XResLevel = 3, rotation = this$Rotation)
    dev.off()
  }

  #==================================================================================================
  SpinSmoothChanged <- function( ... )
  {
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  SpinImageRangeChanged <- function( ... )
  {
    #Grab new roi to plot
    X_left<-svalue(this$Spin_Xmin)
    X_right<-svalue(this$Spin_Xmax)
    Y_bottom<-svalue(this$Spin_Ymin)
    Y_top<-svalue(this$Spin_Ymax)

    #Apply rotation!
    if(this$Rotation >= 0 && this$Rotation < 90)
    {
      x_min <- X_left
      x_max <- X_right
      y_min <- Y_bottom
      y_max <- Y_top
    }
    if(this$Rotation >= 90 && this$Rotation < 180)
    {
      x_min <- Y_bottom
      x_max <- Y_top
      y_min <- this$img$size["y"] - X_right
      y_max <- this$img$size["y"] - X_left
    }
    if(this$Rotation >= 180 && this$Rotation < 270)
    {
      x_min <- this$img$size["x"] - X_right
      x_max <- this$img$size["x"] - X_left
      y_min <- this$img$size["y"] - Y_top
      y_max <- this$img$size["y"] - Y_bottom
    }
    if(this$Rotation >= 270 && this$Rotation < 360)
    {
      x_min <- this$img$size["x"] - Y_top
      x_max <- this$img$size["x"] - Y_bottom
      y_min <- X_left
      y_max <- X_right
    }

    #Keep a copy of image limits
    this$image_range <- c( x_min, x_max, y_min, y_max )

    #Transform Y coord to fit img space
    y_min_aux<-y_min
    y_min<-1 + this$img$size["y"] -y_max
    y_max<-1 + this$img$size["y"] -y_min_aux

    iSel <- ((this$img$pos[, "x"] >= x_min) & (this$img$pos[, "x"] <= x_max))
    iSel <- iSel & ((this$img$pos[, "y"] >= y_min) & (this$img$pos[, "y"] <= y_max))
    this$iSel <- which(iSel)

    # DEBUG PRINTS
    #   cat("\nImage Range Changed!\n")
    #   cat(paste("x_min = ", x_min, "\n"))
    #   cat(paste("x_max = ", x_max, "\n"))
    #   cat(paste("y_min = ", y_min, "\n"))
    #   cat(paste("y_max = ", y_max, "\n"))
    #   cat("iSel:\n")
    #   print(this$.iSel)

    this$PlotMassImageRGB()
  }

  #==================================================================================================
  BtnRotateCCW <- function ( ... )
  {
    this$Rotation <- this$Rotation + 90
    if(this$Rotation == 360)
    {
      this$Rotation <- 0
    }
    this$RotateImage(this$Rotation)
  }

  #==================================================================================================
  BtnRotateCW <- function ( ... )
  {
    this$Rotation <- this$Rotation - 90
    if(this$Rotation == -90)
    {
      this$Rotation <- 270
    }
    this$RotateImage(this$Rotation)
  }

  #==================================================================================================
  RotateImage <- function( angle )
  {
    rotateLabel <- this$Rotation
    if( rotateLabel == 90 ) {rotateLabel<-270}
    else if(rotateLabel == 270){ rotateLabel<-90}
    Lbl_Rotation$set_value(paste("Rotation:", rotateLabel, "º"))

    #Adjust range controls to fit rotation
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()

    if(angle >= 0 && angle < 90 || angle >= 180 && angle < 270)
    {
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmin), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmax), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymin), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymax), min = 1, max = this$img$size["y"] )
    }
    if(angle >= 90 && angle < 180 || angle >= 270 && angle < 360)
    {
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmin), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmax), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymin), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymax), min = 1, max = this$img$size["x"] )
    }

    #Apply rotation to image range
    if(angle >= 0 && angle < 90)
    {
      svalue(this$Spin_Xmin) <- this$image_range[1]
      svalue(this$Spin_Xmax) <- this$image_range[2]
      svalue(this$Spin_Ymin) <- this$image_range[3]
      svalue(this$Spin_Ymax) <- this$image_range[4]
    }
    if(angle >= 90 && angle < 180)
    {
      svalue(this$Spin_Xmin) <- this$img$size["y"] - this$image_range[4] + 1
      svalue(this$Spin_Xmax) <- this$img$size["y"] - this$image_range[3] + 1
      svalue(this$Spin_Ymin) <- this$image_range[1]
      svalue(this$Spin_Ymax) <- this$image_range[2]
    }
    if(angle >= 180 && angle < 270)
    {
      svalue(this$Spin_Xmin) <- this$img$size["x"] - this$image_range[2] + 1
      svalue(this$Spin_Xmax) <- this$img$size["x"] - this$image_range[1] + 1
      svalue(this$Spin_Ymin) <- this$img$size["y"] - this$image_range[4] + 1
      svalue(this$Spin_Ymax) <- this$img$size["y"] - this$image_range[3] + 1
    }
    if(angle >= 270 && angle < 360)
    {
      svalue(this$Spin_Xmin) <- this$image_range[3]
      svalue(this$Spin_Xmax) <- this$image_range[4]
      svalue(this$Spin_Ymin) <- this$img$size["x"] - this$image_range[2] + 1
      svalue(this$Spin_Ymax) <- this$img$size["x"] - this$image_range[1] + 1
    }

    #Re-connect handlers
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    #Plot rotated image
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  ComboBox_XRes_Changed <- function( ... )
  {
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  ROI_Deleted <-function (...)
  {
    ##TODO
  }

  #==================================================================================================
  ROI_Cropped <- function( ... )
  {
    ##TODO
  }

  #==================================================================================================
  ROI_ZoomIn <- function( ... )
  {
    ##TODO
  }

  #==================================================================================================
  ROI_ZoomOut <- function( ... )
  {
    ##TODO
  }

  #==================================================================================================
  ROI_GetSpectra <- function( ... )
  {
    ##TODO
  }

  #==================================================================================================
  IntensityScale_EnableClicked <- function( evt, ...)
  {
    ##TODO remove I dont need it, remove also from the handler action connection
    ###channel <- evt$action #"R" for the red image, "G" for the green image or "Blue" for the blue image

    this$PlotMassImageRGB()

  }

  #Build the GUI
  Top_grp <- ggroup(horizontal = F, container = parent)
  Top_captionGrp <- ggroup(horizontal = T, container = Top_grp)
  lbl_title <- glabel( paste("  Image:",img$name) , container = Top_captionGrp)
  font(lbl_title)<-list(weight = "bold", size = 9)
  addSpring(Top_captionGrp)
  Btn_plot2file<- gbutton("Save to jpeg", container = Top_captionGrp, handler = this$SaveImg2Png)
  Panel_Img<- gpanedgroup(horizontal = T, container = Top_grp,  expand=TRUE )
  spectraListFrame<-gframe("Spectra List", container = Panel_Img,  fill = T, spacing = 5 )
  Grp_Tbl <- ggroup(horizontal = F, container = spectraListFrame,  expand=TRUE, fill = TRUE)
  ID<-0
  X<-0
  Y<-0
  Colour<-"red" #default colour for mean spectra
  Tbl_spotList<-gtable( data.frame(ID,X,Y,Colour), container = Grp_Tbl, multiple = T, chosen.col = 1)
  size( Tbl_spotList )<- c(120, -1)
  ##Set table style using colors
  gtkTreeViewSetGridLines(getToolkitWidget(Tbl_spotList), 3)
  gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0),
                             gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]],
                             background = 3
  )

  render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]]
  render$set( font = "bold")
  gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0), render)

  gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 3)$set(visible = F)

  Grp_BtmTbl <-ggroup(horizontal = F, container =Grp_Tbl)
  Btn_PlotSelSpotList<-gbutton("Plot", container= Grp_BtmTbl,  handler = this$SpectraListSelChange)
  Btn_ClearSpotList<-gbutton("Clear", container= Grp_BtmTbl,  handler = this$BtnClearSpotList)
  Btn_ExportSpotList<-gbutton("Export", container= Grp_BtmTbl,  handler = this$BtnExportSpotList)

  Grp_TopImg <- ggroup(horizontal = F, container = Panel_Img, expand = T)
  Grp_Buttons <- ggroup(horizontal = T, container = Grp_TopImg)
  Lbl_Rotation<- glabel(text = "Rotation: 0º", container = Grp_Buttons)
  Btn_rotate_CCW <- gbutton("", container = Grp_Buttons, handler = this$BtnRotateCCW)
  gtkImageSetFromFile( getToolkitWidget(Btn_rotate_CCW)$image, filename = file.path(system.file(package = "rMSI", "icons"),"Rotate_CCW.png") )
  Btn_rotate_CW <- gbutton("", container = Grp_Buttons, handler = this$BtnRotateCW)
  gtkImageSetFromFile( getToolkitWidget(Btn_rotate_CW)$image, filename = file.path(system.file(package = "rMSI", "icons"),"Rotate_CW.png") )
  Lbl_Smooth<- glabel(text = "Smooth:", container = Grp_Buttons)
  Spin_Smooth<- gspinbutton(from = 0.1, to = 10, digest = 2, by = 0.2 , value = 0.3, handler = this$SpinSmoothChanged, container = Grp_Buttons)
  addSpring(Grp_Buttons)
  Lbl_Xres<- glabel(text = "Resolution:", container = Grp_Buttons)
  Combo_Xres <- gcombobox( items = c("x1","x2","x4"), selected = 2, container = Grp_Buttons, handler = this$ComboBox_XRes_Changed)

  Lbl_XImgRange<- glabel(text = "X range:", container = Grp_Buttons)
  Spin_Xmin<- gspinbutton(from = 1, to =  img$size["x"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Spin_Xmax<- gspinbutton(from = 1, to = img$size["x"], digest = 0, by = 1 , value = img$size["x"], handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Lbl_YImgRange<- glabel(text = "Y range:", container = Grp_Buttons)
  Spin_Ymin<- gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Spin_Ymax<- gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = img$size["y"], handler = this$SpinImageRangeChanged, container = Grp_Buttons)

  Grp_ImgTop<-ggroup( horizontal = T, container =  Grp_TopImg,  fill = T, expand = T)
  imaging_dev <- ggraphics(spacing = 5 )
  size( imaging_dev )<- c(650, 340)
  addHandlerSelectionChanged( imaging_dev, handler = this$OnPixelSelection, action = this)
  add(obj = Grp_ImgTop, child = imaging_dev,  fill = T, expand = T)

  #Red Color Scale
  Grp_RedScale<-ggroup( horizontal = F, container = Grp_ImgTop)
  scaleRed_dev <- ggraphics(spacing = 5 )
  size( scaleRed_dev )<- c(120, 300)
  add(obj = Grp_RedScale, child = scaleRed_dev,  fill = T, expand = T)
  Lbl_RedMz <- glabel("R M/Z", container = Grp_RedScale)
  Btn_RedEnable<-gcheckbox("On", container = Grp_RedScale, use.togglebutton = T, checked = T,  handler = this$IntensityScale_EnableClicked, action = "R")

  #Green Color scale
  Grp_GreenScale<-ggroup( horizontal = F, container = Grp_ImgTop)
  scaleGreen_dev <- ggraphics(spacing = 5 )
  size( scaleGreen_dev )<- c(120, 300)
  add(obj = Grp_GreenScale, child = scaleGreen_dev,  fill = T, expand = T)
  Lbl_GreenMz <- glabel("G M/Z", container = Grp_GreenScale)
  Btn_GreenEnable<-gcheckbox("On", container = Grp_GreenScale, use.togglebutton = T, checked = F, handler = this$IntensityScale_EnableClicked, action = "G")

  #Blue Color scale
  Grp_BlueScale<-ggroup( horizontal = F, container = Grp_ImgTop)
  scaleBlue_dev <- ggraphics(spacing = 5 )
  size( scaleBlue_dev )<- c(120, 300)
  add(obj = Grp_BlueScale, child = scaleBlue_dev,  fill = T, expand = T)
  Lbl_BlueMz <- glabel("B M/Z", container = Grp_BlueScale)
  Btn_BlueEnable<-gcheckbox("On", container = Grp_BlueScale, use.togglebutton = T, checked = F, handler = this$IntensityScale_EnableClicked, action = "B")

  #ROI CTL
  Frame_RoiCtl<-gframe("ROI Controls", container = Grp_TopImg )
  Grp_RoiCtl<-ggroup(horizontal = T, container = Frame_RoiCtl)
  Btn_RoiDelete<-gbutton("Delete", container = Grp_RoiCtl, handler = this$ROI_Deleted)
  Btn_RoiCrop<-gcheckbox("Crop", container = Grp_RoiCtl, handler = this$ROI_Croped, use.togglebutton = T, checked = F)
  Btn_RoiZoomIn<-gbutton("Zoom In", container = Grp_RoiCtl, handler = this$ROI_ZoomIn)
  Btn_RoiZoomOut<-gbutton("Zoom Out", container = Grp_RoiCtl, handler = this$ROI_ZoomOut)
  Btn_RoiGetSpectra<-gbutton("Get Spectra", container = Grp_RoiCtl, handler = this$ROI_GetSpectra)
  addSpring(Grp_RoiCtl)

  ## Set the name for the class
  class(this) <- append(class(this),"MSImagePlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}