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


#MALDI Image reconstruction by selected ion Top Windows
##########################################################

#' Open a MS image from Hdd directly allowing the user to chose it from a file dialog.
#'
#' A file choser dialog is presented to allow the selection of a .tar file containing an rMSI data object.
#' If a ramdisk has been previously created it will be used to speed up the loading process.
#' The image will be presented using the GUI and it will be returned as rMSI object.
#'
#' @param lockExecution if is set to true, R execution will be halted till GUI is closed.
#'
#' @return the rMSI object containing the MS image.
#'
#' @export
OpenMSI<-function( lockExecution = F )
{
  #Load data using GUI
  MSI_obj<-LoadTwoMsImages()

  #Display the MSIWindow
  if( !is.null(MSI_obj$img1_obj) && is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img1_obj, lockExecution = lockExecution)
  }
  if( is.null(MSI_obj$img1_obj) && !is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img2_obj, lockExecution = lockExecution)
  }
  if( !is.null(MSI_obj$img1_obj) && !is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img1_obj, img2 = MSI_obj$img2_obj, lockExecution = lockExecution)
  }
  
  return(list(img1 = MSI_obj$img1_obj, img2 = MSI_obj$img2_obj ))
}


#' Open the GUI to explore a MS image
#'
#' @param img a rMSI data object
#' @param lockExecution if is set to true, R execution will be halted till GUI is closed.
#'
#'  Open the GUI to explore the MS image. A MS image can be provided as a parameter.
#'  Up to two images can be displayed at once using multiple arguments.
#'
#' @export
MSIWindow<-function(img1, img2 = NULL, lockExecution = F)
{
  options(guiToolkit="tcltk") #force to use tcltk
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  spectraWidget <- NULL
  msiWidget1 <- NULL
  msiWidget2 <- NULL
  ChannelRedraw <- list( ReDraw = F, mass = 0, tol = 0  )
  MsiRedraw <- list( Red = ChannelRedraw, Green = ChannelRedraw, Blue = ChannelRedraw  )
  RedrawIsIdle <- F

  #Redraw MS image by timer callback
  ReDrawByTimer <- function(data)
  {
    if( !this$RedrawIsIdle)
    {
      this$RedrawIsIdle <- T
      spectraWidget$ReDrawByTimer(data)
      for( i in 1:3)
      {
        if( this$MsiRedraw[[i]]$ReDraw)
        {
          this$msiWidget1$ImgBuildFun(i, this$MsiRedraw[[i]]$mass, this$MsiRedraw[[i]]$tol)
          this$msiWidget1$PlotMassImageRGB()
          if(!is.null(this$msiWidget2))
          {
            this$msiWidget2$ImgBuildFun(i, this$MsiRedraw[[i]]$mass, this$MsiRedraw[[i]]$tol)
            this$msiWidget2$PlotMassImageRGB()
          }
          this$MsiRedraw[[i]]$ReDraw <- F
        }
      }
      this$RedrawIsIdle <- F
    }
  }

  #Stop gtimer if widget is distroyed
  Widget_Disposed <- function (evt, ...)
  {
    #cat("Stopping spectraWidget draw timer\n")
    this$redrawTimer$stop_timer()
  }

  #A click on mass spectra widgets drives here. From here image recostruction will be called for various widgets
  SpectrumClicked <- function( channel, mass, tol )
  {
    this$MsiRedraw[[channel]]$mass <- mass
    this$MsiRedraw[[channel]]$tol <- tol
    this$MsiRedraw[[channel]]$ReDraw <- T
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  AddSpectra <- function ( mass_axis, intensity_list, color_list, id_list, calling_image_string, normalizations,  mz_min, mz_max )
  {
    for( i in 1:length(intensity_list))
    {
      this$spectraWidget$AddSpectra(mass_axis, intensity_list[[i]]/normalizations[i], col = color_list[[i]], name = paste(calling_image_string,"_ID",as.character(id_list[[i]]), sep = ""))
    }

    if( !is.null(mz_min) && !is.null(mz_max))
    {
      this$spectraWidget$SetPlottedMassRange(mz_min, mz_max)
    }
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  ClearSpectraPlot <- function(id_list, calling_image_string)
  {
    for( i in 1:length(id_list))
    {
      this$spectraWidget$RmSpectra(paste(calling_image_string,"_ID",id_list[i], sep = ""))
    }
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  GetPlotedSpectraInfo <- function( strImgName )
  {
    SpcData <-this$spectraWidget$GetSpectraInfo()
    myImgIds <- c ()
    myImgColors <- c()
    for( i in 1:nrow(SpcData))
    {
      name_id <- unlist(strsplit(as.character(SpcData[i, "names"]), "_ID"))
      if(length(name_id) == 2)
      {
        if( name_id[1] == strImgName )
        {
          myImgIds <- c(myImgIds, as.numeric(name_id[2]) )
          myImgColors <- c(myImgColors, as.character(SpcData[i, "colors"]))
        }
      }
    }

    #Get plotted mass range
    mz_range <- this$spectraWidget$GetPlottedMassRange()

    return(list(ID=myImgIds, color =  myImgColors, mz_min = mz_range$mz_min, mz_max = mz_range$mz_max))
  }

  #GUI builder
  window <- gWidgets2::gwindow ( "MSI Reconstruction" , visible = F )
  Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)
  Grp_Ims <- gWidgets2::gpanedgroup(horizontal = T, container = Grp_Top )

  msiWidget1 <- .MSImagePlotWidget(in_img = img1 , parent_widget = Grp_Ims, AddSpectra_function = this$AddSpectra, GetSpectraInfo_function = this$GetPlotedSpectraInfo, ClearSpectraPlot_function = this$ClearSpectraPlot, meanSpectrumColor = "red", widget_name = "imgLeft")
  if( !is.null(img2))
  {
    msiWidget2 <- .MSImagePlotWidget(in_img = img2 , parent_widget = Grp_Ims, AddSpectra_function = this$AddSpectra, GetSpectraInfo_function = this$GetPlotedSpectraInfo, ClearSpectraPlot_function = this$ClearSpectraPlot, meanSpectrumColor = "blue", widget_name = "imgRight")
  }
  spectraFrame<-gWidgets2::gframe("Spectra Viewer", container = Grp_Top,  fill = T, expand = F, spacing = 2 )
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F,  display_sel_red = T, display_sel_green = T, display_sel_blue = T, display_clearall_button = T, useInternalRedrawTimer = F)

  
  spectraWidget$AddSpectra(  img1$mass, img1$mean, col = "red", name = "imgLeft_ID0")
  if(!is.null(img2))
  {
    spectraWidget$AddSpectra(  img2$mass, img2$mean, col = "blue", name = "imgRight_ID0")
  }

  spectraWidget$ZoomResetClicked()

  #Plot a initial image which is the maximum peak in mean spectrum with a tolereance of 100 ppm of the mass range
  startTol <- 100/1e6*(max(img1$mass)-min(img1$mass))
  startMass = img1$mass[which.max(img1$mean)]
    
  for( i in 1:3)
  {
    SpectrumClicked( channel = i, mass = startMass, tol = startTol)
    spectraWidget$SetSelectedMassTol( channel = i, mass = startMass, tol = startTol )
  }

  #Start the redraw timer
  redrawTimer <- gWidgets2::gtimer(10, this$ReDrawByTimer)
  gWidgets2::addHandlerDestroy( obj = window, handler = this$Widget_Disposed ) #Connect to widget dispose to stop the draw timer
  
  #Set init windows size
  gWidgets2::visible(window) <- TRUE
  gWidgets2::size(window) <- c(1280, 960)
  gWidgets2::svalue(Grp_Top) <- 0.65 #65% of space allocated to MSI image

  ## Set the name for the class
  class(this) <- append(class(this),"MsiWindows")
  gc()

  #A loop to block exectuion
  if(lockExecution)
  {
    while( gWidgets2::isExtant(this$window ) )
    {
      Sys.sleep(0.1)
    }
  }
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
}


