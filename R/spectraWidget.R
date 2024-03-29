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

###A GUI to display spectra in an interactive way

#' Plot Mass Spectra in a interactive way.
#'
#' @param mass the mass axis of spectrum
#' @param intensity a intensity vector of spectrum
#' @param peaks_mass a vector of peak masses to be labeled on spectrum
#' @param peaks_intensity a vector of peak intensities to be labeled on spectrum
#' @param ref_mass a vector of reference masses to be represented as vertical dashed lines
#' @param col the color for the spectrum
#'
#' Simply call this function each time you want to plot a mass spectrum and it will be overlayed with the current spectra.
#' Peaks can be specified in the same call and it will be labeled over the spectrum.
#'
#' @export
#'
plotSpectra<-function( mass = NULL, intensity = NULL, peaks_mass = NULL, peaks_intensity = NULL, ref_mass = NULL, col = "")
{
  if( !exists( x = ".SpectraWidget", mode = "environment") )
  {
    #Spectra plot windows does not exists, create it
    options(guiToolkit="tcltk") #force to use tcltk
    oldWarning<-options()$warn
    options(warn = -1)
    
    window_spectra <- gWidgets2::gwindow ( "Spectra Plot" , visible = FALSE, width = 750, height = 440)
    .SpectraWidget<<-.SpectraPlotWidget( parent_widget = window_spectra )
    gWidgets2::addHandlerDestroy( obj = window_spectra, handler = .plotSpectraWindow_Disposed )

    #window_spectra$widget$present()
    gWidgets2::visible(window_spectra) <- T
    
    #Restore warnings level
    options(warn = oldWarning)
    rm(oldWarning)
  }

  if( !is.null(mass) && !is.null(intensity) )
  {
    .SpectraWidget$AddSpectra( mass_data = mass, intensity_data = intensity,
                             mass_peaks = peaks_mass, intensity_peaks = peaks_intensity, col = col)
  }

  if( !is.null(ref_mass) )
  {
    .SpectraWidget$SetRefMass( ref_mass  )
  }
}

.plotSpectraWindow_Disposed <- function (evt, ...)
{
  rm(.SpectraWidget, envir = .GlobalEnv)
  gc()
}

.SpectraPlotWidget <- function( parent_widget = gWidgets2::gwindow( "Default SpectraPlotWidget" , visible = FALSE ), top_window_widget = NULL,  clicFuntion = NULL, showOpenFileButton = T,
                                display_sel_red = F, display_sel_green = F, display_sel_blue = F, display_sel_spins = T, display_clearall_button = F, max_spectra_limit = 50, useInternalRedrawTimer = T)
{
  options(guiToolkit="tcltk") 
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members

  parent <- parent_widget
  rm(parent_widget)

  if(is.null(top_window_widget))
  {
    top_window <- parent
  }
  else
  {
    top_window <- top_window_widget
  }
  rm(top_window_widget)
  clicFun <- clicFuntion
  rm(clicFuntion)
  AnteriorClickMz <- 0
  plot_device <- NULL
  radio_buttons <- 0
  lbl_mz_coords <- 0
  lbl_in_coords <- 0
  LABEL_LENGTH <- 10
  MouseWheelFunction <- 0
  spectra_data <- list() #Data structure controlled by functions
  ref_mass <- NULL
  mz_lim <- c(0, 1)
  in_lim <- c(0, 1)
  data_mass_range <- c() #c(min(mass_data),max(mass_data)),
  SelIon_mz_R <- NA
  SelIon_tol_R <- NA
  SelIon_mz_G <- NA
  SelIon_tol_G <- NA
  SelIon_mz_B <- NA
  SelIon_tol_B <- NA
  CurrentSelTool <- "Zoom" #Stores the curren state of the sel tool, can be: Zoom, Red, Green and Blue
  MAX_SPECTRA_LIMIT <- max_spectra_limit #Maximum number of spectra that can be added
  ReDraw <- F #Signal when spectra must be redraw
  MAX_MASS_SEL_RANGE <- 200 #Max range of masses to allow selection (in data points)
  PLOT_MARGIN_LEFT <- 77 #Empiric value
  PLOT_MARGIN_RIGHT <- 10 #Empiric value
  ReDrawRedImg <- F
  ReDrawGreenImg <- F
  ReDrawBlueImg <- F

  #Stop gtimer if widget is distroyed
  Widget_Disposed <- function (evt, ...)
  {
    #cat("Stopping spectraWidget draw timer\n")
    this$redrawTimer$stop_timer()
  }

  #Create spectrum object
  #A spectrum is created with a defined mass and intensity. This two params are not modifiable.
  #The original color must be also specified, but this param can be changed after.
  #A spectrum can also provide peaks. This can be added/modified latter.
  CreateSpectrumObj <- function (mass, intensity, color, mass_peaks = NULL, intensity_peaks = NULL, active = T)
  {
    return(list(mass = mass, intensity = intensity, color = color, mass_peaks = mass_peaks, intensity_peaks = intensity_peaks, enabled = active))
  }

  #Set spectrum color
  SetSpectrumColors <- function(spcObj, color)
  {
    spcObj$color <- color
    return(spcObj)
  }

  #Set spectrum peaks
  SetSpectrumPeaks <- function(spcObj, mass_peaks, intensity_peaks )
  {
    spcObj$mass_peaks <- mass_peaks
    spcObj$intensity_peaks <- intensity_peaks
    return(spcObj)
  }

  #Set new color to the internal data list of spectra
  ChangeSpectrumColor <- function(name, color )
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]<-this$SetSpectrumColors(this$spectra_data[[as.character(name)]], color)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set new peaks to the internal data list of spectra
  ChangeSpectrumPeaks <- function(name, peaks_mass, peaks_intensity )
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]<-this$SetSpectrumPeaks(this$spectra_data[[as.character(name)]], peaks_mass, peaks_intensity)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set enabled state of a spectrum, if enabled then it is visble
  SetSpectrumEnabled <- function(name, enabled)
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]$enabled <- enabled
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Return the spectrum enabled state from a given name
  GetSpectrumEnabled <- function(name)
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    return(his$spectra_data[[as.character(name)]]$enabled)
  }

  #Return all names of the current spectra list
  GetSpectraInfo <- function()
  {
    spcNames <- names(this$spectra_data)
    spcColors <- c()
    for( i in 1:length(spcNames))
    {
      spcColors <- c(spcColors, this$spectra_data[[as.character(spcNames[i])]]$color)
    }
    return(data.frame( names = spcNames, colors = spcColors))
  }

  #Return the current plotted mass range
  GetPlottedMassRange <- function()
  {
    return(list( mz_min = this$mz_lim[1], mz_max = this$mz_lim[2]))
  }

  #Set a mew plotted mass range
  SetPlottedMassRange <- function(mz_min, mz_max)
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim[1] <- mz_min
    this$mz_lim[2] <- mz_max
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Add spectrum data
  AddSpectra <- function( mass_data, intensity_data, mass_peaks = NULL, intensity_peaks = NULL, col = "", name = "", add_enabled = T)
  {
    if( length(this$spectra_data) >= this$MAX_SPECTRA_LIMIT)
    {
      gWidgets2::gmessage(paste("The limit number of spectra (",this$MAX_SPECTRA_LIMIT , ") has been reached. Remove some spectrum to add a new one.", sep =""), icon = "error")
      return()
    }

    if(col == "") #Set a default color if color is not provided
    {
      col<-as.character(sample(rainbow(100), 1))
    }
    if(name == "" ) #Provide a unique name if name is not specified
    {
      while(name == "")
      {
        name <- as.character(sample(1:1e3, 1))
        if( length(which( names(this$spectra_data) == name)))
        {
          name <- ""
        }
      }
    }
    this$spectra_data[[as.character(name)]]<-this$CreateSpectrumObj(mass_data, intensity_data, col, mass_peaks, intensity_peaks, add_enabled)
    this$data_mass_range <- c(min(sapply(this$spectra_data, function(x){ return( x$mass[1] ) })),max(sapply(this$spectra_data, function(x){ return( x$mass[length(x$mass)] ) })))

    #Set Spin_massSel range properly
    if(!is.null(this$Spin_massSel))
    {
      gWidgets2::blockHandlers(this$Spin_massSel)
      tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), from = this$data_mass_range[1])
      tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), to = this$data_mass_range[2])
      gWidgets2::unblockHandlers(this$Spin_massSel)
    }

    if(length(this$spectra_data) == 1) #First spectrum added, so zoom properly
    {
      this$mz_lim <- this$data_mass_range
    }
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set ref mass data, ref masses will be ploted as vertical dashed lines
  SetRefMass <- function( mass_data )
  {
    this$ref_mass <- mass_data
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Remove spectrum data
  RmSpectra <- function( name )
  {
    this$spectra_data<- this$spectra_data[-(which(names(this$spectra_data) == as.character(name)))]

    #Recompute mz_limits, data_mass range to scale acording the removed spectra
    if(length(this$spectra_data) > 0)
    {
      this$mz_lim <- c(min(sapply(this$spectra_data, function(x){ return( x$mass[1] ) })),max(sapply(this$spectra_data, function(x){ return( x$mass[length(x$mass)] ) })))
      this$in_lim <- c(0 ,max(sapply(this$spectra_data, function(x){ return( max(x$intensity)) }))*1.1)
      this$data_mass_range <- this$mz_lim

      #Set Spin_massSel range properly
      if(!is.null(this$Spin_massSel))
      {
        gWidgets2::blockHandlers(this$Spin_massSel)
        tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), from = this$data_mass_range[1])
        tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), to = this$data_mass_range[2])
        gWidgets2::unblockHandlers(this$Spin_massSel)
      }
    }
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Clear all spectrum data
  ClearSpectra <- function( )
  {
    this$spectra_data<-list()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Redraw MS image on parent widget in a given channel
  ReDrawParentMSI <- function()
  {
    if(this$ReDrawRedImg) #Red
    {
      this$clicFun(1, this$SelIon_mz_R, this$SelIon_tol_R)
    }
    if(this$ReDrawGreenImg) #Green
    {
      this$clicFun(2, this$SelIon_mz_G, this$SelIon_tol_G)
    }
    if(this$ReDrawBlueImg) #Blue
    {
      this$clicFun(3, this$SelIon_mz_B, this$SelIon_tol_B)
    }

    mz_sel_spin <- NULL
    mz_tol_spin <- NULL
    if(this$CurrentSelTool == "Red" )
    {
      mz_sel_spin <- this$SelIon_mz_R
      mz_tol_spin <- this$SelIon_tol_R
    }

    if(this$CurrentSelTool == "Green" )
    {
      mz_sel_spin <- this$SelIon_mz_G
      mz_tol_spin <- this$SelIon_tol_G
    }

    if(this$CurrentSelTool == "Blue" )
    {
      mz_sel_spin <- this$SelIon_mz_B
      mz_tol_spin <- this$SelIon_tol_B
    }

    #Set Spin_massSel range properly
    if( !is.null(mz_sel_spin) && !is.null(mz_tol_spin) && !is.null(this$Spin_massSel) && !is.null(this$Spin_TolSel))
    {
      gWidgets2::blockHandlers(this$Spin_massSel)
      gWidgets2::blockHandlers(this$Spin_TolSel)
      gWidgets2::svalue(this$Spin_massSel) <- mz_sel_spin
      gWidgets2::svalue(this$Spin_TolSel) <- mz_tol_spin
      tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), increment = mz_tol_spin)
      gWidgets2::unblockHandlers(this$Spin_massSel)
      gWidgets2::unblockHandlers(this$Spin_TolSel)
    }

    this$ReDrawRedImg <- F
    this$ReDrawGreenImg <- F
    this$ReDrawBlueImg <- F
  }

  #Redraw ggraph with interpolation using a timer
  ReDrawByTimer <- function( data )
  {
    if( this$ReDraw )
    {
      redrawPlotWidget(this$plot_device)
    }
    this$ReDraw <- F
  }
  
  #Redraw
  ReDrawPlot <- function()
  {
    #Redraw parent MS image if necessari
    if( !is.null( this$clicFun ))
    {
      ReDrawParentMSI()
    }
    
    #Reduce spectra data size to speed-up ploting
    mass_range <- this$mz_lim
    if(is.null(this$plot_device))
    {
      npoints <- 10000 #in init condition the widget is still not avaialble so I cannot get its size! So I use a default size of 10000
    }
    else
    {
      npoints <- 10*as.numeric(tcltk::tkwinfo("width", this$plot_device)) 
    }
    
    #Visible before plot() forces the target divice for ploting
    par(mar = c(3.1, 5.1, 0.5, 0.5), cex = 0.7, xaxs = "i", yaxs = "i")
    
    #Init Plot
    #i_axt<-pretty(this$in_lim[1]:this$in_lim[2], n = 5)
    i_axt<-seq(from = this$in_lim[1], to = this$in_lim[2], length.out = 5)
    plot(x=0, xlim = this$mz_lim, ylim = this$in_lim, type = "n", xlab = "", ylab ="", yaxt ="n")
    axis(2, at = i_axt, labels = sprintf("%.2e",i_axt), las = 1)
    
    #Draw ref masses as vertical lines
    if(!is.null(this$ref_mass))
    {
      abline(v = this$ref_mass, col = "grey", lty = 2)
    }
    
    #Plot selection range Red
    if(!is.na(this$SelIon_mz_R) && !is.na(this$SelIon_tol_R))
    {
      rect(xleft = this$SelIon_mz_R - this$SelIon_tol_R, xright = this$SelIon_mz_R + this$SelIon_tol_R, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightsalmon", border = "red3")
    }
    
    #Plot selection range Green
    if(!is.na(this$SelIon_mz_G) && !is.na(this$SelIon_tol_G))
    {
      rect(xleft = this$SelIon_mz_G - this$SelIon_tol_G, xright = this$SelIon_mz_G + this$SelIon_tol_G, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightgreen", border = "green3")
    }
    
    #Plot selection range Blue
    if(!is.na(this$SelIon_mz_B) && !is.na(this$SelIon_tol_B))
    {
      rect(xleft = this$SelIon_mz_B - this$SelIon_tol_B, xright = this$SelIon_mz_B + this$SelIon_tol_B, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightblue", border = "blue3")
    }
    
    if(length(this$spectra_data) > 0)
    {
      for(li in 1:length(this$spectra_data))
      {
        if( this$spectra_data[[li]]$enabled)
        {
          
          #Plot the spectrum
          plotData <- ReduceDataPointsC(this$spectra_data[[li]]$mass, this$spectra_data[[li]]$intensity, this$mz_lim[1], this$mz_lim[2], npoints)
          lines(x = plotData$mass, y = plotData$intensity, col= this$spectra_data[[li]]$color)
          if(length(plotData$mass) <= 200)
          {
            points(x = plotData$mass, y = plotData$intensity, col= this$spectra_data[[li]]$color, pch = 20)
          }
          
          #Plot labels
          if( !is.null(this$spectra_data[[li]]$mass_peaks) && !is.null(this$spectra_data[[li]]$intensity_peaks) )
          {
            #cat(paste("Mz peaks:", this$spectra_data[[li]]$mass_peaks, "int peaks:", this$spectra_data[[li]]$intensity_peaks, "\n"))
            pk_lbl <- sprintf("%.4f", this$spectra_data[[li]]$mass_peaks) #4 decimals
            text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = pk_lbl, pos = 3, cex = 0.8, offset = 1.1)
            un_lbl<-sapply(pk_lbl, function(x) { paste(rep("_", nchar(x)), collapse = "") })
            text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = un_lbl, pos = 3, cex = 0.8, offset = 1.1)
            text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = rep("|", length(pk_lbl)), pos = 3, cex = 0.8)
          }
        }
      }
    }
    this$SetStateAccordingSelTool()
  }

  #OpenTXT
  OpenTXT <- function( evt, ... )
  {
    fname<-file.choose()
    spect<-read.table(fname, header = F, sep = " ")
    this$AddSpectra( spect[,1], spect[,2], name = basename(fname))
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Reset Zoom button clicked
  ZoomResetClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim <- this$data_mass_range
    this$in_lim <- c(0 ,max(sapply(this$spectra_data, function(x){ return( max(x$intensity)) }))*1.1)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Auto Zoomin Mz axis
  ZoomMzClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim <- this$data_mass_range
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Auto Zomming Intensity Axis
  ZoomInClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Intensity auto-zoom
  AutoZoomIntensity <- function( )
  {
    if(length(this$spectra_data) == 0) return()
    this$in_lim <- c(0 ,max(unlist(lapply(this$spectra_data, function(x){ max(x$intensity[ which( x$mass >= this$mz_lim[1] & x$mass <= this$mz_lim[2], arr.ind = T ) ]) })))*1.1)
  }

  #Mz zoom in a defined range with autoscaling intensity
  ZoomMzRange <- function(mzLow, mzHigh)
  {
    #Mz zooming
    mzLow<-max(mzLow, this$data_mass_range[1])
    mzHigh<-min(mzHigh, this$data_mass_range[2])
    this$mz_lim<- c(mzLow, mzHigh)

    #Intensity zoom
    this$AutoZoomIntensity()

    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Grab Mouse Selection Changes on plot
  OnSelection <- function( x ,y )
  {
    if(length(this$spectra_data) == 0) return()
    if(this$CurrentSelTool == "Zoom")
    {
      if(abs(x[1] - x[2]) > 0.1)
      {
        #Mz zooming
        this$ZoomMzRange(min(x), max(x))
      }
    }
    else
    {
      top_left <- min(x)
      top_right <- max(x)
      mz_tol <-(top_right - top_left)/2
      mz_sel <- top_left + mz_tol
      mz_tol<-max(mz_tol, 0)

      #Use the tolerance spin if the spectrum was just clicked
      if(mz_tol == 0 && !is.null(this$Spin_TolSel))
      {
        mz_tol <- gWidgets2::svalue(this$Spin_TolSel)
      }

      #Limit selection to 5 Da to avoid selecting large parts of spectra and filling RAM
      if( length(this$spectra_data) == 0)
      {
        cat("No spectra to select\n")
        return()
      }
      
      dpSelL <- which.min(abs(this$spectra_data[[1]]$mass - (mz_sel - mz_tol)))
      dpSelR <- which.min(abs(this$spectra_data[[1]]$mass - (mz_sel + mz_tol)))
      if(( dpSelR - dpSelL ) > this$MAX_MASS_SEL_RANGE )
      {
        cat(paste("Ion selection in a range of", mz_tol*2, "Da has been aborted. To large data sector.\n"))
        return()
      }

      if(this$CurrentSelTool == "Red")
      {
        this$SelIon_mz_R <- mz_sel
        this$SelIon_tol_R <- mz_tol
        this$ReDrawRedImg <- T
      }
      else if ( this$CurrentSelTool == "Green")
      {
        this$SelIon_mz_G <- mz_sel
        this$SelIon_tol_G <- mz_tol
        this$ReDrawGreenImg <- T
      }
      else if ( this$CurrentSelTool == "Blue")
      {
        this$SelIon_mz_B <- mz_sel
        this$SelIon_tol_B <- mz_tol
        this$ReDrawBlueImg <- T
      }
      #Plot selection in spectra
      this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
    }
  }

  #Grab mouse cursor on plot
  OnMouseMotion <- function( x, y )
  {
    if(length(this$spectra_data) == 0) return()

    #Update Labels
    if(x >= this$mz_lim[1] &&
       x <= this$mz_lim[2] &&
       y >= this$in_lim[1] &&
       y <= this$in_lim[2])
    {
      mz_txt<-sprintf(paste("%-",this$LABEL_LENGTH, ".4f", sep = ""), x)
      in_txt<-sprintf(paste("%-",this$LABEL_LENGTH, ".2e", sep = ""), y)
    }
    else
    {
      in_txt<-mz_txt<-paste(rep(" ", this$LABEL_LENGTH), collapse = "")
    }
    this$lbl_mz_coords$set_value(mz_txt)
    this$lbl_in_coords$set_value(in_txt)
  }

  #Zoom on spectra plot handler
  ScrollEventOnSpectra <- function( direction, x, y )
  {
    
    if(length(this$spectra_data) == 0) return()

    dir<- as.double(direction)
    mz_min<-this$data_mass_range[1]
    mz_max<-this$data_mass_range[2]

    if(this$MouseWheelFunction == 0)
    {
      #Mz scrolling
      range<-abs(this$mz_lim[2] - this$mz_lim[1])
      mz_lim<-this$mz_lim - dir*range*0.05
      range<-abs(mz_lim[2] - mz_lim[1])
      if(mz_lim[1] < mz_min)
      {
        #Clip to min
        mz_lim<-c(mz_min, mz_min + range)
      }
      if(mz_lim[2] > mz_max)
      {
        #Clip to max
        mz_lim<-c(mz_max - range, mz_max)
      }
      if(range > abs(mz_max - mz_min))
      {
        #out of range! clip in full spectra
        mz_lim<-c(mz_min, mz_max)
      }
      this$mz_lim<- mz_lim
    }

    else if(this$MouseWheelFunction == 1)
    {
      #Mz zooming
      top_left <- x -  abs(0.1 + dir)*(x - this$mz_lim[1])
      top_right <- x + abs(0.1 + dir)*(this$mz_lim[2] - x)
      top_left<-max(top_left, mz_min)
      top_right<-min(top_right, mz_max)
      this$mz_lim<- c(top_left, top_right)
    }
    else
    {
      #Intensity scaling
      range<-abs(this$in_lim[2] - this$in_lim[1])
      range<-range - dir*range*0.05
      this$in_lim<- c(0, range)
    }

    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Key Press handler
  OnKeyPress <- function( K )
  {
    if( K == "Control_L" || K == "Control_R" )
    {
      if(this$MouseWheelFunction == 0)
      {
        this$MouseWheelFunction<-1
      }
    }
    else if( K == "Shift_L" || K == "Shift_R" )
    {
      if(this$MouseWheelFunction == 0)
      {
        this$MouseWheelFunction<-2
      }
    }
    else if( K == "z" || K == "Z" )
    {
      #z key press
      if(!is.null(this$Btn_ZoomTool))
      {
        tcltk::tkselect(this$Btn_ZoomTool)
        this$ZoomToolSel()
      }
    }
    else if( K == "r" || K == "R" )
    {
      #r key press
      if(!is.null(this$Btn_SelRedTool))
      {
        tcltk::tkselect(this$Btn_SelRedTool)
        this$RedToolSel()
      }
    }
    else if( K == "g" || K == "G" )
    {
      #g key press
      if(!is.null(this$Btn_SelGreenTool))
      {
        tcltk::tkselect(this$Btn_SelGreenTool)
        this$GreenToolSel()
      }
    }
    else if( K == "b" || K == "B" )
    {
      #b key press
      if(!is.null(this$Btn_SelBlueTool))
      {
        tcltk::tkselect(this$Btn_SelBlueTool)
        this$BlueToolSel()
      }
    }
  }

  #Key Release handler
  OnKeyRelease <- function( K )
  {
    if( K == "Control_L" || K == "Control_R" )
    {
      if(this$MouseWheelFunction == 1)
      {
        this$MouseWheelFunction<-0
      }
    }
    else if( K == "Shift_L" || K == "Shift_R" )
    {
      if(this$MouseWheelFunction == 2)
      {
        this$MouseWheelFunction<-0
      }
    }
  }

  #Windows lost focuts, used to restore zoom status
  OnLostFocus <- function( )
  {
    this$MouseWheelFunction<-0
  }

  #Clear all spectra button click
  ClearSpectraClicked <- function( evt, ... )
  {
    this$ClearSpectra()
  }

  #Zoom tool has been selected
  ZoomToolSel <- function( ... )
  {
    if(is.null(this$Btn_ZoomTool))
    {
     return()  
    }
    
    if( getValue_coloredCheckBox(this$Btn_ZoomTool))
    {
      this$CurrentSelTool <- "Zoom"
      if(!is.null(this$Btn_SelRedTool))
      {
        tcltk::tkdeselect(this$Btn_SelRedTool)
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        tcltk::tkdeselect(this$Btn_SelGreenTool)
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        tcltk::tkdeselect(this$Btn_SelBlueTool)
      }

      #Set ion manually selection visibility
      #TODO there is no hide/show implementation for tcltk widgets, so I workarround it by disabling them.. but it is not very elegant
      if(!is.null(this$Spin_massSel))
      {
        gWidgets2::enabled(Lbl_massSel) <- F
        gWidgets2::enabled(Spin_massSel) <- F
        gWidgets2::enabled(Lbl_TolSel) <- F
        gWidgets2::enabled(Spin_TolSel) <- F
      }
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelGreenTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        tcltk::tkselect(this$Btn_ZoomTool)
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Set mass and tolerance spins value and visibility
  SetMassTolSping <- function( mzValue, tolValue, bVisible = T)
  {
    if(is.null(this$Spin_massSel) || is.null(this$Spin_TolSel))
    {
      return()
    }

    #Set Spin_massSel range properly
    gWidgets2::blockHandlers(this$Spin_massSel)
    gWidgets2::blockHandlers(this$Spin_TolSel)
    gWidgets2::svalue(this$Spin_massSel) <- mzValue
    gWidgets2::svalue(this$Spin_TolSel) <- tolValue
    tcltk::tkconfigure(gWidgets2::getToolkitWidget(this$Spin_massSel), increment = tolValue) 
    gWidgets2::unblockHandlers(this$Spin_massSel)
    gWidgets2::unblockHandlers(this$Spin_TolSel)

    #Set ion manually selection visibility
    gWidgets2::enabled(Lbl_massSel) <- bVisible
    gWidgets2::enabled(Spin_massSel) <- bVisible
    gWidgets2::enabled(Lbl_TolSel) <- bVisible
    gWidgets2::enabled(Spin_TolSel) <- bVisible
  }

  #Red tool selected
  RedToolSel <- function( ... )
  {
    if( getValue_coloredCheckBox(this$Btn_SelRedTool))
    {
      this$CurrentSelTool <- "Red"
      if(!is.null(this$Btn_ZoomTool))
      {
        tcltk::tkdeselect(this$Btn_ZoomTool)
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        tcltk::tkdeselect(this$Btn_SelGreenTool)
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        tcltk::tkdeselect(this$Btn_SelBlueTool)
      }

      this$SetMassTolSping(  this$SelIon_mz_R, this$SelIon_tol_R, T)
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- getValue_coloredCheckBox(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelGreenTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        tcltk::tkselect(this$Btn_SelRedTool)
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Green tool selected
  GreenToolSel <- function( ... )
  {
    if( getValue_coloredCheckBox(this$Btn_SelGreenTool))
    {
      this$CurrentSelTool <- "Green"
      if(!is.null(this$Btn_ZoomTool))
      {
        tcltk::tkdeselect(this$Btn_ZoomTool)
      }
      if(!is.null(this$Btn_SelRedTool))
      {
        tcltk::tkdeselect(this$Btn_SelRedTool)
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        tcltk::tkdeselect(this$Btn_SelBlueTool)
      }

      this$SetMassTolSping(  this$SelIon_mz_G, this$SelIon_tol_G, T)
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- getValue_coloredCheckBox(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        tcltk::tkselect(this$Btn_SelGreenTool)
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Blue tool selected
  BlueToolSel <- function( ... )
  {
    if( getValue_coloredCheckBox(this$Btn_SelBlueTool))
    {
      this$CurrentSelTool <- "Blue"
      if(!is.null(this$Btn_ZoomTool))
      {
        tcltk::tkdeselect(this$Btn_ZoomTool)
      }
      if(!is.null(this$Btn_SelRedTool))
      {
        tcltk::tkdeselect(this$Btn_SelRedTool)
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        tcltk::tkdeselect(this$Btn_SelGreenTool)
      }

      this$SetMassTolSping( this$SelIon_mz_B, this$SelIon_tol_B, T)
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- getValue_coloredCheckBox(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- getValue_coloredCheckBox(this$Btn_SelGreenTool) | bTest
      }
      if( bTest == F)
      {
        tcltk::tkselect(this$Btn_SelBlueTool)
      }
    }
    this$SetStateAccordingSelTool()
  }

  SetStateAccordingSelTool <- function ()
  {
    if(!is.null(this$plot_device))
    {
      if (this$CurrentSelTool == "Zoom")
      {
        tcltk::tkconfigure(this$plot_device, cursor = "target")
      }
      else
      {
        tcltk::tkconfigure(this$plot_device, cursor = "based_arrow_down")
      }
    }
  }

  #Set the tool to use externally, can be: Zoom, Red, Green or Blue
  SetActiveTool <- function (tool)
  {
    #Start by deselecting all since tcltk do not trigger signals when selecting programatically
    if(!is.null(this$Btn_ZoomTool)) { tcltk::tkdeselect(this$Btn_ZoomTool)}
    if(!is.null(this$Btn_SelRedTool)) { tcltk::tkdeselect(this$Btn_SelRedTool)}
    if(!is.null(this$Btn_SelGreenTool)) { tcltk::tkdeselect(this$Btn_SelGreenTool)}
    if(!is.null(this$Btn_SelBlueTool)) { tcltk::tkdeselect(this$Btn_SelBlueTool)}
    
    if( tool == "Zoom")
    {
      if(!is.null(this$Btn_ZoomTool))
      {
        tcltk::tkselect(this$Btn_ZoomTool)
        this$ZoomToolSel()
      }
    }
    else if( tool == "Red")
    {
      if(!is.null(this$Btn_SelRedTool))
      {
        tcltk::tkselect(this$Btn_SelRedTool)
        this$RedToolSel()
      }
    }
    else if( tool == "Green")
    {
      if(!is.null(this$Btn_SelGreenTool))
      {
        tcltk::tkselect(this$Btn_SelGreenTool)
        this$GreenToolSel()
      }
    }
    else if( tool == "Blue")
    {
      if(!is.null(this$Btn_SelBlueTool))
      {
        tcltk::tkselect(this$Btn_SelBlueTool)
        this$BlueToolSel()
      }
    }
  }

  MassSelSpinChanged <- function(...)
  {
    if(this$CurrentSelTool == "Red")
    {
      this$SelIon_mz_R <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawRedImg <- T
    }
    else if ( this$CurrentSelTool == "Green")
    {
      this$SelIon_mz_G <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawGreenImg <- T
    }
    else if ( this$CurrentSelTool == "Blue")
    {
      this$SelIon_mz_B <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawBlueImg <- T
    }

    #Set zoom range to see the selected ion centered
    currentRange <- this$mz_lim[2] - this$mz_lim[1]
    this$ZoomMzRange( gWidgets2::svalue(this$Spin_massSel) - (currentRange/2), gWidgets2::svalue(this$Spin_massSel) + (currentRange/2))

    #Plot selection in spectra
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  TolSelSpinChanged <- function(...)
  {
    if(this$CurrentSelTool == "Red")
    {
      this$SelIon_tol_R <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawRedImg <- T
    }
    else if ( this$CurrentSelTool == "Green")
    {
      this$SelIon_tol_G <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawGreenImg <- T
    }
    else if ( this$CurrentSelTool == "Blue")
    {
      this$SelIon_tol_B <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawBlueImg <- T
    }
    #Plot selection in spectra
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  SetSelectedMassTol <- function (channel, mass, tol)
  {
    if(channel == 1)
    {
      this$SelIon_mz_R <- mass
      this$SelIon_tol_R <- tol
    }
    if(channel == 2)
    {
      this$SelIon_mz_G <- mass
      this$SelIon_tol_G <- tol
    }
    if(channel == 3)
    {
      this$SelIon_mz_B <- mass
      this$SelIon_tol_B <- tol
    }
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

    #Build GUI
  Grp_Top <- gWidgets2::ggroup(horizontal = F, container = this$parent, fill = T, expand = T)
  Grp_Buttons<- gWidgets2::ggroup(horizontal = T, container = Grp_Top, fill = F, expand = F)
  Btn_reset_zoom<- gWidgets2::gbutton(text = "Reset Zoom", handler = this$ZoomResetClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_mz<- gWidgets2::gbutton(text = "Auto m/z", handler = this$ZoomMzClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_in<- gWidgets2::gbutton(text = "Auto Intensity", handler = this$ZoomInClicked, action = this, container = Grp_Buttons)
  if(display_clearall_button)
  {
    Btn_RemoveSpectra<- gWidgets2::gbutton(text = "Clear all", handler = this$ClearSpectraClicked, action = this, container = Grp_Buttons)
  }
  
  if(display_sel_red | display_sel_green | display_sel_blue )
  {
    Btn_ZoomTool <- coloredCheckBox(text  = "Zoom", checked = T, handler =  this$ZoomToolSel, container = Grp_Buttons, bold = T)  
  }
  else
  {
    Btn_ZoomTool <- NULL
  }
  
  if( display_sel_red )
  {
    Btn_SelRedTool <- coloredCheckBox(text  = "Sel.Red", checked = F, handler =  this$RedToolSel, container = Grp_Buttons, foreground = "red", bold = T)
  }
  else
  {
    Btn_SelRedTool <- NULL
  }

  if( display_sel_green )
  {
    Btn_SelGreenTool <- coloredCheckBox(text  = "Sel.Green", checked = F, handler =  this$GreenToolSel, container = Grp_Buttons, foreground = "darkgreen", bold = T)
  }
  else
  {
    Btn_SelGreenTool <- NULL
  }

  if( display_sel_blue )
  {
    Btn_SelBlueTool <- coloredCheckBox(text  = "Sel.Blue", checked = F, handler =  this$BlueToolSel, container = Grp_Buttons, foreground = "darkblue", bold = T)
  }
  else
  {
    Btn_SelBlueTool <- NULL
  }

  #Dislay mass and tolerance spin boxes
  if(display_sel_spins)
  {
    Lbl_massSel <- gWidgets2::glabel("m/z:", container = Grp_Buttons)
    Spin_massSel <- gWidgets2::gspinbutton( from = 0, to = 1, value = 0, digits = 4, by = 0.1, container =  Grp_Buttons, handler = this$MassSelSpinChanged)
    Lbl_TolSel <- gWidgets2::glabel("+/-", container = Grp_Buttons)
    Spin_TolSel <- gWidgets2::gspinbutton( from = 0, to = 500, value = 0.1, digits = 4, by = 0.01, container =  Grp_Buttons, handler = this$TolSelSpinChanged)
    gWidgets2::enabled(Lbl_massSel) <- F
    gWidgets2::enabled(Spin_massSel) <- F
    gWidgets2::enabled(Lbl_TolSel) <- F
    gWidgets2::enabled(Spin_TolSel) <- F
  }

  gWidgets2::addSpring(Grp_Buttons)
  if(showOpenFileButton)
  {
    Btn_file_open<-gWidgets2::gbutton(text = "Open spectra TXT", handler = this$OpenTXT, action = this, container = Grp_Buttons)
  }
  rm(showOpenFileButton)

  plot_device <-createPlotWidget(parent = Grp_Top, 
                                      redraw_function = this$ReDrawPlot,
                                      MouseSelection_callback = this$OnSelection,
                                      MouseWheel_callback = this$ScrollEventOnSpectra, 
                                      MouseHover_callback = this$OnMouseMotion, 
                                      initial_width = 800, 
                                      initial_height = 170 ) 

  Grp_BottomLabel<-gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  lbl_help_info<-gWidgets2::glabel(" Mouse wheel -> m/z scroll\n Ctrl + Mouse wheel -> m/z zooming\n Shift + Mouse Wheel -> intensity scaling", container = Grp_BottomLabel)
  gWidgets2::addSpring(Grp_BottomLabel)
  lbl_mz<-gWidgets2::glabel("m/z:", container = Grp_BottomLabel)
  lbl_mz_coords<-gWidgets2::glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)
  lbl_int<-gWidgets2::glabel("Intensity:", container = Grp_BottomLabel)
  lbl_in_coords<-gWidgets2::glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)
  gWidgets2::font(lbl_mz_coords)<-list(family = "monospace", weight = "light", size = 8)
  gWidgets2::font(lbl_in_coords)<-list(family = "monospace", weight = "light", size = 8)
  gWidgets2::font(lbl_mz)<-list(family = "monospace", weight = "light", size = 8)
  gWidgets2::font(lbl_int)<-list(family = "monospace", weight = "light", size = 8)
  gWidgets2::font(lbl_help_info)<-list(family = "monospace", weight = "light", size = 8)

  #TclTk specific signal handlers
  tcltk::tkbind("all", "<KeyPress>", this$OnKeyPress )
  tcltk::tkbind("all", "<KeyRelease>", this$OnKeyRelease )
  tcltk::tkbind(plot_device, "<Leave>", this$OnLostFocus )
  
  #Start the redraw timer
  if(useInternalRedrawTimer)
  {
    redrawTimer <- gWidgets2::gtimer(10, this$ReDrawByTimer)
    gWidgets2::addHandlerDestroy( obj = this$top_window, handler = this$Widget_Disposed ) #Connect to widget dispose to stop the draw timer
  }

  ## Set the name for the class
  class(this) <- append(class(this),"SpectraPlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
