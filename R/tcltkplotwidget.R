#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2022 Pere Rafols Soler
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



#TclTk doc https://www.tcl.tk/man/tcl/TkCmd/contents.html
#' createPlotWidget.
#' 
#' creates a plotting widget to use with gWidgets2tcltk and allows direct mouse points selection.
#'
#' @param parent the gwidgets parent.
#' @param redraw_function a function to call to generate the plot.
#' @param MouseSelection_callback a callback function for the mouse click events.
#' @param MouseWheel_callback a callback function for the mouse wheel events. 
#' @param MouseHover_callback a callback function for the mouse move events.
#' @param initial_width inital plot width in pixels.
#' @param initial_height inital plot height in pixels.
#'
#' @return a plotting widget
#'
createPlotWidget <- function(parent, 
                             redraw_function, 
                             MouseSelection_callback = NULL,
                             MouseWheel_callback = NULL,
                             MouseHover_callback = NULL,
                             initial_width = 480,
                             initial_height = 480
)
{
  img_file <- tempfile(pattern = "rplot", fileext = ".png" )
  
  png(filename =  img_file, width =  initial_width, height = initial_height, res = 100, pointsize = 10, type = "cairo")
  par(xaxs = "i", yaxs = "i") #Force axes range to xlim and ylim provided, this should enable a better mouse clicks control
  try(redraw_function())
  
  #Getting data space coords ranges
  data_range <- par("usr")
  
  #Getting ploting are percentages
  plot_range <- par("plt")
  
  dev.off()
  
  
  #Prepare the image
  image <- tcltk::tkimage.create("photo", file = img_file) 
  
  tkwinviewplot <- tcltk::tkcanvas(gWidgets2::getToolkitWidget(parent),
                                   borderwidth = 0, 
                                   background = "black",
                                   width = as.integer(tcltk::tkimage.width(image)), 
                                   height = as.integer(tcltk::tkimage.height(image)))
  
  #Add a member to store data ranges and map info
  tkwinviewplot$plotting_data <- environment()
  tkwinviewplot$plotting_data$image <- image
  
  tcltk::tkpack(tkwinviewplot)
  tkwinviewplot$plotting_data$selRect <- NULL #No selection rectangle at begining

  #Add the image
  tkwinviewplot$plotting_data$canvasImg <- tcltk::tkcreate(tkwinviewplot, "image", 1, 1, image = tkwinviewplot$plotting_data$image, anchor = "nw" )

  tcltk::tkbind(tkwinviewplot, "<Destroy>", function() tcltk::.Tcl(paste("image delete", tkwinviewplot$plotting_data$image)))
  
  tkwinviewplot$file <- img_file
  tkwinviewplot$fun <- redraw_function
  tkwinviewplot$plotting_data$leftMousePressed <- F #Signal when mouse button is pressed
  
  #Set fill and expand
  tcltk::tkpack.configure(tkwinviewplot, "-fill", "both")
  tcltk::tkpack.configure(tkwinviewplot, "-expand", "1")
  
  setPlotWidgetDataRanges(tkwinviewplot, data_range, plot_range,  initial_width, initial_height)
  
  #Common method to test when an event happened in the plot area and execute a handler
  #Returns NULL if event happened outside plotting area otherwise, returns a vector with the following  format:
  # x_plot, y_plot, x_event, y_event
  #Where plot coords are relative to plotting area and event coords a related to tcltk widget
  testEventinPlottingArea <- function(event_x, event_y)
  {
    x <- as.integer(event_x)
    y <- as.integer(event_y)
    plotCoords <- mapDeviceCoords2PlotCoords(x, y, tkwinviewplot$plotting_data)
    
    ymin <- min(tkwinviewplot$plotting_data$data.ymin, tkwinviewplot$plotting_data$data.ymax)
    ymax <- max(tkwinviewplot$plotting_data$data.ymin, tkwinviewplot$plotting_data$data.ymax)
    
    if( plotCoords[1] >= tkwinviewplot$plotting_data$data.xmin &&
        plotCoords[1] <= tkwinviewplot$plotting_data$data.xmax &&
        plotCoords[2] >= ymin &&
        plotCoords[2] <= ymax )
    {
      return(c(plotCoords, as.integer(event_x), as.integer(event_y)))
    }
    else
    {
      return(NULL)
    }
  }
  
  #Connect callbacks
  leftMousePress <- function(x,y)
  {
    tkwinviewplot$plotting_data$last_press_coords <- testEventinPlottingArea(x, y)
    tkwinviewplot$plotting_data$leftMousePressed <- T
  }
  
  leftMouseRelease <- function(x,y)
  {
    tkwinviewplot$plotting_data$leftMousePressed <- F
    
    #Clear the selection rectangle
    if(!is.null(tkwinviewplot$plotting_data$selRect))
    {
      tcltk::tkdelete(tkwinviewplot, tkwinviewplot$plotting_data$selRect)
      tkwinviewplot$plotting_data$selRect <- NULL
    }
    
    
    if(!is.null(tkwinviewplot$plotting_data$last_press_coords))
    {
      release_coords <- testEventinPlottingArea(x, y)
      
      if(!is.null(release_coords))
      {
        MouseSelection_callback( x = c(tkwinviewplot$plotting_data$last_press_coords[1], release_coords[1]), 
                                 y = c(tkwinviewplot$plotting_data$last_press_coords[2], release_coords[2]))
      }
    }
  }
  
  if(!is.null(MouseSelection_callback))
  {
    tcltk::tkbind(tkwinviewplot, "<ButtonPress-1>", leftMousePress )
    tcltk::tkbind(tkwinviewplot, "<ButtonRelease-1>", leftMouseRelease )
  }
  
  scrollup_callback <- function(x, y)
  {
    pointer_coords <- testEventinPlottingArea(x, y)
    if(!is.null(pointer_coords))
    {
      MouseWheel_callback(direction = 1, x = pointer_coords[1],  y = pointer_coords[2])
    }
  }
  
  scrolldown_callback <- function(x, y)
  {
    pointer_coords <- testEventinPlottingArea(x, y) 
    if(!is.null(pointer_coords))
    {
      MouseWheel_callback(direction = -1, x = pointer_coords[1],  y = pointer_coords[2])
    }
  }
  
  scroll_windows <- function(D, x, y)
  {
    pointer_coords <- testEventinPlottingArea(x, y) 
    if(!is.null(pointer_coords))
    {
      MouseWheel_callback(direction = sign(as.integer(D)), x = pointer_coords[1],  y = pointer_coords[2])
    }
  }
  
  if(!is.null(MouseWheel_callback))
  {
    if (Sys.info()["sysname"] == "Windows")
    {
      tcltk::tkbind(tkwinviewplot, "<MouseWheel>", scroll_windows ) #This is the scroll event for Windows
    }
    else
    {
      tcltk::tkbind(tkwinviewplot, "<Button-4>", scrollup_callback ) #This is the scroll-up event for Unix-like systems
      tcltk::tkbind(tkwinviewplot, "<Button-5>", scrolldown_callback )  #This is the scroll-down event for Unix-like systems
    }
  }
  
  MouseMoved <- function(x,y)
  {
    current_pointer <- testEventinPlottingArea(x, y)
    
    if(!is.null(current_pointer))
    {
      MouseHover_callback(x = current_pointer[1], y = current_pointer[2])
      
      #Rectangle dragging move
      if(tkwinviewplot$plotting_data$leftMousePressed && !is.null(tkwinviewplot$plotting_data$last_press_coords))
      {
        #Delete previous rectangle
        if(!is.null(tkwinviewplot$plotting_data$selRect))
        {
          tcltk::tkdelete(tkwinviewplot, tkwinviewplot$plotting_data$selRect)
        }
        
        #Draw a new rectangle
        tkwinviewplot$plotting_data$selRect <- tcltk::tkcreate(tkwinviewplot, "rect", 
                                                               tkwinviewplot$plotting_data$last_press_coords[3], tkwinviewplot$plotting_data$last_press_coords[4], 
                                                               x, y,
                                                               outline = "#ff0000")
      }
    }
  }
  if(!is.null(MouseHover_callback))
  {
    tcltk::tkbind(tkwinviewplot, "<Motion>", MouseMoved ) 
  }
  
  #Handle resize events to redraw
  plotresize_callback <- function()
  {
    redrawPlotWidget(tkwinviewplot)
  }  
  tcltk::tkbind(tkwinviewplot, "<Configure>", plotresize_callback )
  
  return(tkwinviewplot)
}

#' redrawPlotWidget.
#' 
#' method to trigger the redraw of a tkrplot widget.
#'
#' @param tkplotwidget the widget to redraw.
#'
redrawPlotWidget <- function(tkplotwidget)
{
  win_width <- as.numeric(tcltk::tkwinfo("width", tkplotwidget)) 
  win_height <- as.numeric(tcltk::tkwinfo("height", tkplotwidget))
  
  png(filename =  tkplotwidget$file, width =  win_width, height = win_height, res = 100, pointsize = 10, type = "cairo" )
  par(xaxs = "i", yaxs = "i") #Force axes range to xlim and ylim provided, this should enable a better mouse clicks control
  try(tkplotwidget$fun())
  
  #Getting data space coords ranges
  data_range <- par("usr")
  
  #Getting ploting are percentages
  plot_range <- par("plt")
  
  dev.off()
  
  setPlotWidgetDataRanges(tkplotwidget, data_range, plot_range, win_width, win_height)
  
  tcltk::tkdelete(tkplotwidget, tkplotwidget$plotting_data$canvasImg)
  tkplotwidget$plotting_data$image <- tcltk::tkimage.create("photo", file =  tkplotwidget$file) 
  tkplotwidget$plotting_data$canvasImg <- tcltk::tkcreate(tkplotwidget, "image", 1, 1, image = tkplotwidget$plotting_data$image, anchor = "nw" )
  tcltk::tkitemlower(tkplotwidget, tkplotwidget$plotting_data$canvasImg )
}

#' setPlotWidgetDataRanges.
#' 
#' method to set tkrplot widget plotting data ranges.
#'
#' @param tkplotwidget a tkrplot widget.
#' @param data.range data range obtanied with par().
#' @param plot.range plot range from par().
#' @param imgWin_width widget horizionalt size.
#' @param imgWin_height widget vertical size.
#'
#' @return nothing, the environment is modified directely.
#'
setPlotWidgetDataRanges <- function(tkplotwidget, data.range, plot.range, imgWin_width, imgWin_height )
{
  #Getting data space coords ranges
  tkplotwidget$plotting_data$data.xmin <- data.range[1]
  tkplotwidget$plotting_data$data.xmax <- data.range[2]
  tkplotwidget$plotting_data$data.ymin <- data.range[3]
  tkplotwidget$plotting_data$data.ymax <- data.range[4]
  #cat(paste0("Data Xmin = ", data.xmin, "\n"))
  #cat(paste0("Data Xmax = ", data.xmax, "\n"))
  #cat(paste0("Data Ymin = ", data.ymin, "\n"))
  #cat(paste0("Data Ymax = ", data.ymax, "\n"))
  
  #Getting ploting are percentages
  tkplotwidget$plotting_data$plot.xmin.pixels <- plot.range[1] * imgWin_width
  tkplotwidget$plotting_data$plot.xmax.pixels <- plot.range[2] * imgWin_width
  tkplotwidget$plotting_data$plot.ymin.pixels <- (1-plot.range[3]) * imgWin_height #Includes mapping to tcltk win
  tkplotwidget$plotting_data$plot.ymax.pixels <- (1-plot.range[4]) * imgWin_height #Includes mapping to tcltk win
  
  #cat(paste0("Plot Xmin pixels = ", plot.xmin.pixels, "\n"))
  #cat(paste0("Plot Xmax pixels = ", plot.xmax.pixels, "\n"))
  #cat(paste0("Plot Ymin pixels = ", plot.ymin.pixels, "\n"))
  #cat(paste0("Plot Ymax pixels = ", plot.ymax.pixels, "\n"))
  
  #Calculate coords mapping slopes
  tkplotwidget$plotting_data$mx_device2plot <- (tkplotwidget$plotting_data$data.xmax - tkplotwidget$plotting_data$data.xmin)/(tkplotwidget$plotting_data$plot.xmax.pixels  - tkplotwidget$plotting_data$plot.xmin.pixels)
  tkplotwidget$plotting_data$my_device2plot <- (tkplotwidget$plotting_data$data.ymax - tkplotwidget$plotting_data$data.ymin)/(tkplotwidget$plotting_data$plot.ymax.pixels  - tkplotwidget$plotting_data$plot.ymin.pixels)
  
  #Calculate coords mapping offsets
  tkplotwidget$plotting_data$nx_device2plot <- tkplotwidget$plotting_data$data.xmin - tkplotwidget$plotting_data$mx_device2plot*tkplotwidget$plotting_data$plot.xmin.pixels
  tkplotwidget$plotting_data$ny_device2plot <- tkplotwidget$plotting_data$data.ymin - tkplotwidget$plotting_data$my_device2plot*tkplotwidget$plotting_data$plot.ymin.pixels
  
}

#' mapDeviceCoords2PlotCoords.
#' 
#' Obtain the plot coordinates from pixel coordinates in the plotting device.
#'
#' @param x pixel X coordinate in the plotting device.
#' @param y pixel Y coordinate in the plotting device.
#' @param plotting_data a list object with the plot coords information
#'
#' @return a vector with X,Y coordinates in the data plot space.
#' 
mapDeviceCoords2PlotCoords <- function( x, y, plotting_data)
{
  xplot <- plotting_data$mx_device2plot*x + plotting_data$nx_device2plot
  yplot <- plotting_data$my_device2plot*y + plotting_data$ny_device2plot
  return(c(xplot, yplot))
}
