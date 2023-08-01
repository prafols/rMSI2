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


#' LoadTwoMsImages
#'
#'A GUI to select up to two MS images to load in rMSI data format
#'
#' @return a list with the loaded rMSI objects
#'
LoadTwoMsImages <- function( )
{
  options(guiToolkit="tcltk") #force to use tcltk
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  img1_obj<-NULL
  img2_obj<-NULL
  load_completed<-F

  #Signal handlers
  fileChooseClicked <- function(h, ...)
  {
    fname<-gWidgets2::gfile("Select an MSI file to open", type="open", multi = F, filter =  list("XrMSI , imzML" = list(patterns = c("*.XrMSI", "*.imzML"))), initial.dir = path.expand("~/"))
    if(length(fname) == 0)
    {
      return ()
    }
    gWidgets2::svalue(h$action) <- fname
  }

  cancelClicked <- function(h, ...)
  {
    gWidgets2::dispose(h$obj)
  }

  okClicked <- function(h, ...)
  {

    #Check if images paths are correct and load data
    fpath1<- gWidgets2::svalue( this$txt_img1 )
    fpath2<- gWidgets2::svalue( this$txt_img2 )

    if( file.exists(fpath1) )
    {
      this$img1_obj<-LoadMsiData(data_file = fpath1, imzMLChecksum = F)
    }
    else
    {
      if(fpath1 != "")
      {
        gWidgets2::gmessage("Image 1 file is not valid or does not exists", "Error loading image 1", icon = "error")
      }
    }


    if( file.exists(fpath2) )
    {
      this$img2_obj<-LoadMsiData(data_file = fpath2, imzMLChecksum = F)
    }
    else
    {
      if(fpath2 != "")
    {
        gWidgets2::gmessage("Image 2 file is not valid or does not exists", "Error loading image 2", icon = "error")
      }
    }

    if( !file.exists(fpath1) && !file.exists(fpath2) )
    {
      gWidgets2::gmessage("At least one file must be selected.", "Error loading images", icon = "error")
    }
    else
    {
      this$load_completed <- T
      gWidgets2::dispose(h$obj)
    }
  }

  #Build GUI
  win_tp<-gWidgets2::gwindow("Load MSI data", visible = F, width = 600, height = 200, parent = c(0.5*as.integer(tcltk::tkwinfo("screenwidth", ".")) - 300, 0.5* as.integer(tcltk::tkwinfo("screenheight", ".")) - 100))
  box_tp<-gWidgets2::ggroup(horizontal = F,container = win_tp, spacing = 10)
  box_title <- gWidgets2::ggroup(horizontal = T,container = box_tp, spacing = 10, expand = T, fill = T)
  lbl_tl<-gWidgets2::glabel("Select up to two MS images to load in .imzML or .XrMSI format", container = box_title)
  gWidgets2::addSpring(box_title)

  frm_img1<-gWidgets2::gframe(spacing = 5, container = box_tp)
  box_img1<-gWidgets2::ggroup(horizontal = T, container = frm_img1, spacing = 10, expand = T, fill = T)
  lbl_img1<-gWidgets2::glabel("  Image 1:", container = box_img1)
  txt_img1<-gWidgets2::gedit(width = 60, container = box_img1, expand = T, fill = T)
  btn_img1<-gWidgets2::gbutton("Select file", container = box_img1, handler = fileChooseClicked, action = txt_img1)

  frm_img2<-gWidgets2::gframe(spacing = 5, container = box_tp)
  box_img2<-gWidgets2::ggroup(horizontal = T, container = frm_img2, spacing = 10, expand = T, fill = T)
  lbl_img2<-gWidgets2::glabel("  Image 2:", container = box_img2)
  txt_img2<-gWidgets2::gedit(width = 60, container = box_img2, expand = T, fill = T)
  btn_img2<-gWidgets2::gbutton("Select file", container = box_img2, handler = fileChooseClicked, action = txt_img2)

  frm_btns <- gWidgets2::gframe( spacing = 5, container = box_tp)
  box_btns <- gWidgets2::ggroup(horizontal = T, spacing = 10, container = frm_btns, expand = T, fill = T)
  gWidgets2::addSpring(box_btns)
  btn_cancel <- gWidgets2::gbutton("Cancel", container = box_btns, handler = cancelClicked)
  btn_proced <- gWidgets2::gbutton("OK", container = box_btns, handler = okClicked)
  
  gWidgets2::visible(win_tp)<-T
  ## Set the name for the class
  class(this) <- append(class(this),"LoadTwoMsImages")
  gc()
  
  #Do not return until this windos is disposed...
  while(gWidgets2::isExtant(this$win_tp ))
  {
    Sys.sleep(0.1)
  }

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}

