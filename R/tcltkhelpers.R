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

#' gcombobox_rMSI.
#' 
#' Build on top of gwidgets2tcltk implementation to properly set the widget width.
#'
#' @param items 
#' @param selected 
#' @param container 
#' @param handler 
#'
#' @return a gcombobox widget
#'
gcombobox_rMSI <- function(items, selected = 1, container = NULL, handler = NULL)
{
  char_length <- max(unlist(lapply(items, nchar))) + 2
  combo_widget <- gWidgets2::gcombobox( items = items, selected = selected, container = container, handler = handler)  
  tcltk::tkconfigure(gWidgets2::getToolkitWidget(combo_widget), width = char_length)
  return(combo_widget)
}

#' gspinbutton_rMSI.
#' 
#' Build on top of gwidgets2tcltk implementation to properly set the widget width.
#'
#' @param from 
#' @param to 
#' @param by
#' @param value
#' @param digits
#' @param container 
#' @param handler 
#'
#' @return a gcombobox widget
#'
gspinbutton_rMSI <- function(from = 0, to = 10, by = 1, value = from, digits = 0, container = NULL, handler = NULL)
{
  char_length <- digits + 4
  spin_widget <- gWidgets2::gspinbutton( from = from, to = to, by = by, value = value, digits = digits, container = container, handler = handler)  
  tcltk::tkconfigure(gWidgets2::getToolkitWidget(spin_widget), width = char_length)
  return(spin_widget)
}

#' gbutton_icon.
#' 
#' button with an image as icon.
#'
#' @param image_file 
#' @param container 
#' @param handler 
#'
#' @return a gbutton widget
#' 
gbutton_icon <- function(image_file, container, handler)
{
  icon_img <- tcltk::tkimage.create("photo", file = image_file)
  btn <- gWidgets2::gbutton("", container = container, handler = handler)
  tcltk::tkconfigure(gWidgets2::getToolkitWidget(btn), image = icon_img, width = 3 )
  return(btn)
}

#' coloredCheckBox.
#' 
#' A speciall checkbox using only tcltk to set colors and font types.
#'
#' @param text 
#' @param checked 
#' @param handler 
#' @param container 
#' @param foreground 
#' @param background 
#' @param bold 
#'
#' @return
#'
coloredCheckBox <- function(text, checked = F, handler = NULL, container, foreground = NULL, background = NULL, bold = F)
{
  frm <- gWidgets2::ggroup(container = container, horizontal = T, spacing = 0, expand = F, fill = F)
  
  chb <- tcltk::tkcheckbutton(gWidgets2::getToolkitWidget(frm), text = text)
  chb$value <- tcltk::tclVar("0")
  chb$frame <- frm #Store the frame along with the widget to allow gWidgets2 commands on the frame
  tcltk::tkconfigure(chb, variable = chb$value)
  
  if(!is.null(foreground))
  {
    tcltk::tkconfigure(chb, foreground = foreground, activeforeground = foreground)
  }
  
  if(!is.null(background))
  {
    tcltk::tkconfigure(chb, background = background)
  }
  
  if(bold)
  {
    tcltk::tkconfigure(chb, font = tcltk::tkfont.create(weight="bold"))
  }
  
  if(checked)
  {
    tcltk::tkselect(chb)  
  }
  else
  {
    tcltk::tkdeselect(chb)
  }
  
  tcltk::tkpack(chb)
  
  if(!is.null(handler))
  {
    tcltk::tkbind(chb, "<ButtonRelease-1>", handler )
  }
  
  return(chb)
}

#' getValue_coloredCheckBox.
#' 
#' get the value of a tcltk checkbox as boolean.
#'
#' @param widget 
#'
#' @return checkboz current value
#'
getValue_coloredCheckBox <- function(widget)
{
  return( as.logical(as.integer(tcltk::tclvalue(widget$value))) )
}
