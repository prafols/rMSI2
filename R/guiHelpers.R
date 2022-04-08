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

#' createPlotWidget.
#' wrapper function to create tkrplots widgets integrated with gwidgets2
#'
#' @param parent gWidgets2 parent widget.
#' @param redraw_function function to draw the plot.
#'
#' @return a plot widget.
#'
createPlotWidget <- function(parent, redraw_function )
{
  tkwinviewplot <- tkrplot::tkrplot(gWidgets2::getToolkitWidget(parent), redraw_function,  hscale = 1, vscale = 1 )
  tcltk::tkpack(tkwinviewplot)
  return(tkwinviewplot)
}

#' redrawPlotWidget.
#' redraws the contents of a plot widget
#'
#' @param tkplotwidget the plot widget to redraw.
#'
redrawPlotWidget <- function(tkplotwidget)
{
  invisible(tkrplot::tkrreplot(tkplotwidget))
}
