#########################################################################
#     rMSIproc - R package for MSI data processing
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

FileBrowseWidget <- function( parent_widget, sLabel = "File:", dirSel = F, multiSel = T, fFilter = "txt", setdir_fun = NULL, getdir_fun = NULL )
{
  options(guiToolkit="tcltk") #force to use tcltk
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ## Class data members 
  directorySelection <- dirSel
  fileFilter <- fFilter
  sPath <- c()
  MultiFileSelect <- multiSel
  GetWorkDir <- getdir_fun
  SetWorkDir <- setdir_fun
  
  ## Class methods
  BrowseDialog <- function(evt, ...)
  {
    if(this$directorySelection)
    {
      sTitle <- "Select a directory"
      sType <- "selectdir"
      lFilter <- list("All dirs" = list(patterns = c("*")))
    }
    else
    {
      sTitle <- "Select a file"
      sType <- "open"
      lFilter <- list(list(patterns = c(paste("*.", this$fileFilter, sep =""), paste("*.", toupper(this$fileFilter), sep =""))))
      names(lFilter) <- paste(as.character(gWidgets2::svalue(this$lbl)), "*.", this$fileFilter)
    }
    
    init_dir <- path.expand("~/")
    if(!is.null(this$GetWorkDir))
    {
      init_dir <- this$GetWorkDir()
    }
    mPath <- gWidgets2::gfile( text = sTitle, type = sType, multi = (!this$directorySelection) & this$MultiFileSelect, filter = lFilter, initial.dir = init_dir  ) 
    if( length( mPath) > 0)
    {
      this$sPath <- mPath
      if(this$MultiFileSelect && !this$directorySelection)
      {
        mPathTxt <- ""
        for( i in 1:length(this$sPath))
        {
          mPathTxt <- paste(mPathTxt, basename(this$sPath[i]),sep ="")
          if( i < length(this$sPath))
          {
            mPathTxt <- paste(mPathTxt, "\n", sep ="")
          }
        }
        gWidgets2::svalue(this$entry_files) <- mPathTxt
        gWidgets2::svalue(this$entry_dir) <-  dirname(this$sPath[1])
      }
      else
      {
        gWidgets2::svalue(this$entry_dir) <-  this$sPath
      }
    }
    if(!is.null(this$SetWorkDir))
    {
      if(dir.exists(gWidgets2::svalue(this$entry_dir)))
      {
        this$SetWorkDir( gWidgets2::svalue(this$entry_dir)) #It is already a dir so not using the dirname funciton
      }
      else
      {
        this$SetWorkDir( dirname(gWidgets2::svalue(this$entry_dir)) )
      }
    }
  }
  
  ##Public Methods
  SetDirSelection <- function( isDirSelection )
  {
    this$directorySelection <- isDirSelection
    gWidgets2::visible(this$entry_files) <- !isDirSelection
  }
  
  SetFileFilter <- function ( fFilter )
  {
    this$fileFilter <- fFilter
  }
  
  SetLabel <- function (sLabel)
  {
    gWidgets2::svalue(this$lbl) <- sLabel
  }
  
  GetPath <- function()
  {
    return(this$sPath)
  }
  
  SetEnabled <- function( enabled )
  {
    gWidgets2::enabled( this$frm ) <- enabled  
  }
  
  ClearPath <- function()
  {
    this$sPath <- c()
    gWidgets2::svalue(this$entry_dir) <-  ""
    gWidgets2::svalue(this$entry_files) <- ""
  }
    
  ## Class constructor
  frm <- gWidgets2::gframe( container = parent_widget, fill = T, expand = T)
  box <- gWidgets2::ggroup(horizontal = T, container = frm, expand = T, fill = T)
  box_lbl <- gWidgets2::ggroup(horizontal = F, container = box)
  lbl <- gWidgets2::glabel(sLabel, container = box_lbl)

  box_entry <- gWidgets2::ggroup(horizontal = F, container = box, expand = T, fill = T)
  entry_dir<-gWidgets2::gedit(container = box_entry) #The width argument has no effect on tcltk 
  
  #set it not editable by disconnecting the key event
  tcltk::tkbind(gWidgets2::getToolkitWidget(entry_dir), "<KeyPress>", "break")
  
  
  if( multiSel && !dirSel )
  {
    entry_files<-gWidgets2::gtext(container = box_entry) #The width argument has no effect on tcltk 
    #the width of the text entry defaults to 80 even if setting it to 50 with the following tcl command
    #tcltk::tcl(gWidgets2::getToolkitWidget(entry_files), "configure", "-width", 50)
    tcltk::tkbind(gWidgets2::getToolkitWidget(entry_files), "<KeyPress>", "break")
  }
  box_btn <- gWidgets2::ggroup(horizontal = F, container = box)
  btn <- gWidgets2::gbutton("Browse", handler = this$BrowseDialog, container = box_btn )
  
  ## Set the name for the class
  class(this) <- append(class(this),"FileBrowseWidget")
  gc()
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  return(this)
}
