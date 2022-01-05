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

ImportWizardGui <- function()
{
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ##Class data members
  abort_GUI <- T
  ParamList <- list( procparams = ProcessingParameters(), datadesc = ImzMLDataDescription(), numberOfThread = max(parallel::detectCores() - 2, 2), memoryPerThreadMB = 100) #Start with defaults!
  workdir <- path.expand("~/")
  
    
  SetWorkingDir <- function( sDir )
  {
    this$workdir <- sDir
    
    MSIfilesdatapath <- this$browseMSIFile$GetPath()
    DataOutpath <- this$browseOut$GetPath()
    if( length(MSIfilesdatapath) > 0 && !is.null(DataOutpath))
    {
      gWidgets2::enabled(btn_run) <- T 
    }
  }
  
  GetWorkingDir <- function()
  {
    return(this$workdir)
  }
  
  CheckIfDataIsSet <- function( )
  {
    #Data source and output list
    MSIfilesdatapath <- this$browseMSIFile$GetPath()
    
    if( length(MSIfilesdatapath) == 0)
    {
      #Error no file selected
      gWidgets2::gmessage("Error: No file selected.", icon = "error") 
      return(FALSE)
    }
    
    for(i in 1:length(MSIfilesdatapath))
    {
      if( !file.exists(MSIfilesdatapath[i]) )
      {
        #Error file does not exists
        gWidgets2::gmessage(paste("Error: file", MSIfilesdatapath[i], "\nDoes not exists."), icon = "error")
        return(FALSE)
      }
    }
    
    DataOutpath <- this$browseOut$GetPath()
    if( is.null(DataOutpath) )
    {
      gWidgets2::gmessage("Error: No output directory selected.", icon = "error")
      return(FALSE)
    }
    if( !dir.exists(DataOutpath))
    {
      gWidgets2::gmessage("Error: output directory does not exists.", icon = "error")
      return(FALSE)
    }
    return(TRUE)
  }
  
  
  UpdateDataDescriptions <- function()
  {
    #Data source and output list
    MSIfilesdatapath <- this$browseMSIFile$GetPath()

    for(i in 1:length(MSIfilesdatapath))
    {
      if( this$xmlRoiFiles[i] == "")
      {
        this$ParamList$datadesc$appendImzMLDataPath(path_imzML = MSIfilesdatapath[i] )
      }
      else
      {
        this$ParamList$datadesc$appendImzMLDataPath(path_imzML = MSIfilesdatapath[i], subimage_roi_xml =this$xmlRoiFiles[i] )
      }
    }
    
    #Set the output data path
    this$ParamList$datadesc$setOutputPath(this$browseOut$GetPath())
  }
  
  UpdateParamListStruct <- function()
  {
    #Smoothing list
    this$ParamList$procparams$preprocessing$smoothing$enable <- gWidgets2::svalue(this$check_smoothing)
    this$ParamList$procparams$preprocessing$smoothing$kernelSize <- as.integer(gWidgets2::svalue(this$ratio_SGkernSize))
    
    #Alignment list
    this$ParamList$procparams$preprocessing$alignment$enable <- gWidgets2::svalue(this$check_alignment) 
    this$ParamList$procparams$preprocessing$alignment$iterations <- as.integer(gWidgets2::svalue(this$spin_AlignIterations))
    this$ParamList$procparams$preprocessing$alignment$maxShiftppm <- as.double(gWidgets2::svalue(this$spin_AlignMaxDev))
    this$ParamList$procparams$preprocessing$alignment$bilinear <- as.logical(gWidgets2::svalue(this$check_AlignBilinear))
    this$ParamList$procparams$preprocessing$alignment$refLow <- as.double(gWidgets2::svalue(this$spin_AlignRefLow))
    this$ParamList$procparams$preprocessing$alignment$refMid <- as.double(gWidgets2::svalue(this$spin_AlignRefMid))
    this$ParamList$procparams$preprocessing$alignment$refHigh <- as.double(gWidgets2::svalue(this$spin_AlignRefHigh))
    this$ParamList$procparams$preprocessing$alignment$overSampling <- as.integer(gWidgets2::svalue(this$spin_AlignOverSampling))
    
    #Calibration list
    this$ParamList$procparams$preprocessing$massCalibration <- gWidgets2::svalue(this$check_calibration)
    
    #Peak picking list
    this$ParamList$procparams$preprocessing$peakpicking$enable <- gWidgets2::svalue(this$check_peakpicking)
    this$ParamList$procparams$preprocessing$peakpicking$SNR <- as.double(gWidgets2::svalue(this$spin_SNR))
    this$ParamList$procparams$preprocessing$peakpicking$WinSize <- as.integer(gWidgets2::svalue(this$spin_peakWin))
    this$ParamList$procparams$preprocessing$peakpicking$overSampling <- as.integer(gWidgets2::svalue(this$spin_peakOversample))
    
    #Peak binning list
    this$ParamList$procparams$preprocessing$peakbinning$enable <- gWidgets2::svalue(this$check_peakbinnig)
    this$ParamList$procparams$preprocessing$peakbinning$tolerance <- as.double(gWidgets2::svalue(this$spin_binTolerance))
    
    if(gWidgets2::svalue(this$ratio_binningUnits) == "[ppm]")
    {
      this$ParamList$procparams$preprocessing$peakbinning$tolerance_in_ppm <- T
    }
    else
    {
      this$ParamList$procparams$preprocessing$peakbinning$tolerance_in_ppm <- F
    }
    this$ParamList$procparams$preprocessing$peakbinning$binFilter <- as.double(gWidgets2::svalue(this$spin_binFilter)/100)
    
    #Number of threads
    this$ParamList$numberOfThread <- as.integer( gWidgets2::svalue( this$spin_nThreads ))
    this$ParamList$memoryPerThreadMB <- as.integer( gWidgets2::svalue( this$spin_memXThread ) )
    
    #If datasets must be merged or not
    this$ParamList$procparams$preprocessing$merge <- gWidgets2::svalue(this$check_datamerge)
  }  
  
  RunClicked <- function(h, ...)
  {
    this$abort_GUI <- F
    UpdateParamListStruct()
    #Set the data_is_peaklist field
    this$ParamList$datadesc$setImzMLIsPeakList(gWidgets2::svalue(this$check_imzMLisPeakList))
    gWidgets2::dispose(h$obj)
  }
  
  ButtonLoadParamsClicked <- function(h, ...)
  {
    #Show a file selector widget to load params
    mPath <- gWidgets2::gfile( text = "Load rMSI2 processing parameters", type = "open", multi = F, filter = list("rMSI proc parameters" = list(patterns = c("*.params"))), initial.dir = this$workdir ) 
    if( length( mPath) > 0)
    {
      this$ParamList$procparams <- LoadProcParams( mPath )
      
      #Smoothing list
      gWidgets2::svalue(this$check_smoothing) <- this$ParamList$procparams$preprocessing$smoothing$enable
      gWidgets2::svalue(this$ratio_SGkernSize) <- this$ParamList$procparams$preprocessing$smoothing$kernelSize
      
      #Alignment list
      gWidgets2::svalue(this$check_alignment) <- this$ParamList$procparams$preprocessing$alignment$enable
      gWidgets2::svalue(this$spin_AlignIterations) <- this$ParamList$procparams$preprocessing$alignment$iterations 
      gWidgets2::svalue(this$spin_AlignMaxDev) <- this$ParamList$procparams$preprocessing$alignment$maxShiftppm 
      gWidgets2::svalue(this$check_AlignBilinear) <- this$ParamList$procparams$preprocessing$alignment$bilinear
      gWidgets2::svalue(this$spin_AlignRefLow) <- this$ParamList$procparams$preprocessing$alignment$refLow
      gWidgets2::svalue(this$spin_AlignRefMid) <- this$ParamList$procparams$preprocessing$alignment$refMid
      gWidgets2::svalue(this$spin_AlignRefHigh) <- this$ParamList$procparams$preprocessing$alignment$refHigh
      gWidgets2::svalue(this$spin_AlignOverSampling) <- this$ParamList$procparams$preprocessing$alignment$overSampling
      
      #Calibration list
      gWidgets2::svalue(this$check_calibration) <- this$ParamList$procparams$preprocessing$massCalibration
      
      #Peak picking list
      gWidgets2::svalue(this$check_peakpicking) <- this$ParamList$procparams$preprocessing$peakpicking$enable
      gWidgets2::svalue(this$spin_SNR) <- this$ParamList$procparams$preprocessing$peakpicking$SNR 
      gWidgets2::svalue(this$spin_peakWin) <- this$ParamList$procparams$preprocessing$peakpicking$WinSize 
      gWidgets2::svalue(this$spin_peakOversample) <- this$ParamList$procparams$preprocessing$peakpicking$overSampling
      
      #Peak binning list
      gWidgets2::svalue(this$check_peakbinnig) <- this$ParamList$procparams$preprocessing$peakbinning$enable
      gWidgets2::svalue(this$spin_binTolerance) <- this$ParamList$procparams$preprocessing$peakbinning$tolerance 
      if(this$ParamList$procparams$preprocessing$peakbinning$tolerance_in_ppm )
      {
        gWidgets2::svalue(this$ratio_binningUnits) <- "[ppm]"
      }
      else
      {
        gWidgets2::svalue(this$ratio_binningUnits) <- "[scans]"
      }
      
      gWidgets2::svalue(this$spin_binFilter) <- this$ParamList$procparams$preprocessing$peakbinning$binFilter * 100
    }
  }
  
  ButtonSaveParamsClicked <- function(h, ...)
  {
    UpdateParamListStruct()
    mPath <- gWidgets2::gfile( text = "Save rMSI2 processing parameters", type = "save", multi = F, filter = list("rMSI proc parameters" = list(patterns = c("*.params"))), initial.dir = this$workdir ) 
    if( length( mPath) > 0)
    {
      StoreProcParams(mPath, this$ParamList$procparams)
    }
  }
  
  #Checkboxes enabled status
  ChkBoxImzMLisPeakList <- function(...)
  {
    gWidgets2::enabled(this$frm_SGKernSize) <- !gWidgets2::svalue(this$check_imzMLisPeakList)
    gWidgets2::enabled(this$frm_alignment) <- !gWidgets2::svalue(this$check_imzMLisPeakList)
    gWidgets2::enabled(this$frm_calibration) <- !gWidgets2::svalue(this$check_imzMLisPeakList)
    gWidgets2::enabled(this$frm_peakpick) <- !gWidgets2::svalue(this$check_imzMLisPeakList)
    
    if( gWidgets2::svalue(this$check_imzMLisPeakList) )
    {
      #Force to enable peak-binning when data is imzML peaklist
      gWidgets2::svalue(this$check_peakbinnig) <- T
      gWidgets2::enabled(this$check_peakbinnig) <- F
    }
    else
    {
      gWidgets2::enabled(this$check_peakbinnig) <- T 
    }
  }
  
  #Checkboxes enabled status
  ChkBoxSmoothingChanged <- function(...)
  {
    gWidgets2::enabled( this$lblSg) <- gWidgets2::svalue(this$check_smoothing)
    gWidgets2::enabled( this$ratio_SGkernSize) <- gWidgets2::svalue(this$check_smoothing)
  }
  
  ChkBoxAlignmentChanged <- function(...)
  {
    gWidgets2::enabled( this$spin_AlignIterations) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$spin_AlignMaxDev) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$check_AlignBilinear) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$spin_AlignRefLow) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$spin_AlignRefMid) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$spin_AlignRefHigh) <- gWidgets2::svalue(this$check_alignment)
    gWidgets2::enabled( this$spin_AlignOverSampling) <- gWidgets2::svalue(this$check_alignment)
  }
  
  ChkBoxCalibrationChanged <- function(...)
  {
    #Nothing here currently, keeping it for future features...
  }
  
  ChkBoxPeakPickingChanged <- function(...)
  {
    gWidgets2::enabled( this$spin_SNR) <- gWidgets2::svalue(this$check_peakpicking)
    gWidgets2::enabled( this$spin_peakWin) <- gWidgets2::svalue(this$check_peakpicking)
    gWidgets2::enabled( this$spin_peakOversample) <- gWidgets2::svalue(this$check_peakpicking)
    
    if( !gWidgets2::svalue(this$check_peakpicking) && !gWidgets2::svalue(this$check_imzMLisPeakList) )
    {
      #Force to disable peak-binning when not peak-picking or data is not imzML peaklist
      gWidgets2::svalue(this$check_peakbinnig) <- F
    }
  }
  
  ChkBoxPeakBinningChanged <- function(...)
  {
    gWidgets2::enabled( this$spin_binTolerance) <- gWidgets2::svalue(this$check_peakbinnig)
    gWidgets2::enabled( this$spin_binFilter) <- gWidgets2::svalue(this$check_peakbinnig)
    gWidgets2::enabled( this$frm_binUnits) <- gWidgets2::svalue(this$check_peakbinnig)
  }
  
  mainWin<-gWidgets2::gwindow(title = "MSI data import and process wizard", visible = F)
  box_mainH <- gWidgets2::ggroup(horizontal = T, container = mainWin, expand = T, fill = T)
  box_mainV <- gWidgets2::ggroup(horizontal = F, container = box_mainH, expand = , fill = T)
  
  #Param management
  frm_paramManager <-  gWidgets2::gframe( "Processing Parameters Management", container =  box_mainV, expand = F, fill = F)
  box_paramManager <- gWidgets2::ggroup( horizontal = T, container = frm_paramManager, expand = F, fill = F)
  btn_paramLoad <- gWidgets2::gbutton("Load parameters", handler = this$ButtonLoadParamsClicked, container = box_paramManager)
  btn_paramSave <- gWidgets2::gbutton("Save parameters", handler = this$ButtonSaveParamsClicked, container = box_paramManager)
  
  #Data Input box
  frm_dataInput <-  gWidgets2::gframe( "Data Source", container =  box_mainV, expand = T, fill = T)
  box_dataInput <- gWidgets2::ggroup( horizontal = F, container = frm_dataInput, expand = T, fill = T)
  check_imzMLisPeakList <- gWidgets2::gcheckbox("ImzML files contain peaks-lists", checked = F, container = box_dataInput, handler = this$ChkBoxImzMLisPeakList)
  browseMSIFile <- FileBrowseWidget( box_dataInput, setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  #Set initial state properly
  this$browseMSIFile$ClearPath()
  this$browseMSIFile$SetFileFilter("imzML")
  this$browseMSIFile$SetLabel("imzML files:")
  
  #Pre-processing box
  frm_preProcessing <- gWidgets2::gframe( "Pre-Processing parameters", container = box_mainH, expand = T, fill = T )
  box_procH <- gWidgets2::ggroup(horizontal = T, container = frm_preProcessing, expand = T, fill = T, spacing = 20)
  box_proc1 <- gWidgets2::ggroup(horizontal = F, container = box_procH, expand = T, fill = T, spacing = 20)
  box_proc2 <- gWidgets2::ggroup(horizontal = F, container = box_procH, expand = T, fill = T, spacing = 20)
  drawLabelSpin <- function( parent_widget, sText, minVal, maxVal, defaultVal, decPlaces = 0, increments = 1 )
  {
    box_spin <- gWidgets2::ggroup(horizontal = T, container = parent_widget)
    gWidgets2::glabel(sText, container = box_spin )
    gWidgets2::addSpring(box_spin)
    return (gWidgets2::gspinbutton(from = minVal, to = maxVal, by = increments, value = defaultVal, container = box_spin, digits = decPlaces)  )
  }
  
  #Smoothing params
  frm_SGKernSize <- gWidgets2::gframe("Savitzky-Golay Smoothing:", container = box_proc1)
  box_smoothing <- gWidgets2::ggroup(horizontal = F, container = frm_SGKernSize)
  check_smoothing <- gWidgets2::gcheckbox("Enable smoothing", checked = this$ParamList$procparams$preprocessing$smoothing$enable, container = box_smoothing, handler = this$ChkBoxSmoothingChanged )
  lblSg <- gWidgets2::glabel("Savitzky-Golay kernel size:", container = box_smoothing)
  smoothingPosibleValues <- c(5,7,9,11,13,15)
  smoothingSelValue <- which.min(abs(smoothingPosibleValues - this$ParamList$procparams$preprocessing$smoothing$kernelSize))
  ratio_SGkernSize<- gWidgets2::gradio(smoothingPosibleValues, selected = smoothingSelValue, horizontal = T, container = box_smoothing)
    
  #Alignment params
  frm_alignment <- gWidgets2::gframe("Alignment", container = box_proc1, spacing = 10)
  box_alignment <- gWidgets2::ggroup(horizontal = F, container = frm_alignment)
  check_alignment <- gWidgets2::gcheckbox("Enable alignment", checked = this$ParamList$procparams$preprocessing$alignment$enable, container = box_alignment, handler = this$ChkBoxAlignmentChanged)
  check_AlignBilinear <- gWidgets2::gcheckbox("Bilinear mode", checked = this$ParamList$procparams$preprocessing$alignment$bilinear, container = box_alignment )
  spin_AlignIterations <- drawLabelSpin(box_alignment, "Iterations:", 1, 5, this$ParamList$procparams$preprocessing$alignment$iterations) 
  spin_AlignMaxDev <- drawLabelSpin(box_alignment, "Max Shift [ppm]:", 10, 1000, this$ParamList$procparams$preprocessing$alignment$maxShiftppm)
  spin_AlignRefLow <- drawLabelSpin(box_alignment, "Ref. Bottom:", 0, 1, this$ParamList$procparams$preprocessing$alignment$refLow, decPlaces = 2, increments = 0.05)
  spin_AlignRefMid <- drawLabelSpin(box_alignment, "Ref. Center:", 0, 1, this$ParamList$procparams$preprocessing$alignment$refMid, decPlaces = 2, increments = 0.05)
  spin_AlignRefHigh <- drawLabelSpin(box_alignment, "Ref. Top:", 0, 1, this$ParamList$procparams$preprocessing$alignment$refHigh, decPlaces = 2, increments = 0.05)
  spin_AlignOverSampling <- drawLabelSpin(box_alignment, "Over-Sampling:", 1, 10, this$ParamList$procparams$preprocessing$alignment$overSampling, decPlaces = 0)
  
  #Mass Calibration params
  frm_calibration <- gWidgets2::gframe("Mass Calibration", container = box_proc1, spacing = 10)
  box_calibration <- gWidgets2::ggroup(horizontal = F, container = frm_calibration)
  check_calibration <- gWidgets2::gcheckbox("Enable calibration", checked = this$ParamList$procparams$preprocessing$massCalibration, container = box_calibration, handler = this$ChkBoxCalibrationChanged)
  
  #Peak picking params
  frm_peakpick <- gWidgets2::gframe("Peak-Picking", container = box_proc2, spacing = 10)
  box_peakpick <- gWidgets2::ggroup(horizontal = F, container = frm_peakpick)
  check_peakpicking <- gWidgets2::gcheckbox("Enable peak-picking", checked = this$ParamList$procparams$preprocessing$peakpicking$enable, container = box_peakpick, handler = this$ChkBoxPeakPickingChanged)
  spin_SNR <- drawLabelSpin(box_peakpick, "SNR Threshold:", 1, 100, this$ParamList$procparams$preprocessing$peakpicking$SNR)
  spin_peakWin <- drawLabelSpin(box_peakpick, "Detector window size:", 5, 200, this$ParamList$procparams$preprocessing$peakpicking$WinSize)
  spin_peakOversample <- drawLabelSpin(box_peakpick, "Peak shape over-sampling:", 1, 50, this$ParamList$procparams$preprocessing$peakpicking$overSampling)
  
  #Peak binning params
  frm_peakbinning <- gWidgets2::gframe("Peak-Binning", container = box_proc2, spacing = 10)
  box_peakbinning <- gWidgets2::ggroup(horizontal = F, container = frm_peakbinning)
  check_peakbinnig <- gWidgets2::gcheckbox("Enable peak-binning", checked = this$ParamList$procparams$preprocessing$peakbinning$enable, container = box_peakbinning, handler = this$ChkBoxPeakBinningChanged)
  spin_binTolerance <- drawLabelSpin(box_peakbinning, "Peak-Bin Tolerance:", 1, 1000, this$ParamList$procparams$preprocessing$peakbinning$tolerance, decPlaces = 0)
  frm_binUnits <- gWidgets2::gframe("Binning Tolerance Units:", container = box_peakbinning)
  if(this$ParamList$procparams$preprocessing$peakbinning$tolerance_in_ppm)
  {
    peakBinningToleranceModeSelection <- 1
  }
  else
  {
    peakBinningToleranceModeSelection <- 2
  }
  ratio_binningUnits <- gWidgets2::gradio(c("[ppm]", "[scans]"), container = frm_binUnits, selected = peakBinningToleranceModeSelection, horizontal = T)
  spin_binFilter <- drawLabelSpin(box_peakbinning, "Peak Filter [%]:", 1, 100, this$ParamList$procparams$preprocessing$peakbinning$binFilter * 100)

  #Number of processing threads
  frm_procThreads <- gWidgets2::gframe("Processing Threads", container = box_proc2, spacing = 10)
  threads_Vbox <- gWidgets2::ggroup(horizontal = F, container = frm_procThreads)
  spin_nThreads <- drawLabelSpin(threads_Vbox, "Max Threads:", 1, parallel::detectCores(), ParamList$numberOfThread, decPlaces = 0)
  spin_memXThread <- drawLabelSpin(threads_Vbox, "Memory for each thread [MB]:", 10, 1024, ParamList$memoryPerThreadMB, decPlaces = 0)
  
  #Data output box
  gWidgets2::addSpring(box_mainV)
  frm_dataOutput <- gWidgets2::gframe( "Data Output", container = box_mainV, expand = T, fill = T)
  browseOut <- FileBrowseWidget( frm_dataOutput, sLabel = "Output directory:", dirSel = T, fFilter = "txt", setdir_fun = SetWorkingDir, getdir_fun = GetWorkingDir )
  
  #Data merge
  check_datamerge <- gWidgets2::gcheckbox("Merged processing", checked =  this$ParamList$procparams$getMergedProcessing(), container = box_proc2)
  
  #Run button
  btn_run <- gWidgets2::gbutton("Process Data", handler = this$RunClicked, container = box_proc2)
  gWidgets2::enabled(btn_run) <- F
  
  gWidgets2::visible(mainWin) <- T
  
  ## Set the name for the class
  class(this) <- append(class(this),"DataProcessWindow")
  gc()
  
  #Do not return until this window is disposed...
  while(gWidgets2::isExtant(this$mainWin ))
  {
    Sys.sleep(0.1)
  }
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  if(this$abort_GUI)
  {
    #User aborted the GUI
    return(NULL)
  }  
  
  #Ask the user for the XML's containing the ROI files
  if(!this$CheckIfDataIsSet())
  {
    #Data description contains errors
    return(NULL)
  }
  
  data_imzML_files <- browseMSIFile$GetPath()
  xmlRoiFiles <- XmlRoiSelectionDialog(basename(data_imzML_files), init_dir = dirname(data_imzML_files[1] ) )
  if(is.null(xmlRoiFiles))
  {
    #Process aborted by user
    cat("Processing aborted\n")
    return(NULL)
  }
  else
  {
    #Append data and ROI's
    UpdateDataDescriptions()
  }
  
  #Return a structured list with all input data
  return( ParamList )
}

XmlRoiSelectionDialog <- function( img_names, init_dir = getwd())
{
  oldWarning<-options()$warn
  options(warn = -1)
  
  this <- environment()
  initial_dir <- init_dir
  abort_process <- T
  xmlList_subimg <- rep("", length(img_names))
  
  browseButtonClicked <- function (evt, ...)
  {
    mPath <- gWidgets2::gfile( text = "Select a ROI XML file",
                               type = "open", 
                               multi = F,
                               filter = list("XML file"  = list(patterns = c("*xml", "*XML"))), 
                               initial.dir = this$initial_dir)
    if( length( mPath) > 0)
    {
      if(evt$action$source == "subimg")
      {
        #Setting the subimg widget
        RGtk2::gtkEntrySetText(gWidgets2::getToolkitWidget(this$selFilesEntry_list[[evt$action$img]]$subImg), basename(mPath)) 
        this$xmlList_subimg[evt$action$img] <- mPath
      }
      
      else
      {
        stop("Error: browseButtonClicked() in XmlRoiSelectionDialog(). Invaldi source field\n")
      }
      this$initial_dir <- dirname(mPath)
    }
  }
  
  imgRoiBrowseWidget <- function( text, parent_widget, index )
  {
    frm <- gWidgets2::gframe( container = parent_widget)
    hboxSubImg <- gWidgets2::ggroup(horizontal = T, container = frm, expand = T)
    lblSubImg <- gWidgets2::glabel(paste0(text, "   ROI's sub-images:"), container = hboxSubImg)
    gWidgets2::addSpring(hboxSubImg)
    entrySubImg <- gWidgets2::gedit( container = hboxSubImg, width = 30 )
    btnSubImg <- gWidgets2::gbutton("Select", handler = this$browseButtonClicked, action = list(img = index, source = "subimg"), container = hboxSubImg)
    
    return(list(subImg = entrySubImg))
  }
  
  OkButtonClicked <- function (h, ...)
  {
    this$abort_process <- F
    gWidgets2::dispose(h$obj)
  }
  
  AbortButtonClicked <- function (h, ...)
  {
    gWidgets2::dispose(h$obj)
  }
  
  dlgMain <- gWidgets2::gwindow("Select ROI's for each image (optional)")
  gWidgets2::size(dlgMain) <- c(800, 480)
  main_vBox <- gWidgets2::ggroup(horizontal = F, container = dlgMain)
  
  #An informative label
  lbl_info <- gWidgets2::glabel( "Add a ROI XML file for each image. Or just click Ok to proceed without ROI information", container = main_vBox)
  
  #Fill the list of image
  xml_vBox <- gWidgets2::ggroup(horizontal = F, use.scrollwindow = T, container = main_vBox, expand=T)
  selFilesEntry_list <- list()
  for( i in 1:length(img_names))
  {
    selFilesEntry_list[[i]]<-imgRoiBrowseWidget(img_names[i], xml_vBox, i)
  }
  
  #The buttons
  btn_hBox <- gWidgets2::ggroup(horizontal = T, container = main_vBox)
  gWidgets2::addSpring(btn_hBox)
  btn_Ok <- gWidgets2::gbutton("Ok", handler = this$OkButtonClicked, container = btn_hBox)
  btn_Abort <- gWidgets2::gbutton("Abort", handler = this$AbortButtonClicked, container = btn_hBox)
  
  gWidgets2::visible(dlgMain) <- T
  
  ## Set the name for the class
  class(this) <- append(class(this),"RoiSelWindow")
  gc()
  
  #Do not return until this window is disposed...
  while(gWidgets2::isExtant(this$dlgMain ))
  {
    Sys.sleep(0.1)
  }

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
  
  if(this$abort_process)
  {
    #Process aborted by user
    return(NULL)
  }
  
  return( this$xmlList_subimg )
}

