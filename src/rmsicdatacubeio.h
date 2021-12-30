/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2014 Pere Rafols Soler
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/
#ifndef RMSI_DATA_CUBE_IO_H
  #define RMSI_DATA_CUBE_IO_H

#include <string>
#include <Rcpp.h>
#include "imzMLBin.h"
#include "peakpicking.h" //needed for peak list definition

/********************************************************************************
 *  CrMSIDataCubeIO: A C++ class to extract datacubes from multiple imzML files
 *  
 *  This class provides a simple interface to load and save data from imzML files
 *  using the rMSIproc multithread processing approach.
 ********************************************************************************/

#define MAX_DATACUBES_ROWS_WITH_PEAKLISTS_ONLY 1000
//When reading peak lists without the spectral interpolation neigther a common mass axis the max number of rows used in each data cube is set by this define.
//This is the case for peak binning first stage where only peak lists are accessed. 
//The default value of a thousand is a good trade of to maximize parallelization.

typedef enum DataCubeIOMode
{
  DATA_READ, //Read spectral data with interpolation to the common mass axis
  DATA_STORE, //Store spectral
  PEAKLIST_READ, //Read a peaklist (interpolated data will be empty) 
  PEAKLIST_STORE, //Store a peaklist
  DATA_AND_PEAKLIST_READ //Read spectral data with interpolation to the common mass axis along with its corresponding peak list... 
}DataCubeIOMode;

class CrMSIDataCubeIO
{
  public:
    
    // Constructor Arguments:
    // - massAxis: The common mass axis for all the data.
    // - cubeMemoryLimitMB: Memory limit for the interpolated spectra in a cube, thus the acutal used memory can be higher due to stored data as-is in the imzML.
    // - dataModeEnum: set the data store mode.
    CrMSIDataCubeIO(Rcpp::NumericVector massAxis, double cubeMemoryLimitMB, DataCubeIOMode dataModeEnum, Rcpp::String imzMLOutputPath = "");
    ~CrMSIDataCubeIO();
    
    //Struct to define a whole data cube in memory
    typedef struct
    {
      int cubeID;
      int ncols;
      int nrows;
      imzMLSpectrum *dataOriginal; //Pointer to multiple imzMLSpectrum structs 
      PeakPicking::Peaks **peakLists; //Pointer to the peaklists assosiated with a datacube
      double **dataInterpolated;   
    } DataCube;
    
    //Appends an image to be processed.
    // - rMSIObj: The rMSI object describing a complete MSI image.
    // - outputImzMLuuid: UUID for the output imzML file. Only used if storeData is set in the constructor.
    // - outputImzMLsuffix: suffix for the output imzML filenmaes. Only used if storeData is set in the constructor.
    void appedImageData(Rcpp::List rMSIOoj,  std::string outputImzMLuuid = "", std::string outputImzMLfname = "");
    
    //Loads a data cube specified by iCube into data_ptr
    //WARNING: This is not a thread-safe method, it must be only used on the main thread who'll take care of loading data from HDD and copy it to other thread mem space.
    //It returns a pointer to an allocated structure containing the datacube.
    //It is responisability of user to free memmory of the loaded datacube using freeDataCube() function.
    DataCube *loadDataCube( int iCube);
    void freeDataCube(DataCube *data_ptr);
    
    //Stores a datacube to the path assosiated with its ID
    void storeDataCube(int iCube, DataCube *data_ptr);
    
    //Return the total number of cubes in the ramdisk
    int getNumberOfCubes();
    
    //Return the length of the common mass axis
    int getMassAxisLength();
    
    //Return the total number of pixels in each cube
    int getNumberOfPixelsInCube(int iCube);
    
    //Return the image index (imzMLreader/writer for a given pixels specified as cube id and pixel row)
    int getImageIndex(int iCube, int cubeRow);
    
    //Return the pixel id in an image (imzMLreader/writer for a given pixels specified as cube id and pixel row)
    int getPixelId(int iCube, int cubeRow);
    
    //Return the corresponding peak Matrix row for a specific row in an cube
    unsigned int getPeakMatrixRow(int iCube, int cubeRow);
    
    //Return a Data Frame with the imzML offsets of a specified imzMLWriter
    Rcpp::DataFrame get_OffsetsLengths(unsigned int index);
    
    //Return the number of images attached to the data cube object
    unsigned int get_images_count();

    //Return the average spectrum for a specified imzMLWriter
    Rcpp::NumericVector get_AverageSpectrum(unsigned int index);
    
    //Return the base spectrum for a specified imzMLWriter
    Rcpp::NumericVector get_BaseSpectrum(unsigned int index);
    
    //Return true if a peaklist is the imzMLPeaksReaders array is in rMSI format (include snr and area)
    bool get_all_peakLists_are_rMSIformated();

  private:
    DataCubeIOMode dataMode; //An enum to set the data mode.
    std::string dataOutputPath; //A path to save output imzML data
    unsigned int cubeMaxNumRows; //The maximum rows in a datacube calculated from the maximum memory allowed by each cube and the mass axis length.
    Rcpp::NumericVector mass; //A common mass axis for all images to process
    std::vector<ImzMLBinRead*> imzMLReaders;  //Pointers to multiple imzMLReadrs initialized with openIbd = false to avoid exiding the maximum open files.
    std::vector<ImzMLBinWrite*> imzMLWriters; //Pointers to multiple imzMLWriters initialized with openIbd = false to avoid exiding the maximum open files.
    std::vector<ImzMLBinRead*> imzMLPeaksReaders;  //Pointers to multiple imzMLReadrs for peaklist reading initialized with openIbd = false to avoid exiding the maximum open files.
    std::vector<Rcpp::NumericVector> acumulatedSpectrum; //A vector to contain all the average spectra (there is one for each imzMLWriter)
    std::vector<Rcpp::NumericVector> baseSpectrum; //A vector to contain all the base spectra (there is one for each imzMLWriter)  
    
    unsigned int next_peakMatrix_row; //A counter to follow added peak matrix rows
    
    //Struct to internally handle data cube accessors
    typedef struct
    {
      int pixel_ID; //ID of a pixel in the imzML file
      int imzML_ID; //integer pointing to a imzML Handle (reader and writer) for a pixel ID, -1 value is allowed to indicate to an unallocated imzML
      unsigned int peakMatrix_row; //The row index (C style begining with zero) the pixel will ocuppy in the peak matrix
    } PixelDescription;
    
    typedef std::vector<PixelDescription> DataCubeDescription; //Each data cube is defined as standard vector of pixel descriptions
    
    std::vector<DataCubeDescription> dataCubesDesc; //A vector to describe all data cubes in the data set
};

#endif
