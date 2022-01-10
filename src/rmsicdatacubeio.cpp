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

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 
#include <stdexcept>
#include "rmsicdatacubeio.h"
using namespace Rcpp;

CrMSIDataCubeIO::CrMSIDataCubeIO(Rcpp::NumericVector massAxis, double cubeMemoryLimitMB, DataCubeIOMode dataModeEnum, Rcpp::String imzMLOutputPath)
  :mass(massAxis), dataMode(dataModeEnum), dataOutputPath(imzMLOutputPath.get_cstring()), next_peakMatrix_row(0)
{
  if(mass.length() > 0)
  {
    cubeMaxNumRows = std::ceil((1024*1024*cubeMemoryLimitMB)/(double)(8*mass.length()));
  }
  else
  {
    cubeMaxNumRows = MAX_DATACUBES_ROWS_WITH_PEAKLISTS_ONLY;
  }
}

CrMSIDataCubeIO::~CrMSIDataCubeIO()
{
  for( int i = 0; i < imzMLReaders.size(); i++)
  {
    delete imzMLReaders[i];
  }
  
  for( int i = 0; i < imzMLWriters.size(); i++)
  {
    delete imzMLWriters[i];
  }
  
  for( int i = 0; i < imzMLPeaksReaders.size(); i++)
  {
    delete imzMLPeaksReaders[i];
  }
}

void CrMSIDataCubeIO::appedImageData(Rcpp::List rMSIobj, std::string outputImzMLuuid, std::string outputImzMLfname)
{
  if(outputImzMLuuid.empty() && (dataMode == DataCubeIOMode::DATA_STORE || dataMode == DataCubeIOMode::PEAKLIST_STORE ))
  {
    throw std::runtime_error("Error: A valid UUID must be provided for the stored imzML file.\n");
  }
  
  if(outputImzMLfname.empty() && (dataMode == DataCubeIOMode::DATA_STORE || dataMode == DataCubeIOMode::PEAKLIST_STORE ))
  {
    throw std::runtime_error("Error: A filename suffix must be provided for the stored imzML file.\n");
  }
  
  if(outputImzMLuuid.length() != 32 && (dataMode == DataCubeIOMode::DATA_STORE || dataMode == DataCubeIOMode::PEAKLIST_STORE ))
  {
    throw std::runtime_error("Error: The UUID must contain 16 bytes.\n");
  }
  
  //get the rMSI object data field
  List data = rMSIobj["data"];
  List imzML;
  DataFrame imzMLrun;
  
  //Set the imzML reader
  if(dataMode != DataCubeIOMode::PEAKLIST_READ)
  {
    imzML = data["imzML"];
    imzMLrun = as<DataFrame>(imzML["run"]);
    std::string sFilePath = as<std::string>(data["path"]);
    std::string sFnameImzML = as<std::string>(imzML["file"]);
    std::string sFnameImzMLInput = sFilePath + "/" + sFnameImzML + ".ibd";
    
    //Rcpp::Rcout << "CrMSIDataCubeIO::appedImageData()--> sFnameImzML = "<< sFnameImzML << std::endl; //DEBUG line
  
    //Append the imzML Reader
    imzMLReaders.push_back(new ImzMLBinRead(sFnameImzMLInput.c_str(), 
                                       imzMLrun.nrows(), 
                                       as<String>(imzML["mz_dataType"]),
                                       as<String>(imzML["int_dataType"]) ,
                                       (as<bool>(imzML["continuous_mode"])),
                                       false, //Do not call the file open() on constructor
                                       false //Used to read spectral data so just keep the default
                                       )); 
    
    //If data is in continuous mode but resampling is needed, then read the imzML in processed mode to enable interpolation.
    imzMLReaders.back()->setCommonMassAxis(mass.length(), mass.begin());
  
    NumericVector imzML_mzLength = imzMLrun["mzLength"];
    NumericVector imzML_mzOffsets = imzMLrun["mzOffset"];
    NumericVector imzML_intLength = imzMLrun["intLength"];
    NumericVector imzML_intOffsets = imzMLrun["intOffset"];
    imzMLReaders.back()->set_mzLength(&imzML_mzLength);  
    imzMLReaders.back()->set_mzOffset(&imzML_mzOffsets);
    imzMLReaders.back()->set_intLength(&imzML_intLength);
    imzMLReaders.back()->set_intOffset(&imzML_intOffsets);
  }
  
  //Create the imzMLPeaksReaders corresponding to the current imzMLread
  DataFrame peaksimzMLrun;
  if(dataMode == DataCubeIOMode::PEAKLIST_READ || dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ)
  {
    if( data["peaklist"] == R_NilValue)
    {
      throw std::runtime_error("Error: peak list data not available\n");
    }
    
    List peakListimzML = data["peaklist"];
    peaksimzMLrun = as<DataFrame>(peakListimzML["run_data"]);
    if(dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ)
    {
      if(peaksimzMLrun.nrows() != imzMLrun.nrows())
      {
        throw std::runtime_error("Error: peak list number of pixels do not match with spectral data number of pixels\n");
      }
    }
    
    std::string speaksFilePath = as<std::string>(peakListimzML["path"]);
    std::string speaksFnameImzML = as<std::string>(peakListimzML["file"]);
    std::string speaksFnameImzMLInput = speaksFilePath + "/" + speaksFnameImzML + ".ibd";
    
    //Append the imzML peaks Reader
    imzMLPeaksReaders.push_back(new ImzMLBinRead(speaksFnameImzMLInput.c_str(), 
                                                 peaksimzMLrun.nrows(), 
                                                 as<String>(peakListimzML["mz_dataType"]),
                                                 as<String>(peakListimzML["int_dataType"]),
                                                 false,
                                                 false, //Do not call the file open() on constructor
                                                 as<bool>(peakListimzML["rMSIpeakList"])//true if the peak list is in rMSI format
                                                 )); 
    
    NumericVector imzML_mzLength = peaksimzMLrun["mzLength"];
    NumericVector imzML_mzOffsets = peaksimzMLrun["mzOffset"];
    NumericVector imzML_intLength = peaksimzMLrun["intLength"];
    NumericVector imzML_intOffsets = peaksimzMLrun["intOffset"];
    imzMLPeaksReaders.back()->set_mzLength(&imzML_mzLength);  
    imzMLPeaksReaders.back()->set_mzOffset(&imzML_mzOffsets);
    imzMLPeaksReaders.back()->set_intLength(&imzML_intLength);
    imzMLPeaksReaders.back()->set_intOffset(&imzML_intOffsets);
  }
  
  //Create the imzMLwriter corresponding to the current imzMLreader
  if(dataMode == DataCubeIOMode::DATA_STORE || dataMode == DataCubeIOMode::PEAKLIST_STORE )
  {
    std::string sFnameImzMLOutput= dataOutputPath + "/" + outputImzMLfname + ".ibd"; 
    
    if(dataMode == DataCubeIOMode::DATA_STORE)
    {
      //Outputing spectral data
      imzMLWriters.push_back(new ImzMLBinWrite(sFnameImzMLOutput.c_str(),
                                            imzMLrun.nrows(), 
                                            as<String>(imzML["mz_dataType"]),
                                            as<String>(imzML["int_dataType"]) ,
                                            as<bool>(imzML["continuous_mode"]),
                                            true, //Data cubes will be sotred using the sequential mode
                                            false)); //Dont call the file open() on constructor to allow directely writing the uuid, then I'll close it!
      
      //Only perform average accumulations when processing spectral data
      acumulatedSpectrum.push_back(Rcpp::NumericVector(mass.length()));
      baseSpectrum.push_back(Rcpp::NumericVector(mass.length()));
    }
    
    if(dataMode == DataCubeIOMode::PEAKLIST_STORE)
    {
      //Outputing peaklist data
      imzMLWriters.push_back(new ImzMLBinWrite(sFnameImzMLOutput.c_str(),
                                               imzMLrun.nrows(), 
                                               "double",
                                               "double",
                                               false,
                                               true, //Data cubes will be sotred using the sequential mode
                                               false)); //Dont call the file open() on constructor to allow directely writing the uuid, then I'll close it!
    }
    
    imzMLWriters.back()->open(true); //Open and truncate the file
    imzMLWriters.back()->writeUUID(outputImzMLuuid);
    imzMLWriters.back()->close(); 
  }
  
  //Initialize the cube description if this is the first call to appedImageData()
  if(dataCubesDesc.size() == 0)
  {
    dataCubesDesc.push_back(DataCubeDescription());
  }
  
  unsigned int iters;
  if( dataMode == DataCubeIOMode::PEAKLIST_READ)
  {
    iters = peaksimzMLrun.nrows(); 
  }
  else
  {
    iters = imzMLrun.nrows();
  }
  for( unsigned int i = 0; i < iters; i++) //For each pixel in imzML file
  {
    if(dataCubesDesc.back().size() == cubeMaxNumRows)
    {
      //No more space left in the cube, so create a new cube
      dataCubesDesc.push_back(DataCubeDescription());
    }
    
    dataCubesDesc.back().push_back(PixelDescription());
    dataCubesDesc.back().back().pixel_ID = i;
    if( dataMode == DataCubeIOMode::PEAKLIST_READ)
    {
      dataCubesDesc.back().back().imzML_ID = imzMLPeaksReaders.size()-1;
    }
    else
    {
      dataCubesDesc.back().back().imzML_ID = imzMLReaders.size()-1;
    }
    dataCubesDesc.back().back().peakMatrix_row = next_peakMatrix_row;
    next_peakMatrix_row++;
  }
}

CrMSIDataCubeIO::DataCube *CrMSIDataCubeIO::loadDataCube(int iCube)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> iCube=" << iCube << std::endl; //DEBUG line
  
  DataCube *data_ptr = new DataCube();
  data_ptr->cubeID = iCube;
  data_ptr->nrows = dataCubesDesc[iCube].size();
  data_ptr->ncols = mass.length();
  
  if(dataMode != DataCubeIOMode::PEAKLIST_READ)
  {
    data_ptr->dataOriginal = new imzMLSpectrum[data_ptr->nrows];
    data_ptr->dataInterpolated = new double*[data_ptr->nrows];
  }
  
  if(dataMode == DataCubeIOMode::PEAKLIST_STORE || dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ || dataMode == DataCubeIOMode::PEAKLIST_READ)
  {
    data_ptr->peakLists = new PeakPicking::Peaks*[data_ptr->nrows];
  }
  
  //Data reading
  int current_imzML_id;
  int previous_imzML_id = -1; //Start previous as -1 to indicate an unallocated imzML
  for(unsigned int i = 0; i < data_ptr->nrows; i++) //For each spectrum belonging to the selected datacube
  {
    current_imzML_id = dataCubesDesc[iCube][i].imzML_ID;
    
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> current_imzML_id=" << current_imzML_id << std::endl; //DEBUG line!
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> mzLength(0)=" << imzMLReaders[current_imzML_id]->get_mzLength(0)  << std::endl; //DEBUG line
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> ibd file=" << imzMLReaders[current_imzML_id]->getIbdFilePath()  << std::endl; //DEBUG line
    
    if(current_imzML_id != previous_imzML_id)
    {
      if(previous_imzML_id != -1)
      {
        if( (dataMode == DataCubeIOMode::PEAKLIST_READ) || (dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ))
        {
          imzMLPeaksReaders[previous_imzML_id]->close();
        }
        if(dataMode != DataCubeIOMode::PEAKLIST_READ)
        {
          imzMLReaders[previous_imzML_id]->close();
        }
      }
      
      if( (dataMode == DataCubeIOMode::PEAKLIST_READ) || (dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ))
      {
        imzMLPeaksReaders[current_imzML_id]->open();
      }
      if(dataMode != DataCubeIOMode::PEAKLIST_READ)
      {
        imzMLReaders[current_imzML_id]->open();
      }
    }
    
    if(dataMode != DataCubeIOMode::PEAKLIST_READ)
    {
      data_ptr->dataInterpolated[i] = new double[data_ptr->ncols];
      data_ptr->dataOriginal[i] = imzMLReaders[current_imzML_id]->ReadSpectrum(dataCubesDesc[iCube][i].pixel_ID, //pixel id to read
                                                  0, //unsigned int ionIndex
                                                  mass.length(),//unsigned int ionCount
                                                  data_ptr->dataInterpolated[i], //Store data directely at the datacube mem
                                                  false //Disable auto-interpolation
                                                  );
    }
    
    if(dataMode == DataCubeIOMode::PEAKLIST_READ || dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ)
    {
      //Load the peak list in reading mode
      data_ptr->peakLists[i] = imzMLPeaksReaders[current_imzML_id]->ReadPeakList(dataCubesDesc[iCube][i].pixel_ID); 
    }
                                  
    previous_imzML_id = current_imzML_id;
  }
  
  //Force to close the last opened imzML
  if( (dataMode == DataCubeIOMode::PEAKLIST_READ) || (dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ))
  {
    imzMLPeaksReaders[current_imzML_id]->close();
  }
  if(dataMode != DataCubeIOMode::PEAKLIST_READ)
  {
    imzMLReaders[current_imzML_id]->close();
  }
  
  return data_ptr;
}

void CrMSIDataCubeIO::freeDataCube(DataCube *data_ptr)
{
  for( int i = 0; i < data_ptr->nrows; i++ )
  {
    if(dataMode != DataCubeIOMode::PEAKLIST_READ)
    {
      delete[] data_ptr->dataInterpolated[i];
    }
    if(dataMode == DataCubeIOMode::PEAKLIST_STORE || dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ || dataMode == DataCubeIOMode::PEAKLIST_READ)
    {
      delete data_ptr->peakLists[i];
    }
  }
  
  if(dataMode != DataCubeIOMode::PEAKLIST_READ)
  {
    delete[] data_ptr->dataOriginal;
    delete[] data_ptr->dataInterpolated;
  }
  
  if(dataMode == DataCubeIOMode::PEAKLIST_STORE || dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ || dataMode == DataCubeIOMode::PEAKLIST_READ)
  {
    delete[] data_ptr->peakLists;
  }
  delete data_ptr;
}

void CrMSIDataCubeIO::interpolateDataCube(DataCube *data_ptr)
{
  if(dataMode != DataCubeIOMode::PEAKLIST_READ)
  {
    if(data_ptr->cubeID >= dataCubesDesc.size())
    {
      throw std::runtime_error("Error: DataCube index out of range\n");
    }
    
    int current_imzML_id;
    for(unsigned int i = 0; i < data_ptr->nrows; i++) //For each spectrum belonging to the selected datacube
    {
      current_imzML_id = dataCubesDesc[data_ptr->cubeID][i].imzML_ID;
      imzMLReaders[current_imzML_id]->InterpolateSpectrum( &(data_ptr->dataOriginal[i]), 0, mass.length(), data_ptr->dataInterpolated[i]);
    }
  }
}

void CrMSIDataCubeIO::storeDataCube(DataCube *data_ptr) 
{
  if(data_ptr->cubeID >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  int current_imzML_id;
  int previous_imzML_id = -1; //Start previous as -1 to indicate an unallocated imzML
  
  for(unsigned int i = 0; i < data_ptr->nrows; i++) //For each spectrum belonging to the selected datacube
  {
    current_imzML_id = dataCubesDesc[data_ptr->cubeID][i].imzML_ID;
    
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> current_imzML_id=" << current_imzML_id << std::endl; //DEBUG line!
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> mzLength(0)=" << imzMLReaders[current_imzML_id]->get_mzLength(0)  << std::endl; //DEBUG line
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> ibd file=" << imzMLReaders[current_imzML_id]->getIbdFilePath()  << std::endl; //DEBUG line
    
    if(current_imzML_id != previous_imzML_id)
    {
      if(previous_imzML_id != -1)
      {
        imzMLWriters[previous_imzML_id]->close();
      }
      imzMLWriters[current_imzML_id]->open();
    }
    
    //Store Spectral data
    if( dataMode == DataCubeIOMode::DATA_STORE)
    {
      //Calc average and base spectrum
      for(unsigned int j = 0; j < mass.length(); j++)
      {
        acumulatedSpectrum[current_imzML_id][j] += data_ptr->dataInterpolated[i][j];  
        baseSpectrum[current_imzML_id][j] = data_ptr->dataInterpolated[i][j] > baseSpectrum[current_imzML_id][j] ? data_ptr->dataInterpolated[i][j] : baseSpectrum[current_imzML_id][j];
      }
      
      if(imzMLWriters[current_imzML_id]->get_continuous())
      {
        //Continuous mode write
        imzMLWriters[current_imzML_id]->writeMzData(mass.length(), mass.begin()); //The mass axis will be writtem only once
        imzMLWriters[current_imzML_id]->writeIntData(mass.length(), data_ptr->dataInterpolated[i]);
      }
      else
      {
        //Processed mode write
        imzMLWriters[current_imzML_id]->writeMzData(data_ptr->dataOriginal[i].imzMLmass.size(), data_ptr->dataOriginal[i].imzMLmass.data());
        imzMLWriters[current_imzML_id]->writeIntData(data_ptr->dataOriginal[i].imzMLintensity.size(), data_ptr->dataOriginal[i].imzMLintensity.data());
      }
    }
    
    //Store Peak lists
    if( dataMode == DataCubeIOMode::PEAKLIST_STORE)
    {
      //Processed mode write with phantom data for snr, area and bin size
      //Rcpp::Rcout<<"iCube = " <<iCube<< "  i = " <<i<<"  mass.size() = "<<data_ptr->peakLists[i]->mass.size()<< "  intensity.size() = " << data_ptr->peakLists[i]->intensity.size() <<"\n";
      imzMLWriters[current_imzML_id]->writePeakList(data_ptr->peakLists[i]->mass.size(), 
                                                    data_ptr->peakLists[i]->mass.data(),
                                                    data_ptr->peakLists[i]->intensity.data(),
                                                    data_ptr->peakLists[i]->area.data(),
                                                    data_ptr->peakLists[i]->SNR.data(),
                                                    data_ptr->peakLists[i]->binSize.data() );
    }
    
    previous_imzML_id = current_imzML_id;
  }

  imzMLWriters[current_imzML_id]->close(); //Force to close the last opened imzML 
}

int CrMSIDataCubeIO::getNumberOfCubes()
{
  return dataCubesDesc.size();
}

int CrMSIDataCubeIO::getMassAxisLength()
{
  return mass.length();
}

int CrMSIDataCubeIO::getNumberOfPixelsInCube(int iCube)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  return dataCubesDesc[iCube].size();
}

Rcpp::DataFrame CrMSIDataCubeIO::get_OffsetsLengths(unsigned int index)
{
  if(index >= imzMLWriters.size())
  {
    throw std::runtime_error("Error: imzMLWriter index out of range\n");
  }
  return imzMLWriters[index]->get_OffsetsLengths();
}

unsigned int CrMSIDataCubeIO::get_images_count()
{
 return imzMLReaders.size();
}

int CrMSIDataCubeIO::getImageIndex(int iCube, int cubeRow)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  if( cubeRow >= dataCubesDesc[iCube].size())
  {
    throw std::runtime_error("Error: DataCube row out of range\n");
  }
  return dataCubesDesc[iCube][cubeRow].imzML_ID;
}

int CrMSIDataCubeIO::getPixelId(int iCube, int cubeRow)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  if( cubeRow >= dataCubesDesc[iCube].size())
  {
    throw std::runtime_error("Error: DataCube row out of range\n");
  }
  return dataCubesDesc[iCube][cubeRow].pixel_ID;
}

unsigned int CrMSIDataCubeIO::getPeakMatrixRow(int iCube, int cubeRow)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  if( cubeRow >= dataCubesDesc[iCube].size())
  {
    throw std::runtime_error("Error: DataCube row out of range\n");
  }
  return dataCubesDesc[iCube][cubeRow].peakMatrix_row;
}

Rcpp::NumericVector CrMSIDataCubeIO::get_AverageSpectrum(unsigned int index)
{
  if(index >= acumulatedSpectrum.size())
  {
    throw std::runtime_error("Error: acumulatedSpectrum index out of range\n");
  }
  
  Rcpp::NumericVector AverageSpectrum(mass.length());
  for( unsigned int i = 0; i < mass.length(); i++)
  {
    AverageSpectrum[i] = acumulatedSpectrum[index][i] / ((double) (imzMLWriters[index]->get_number_of_pixels()));
  }
  
  return AverageSpectrum;
}

Rcpp::NumericVector CrMSIDataCubeIO::get_BaseSpectrum(unsigned int index)
{
  if(index >= baseSpectrum.size())
  {
    throw std::runtime_error("Error: baseSpectrum index out of range\n");
  }
  
  return baseSpectrum[index];
}

bool CrMSIDataCubeIO::get_all_peakLists_are_rMSIformated()
{
  bool rMSIpeakListFormated = true;
  if( (dataMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ) || (dataMode == DataCubeIOMode::PEAKLIST_READ) )
  {
    for(int i = 0; i < imzMLPeaksReaders.size(); i++)
    {
      rMSIpeakListFormated &= imzMLPeaksReaders[i]->get_rMSIPeakListFormat(); 
    }
    
  }
  else
  {
    throw std::runtime_error("Error: get_all_peakLists_are_rMSIformated() can only be called when data mode is set to DATA_AND_PEAKLIST_READ or PEAKLIST_READ\n");
  }
    
  return rMSIpeakListFormated;
}

