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
#include <cmath>
#include "MTNormalizationMeanSpectra.h"
using namespace Rcpp;

MTNormalizationMeanSpectra::MTNormalizationMeanSpectra(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis), 
  rMSIObj_lst(rMSIObj_list)
{
  averageSpectrum.resize(rMSIObj_lst.length());
  baseSpectrum.resize(rMSIObj_lst.length());
  num_of_pixels.resize(rMSIObj_lst.length());
  
  for( int i = 0; i < rMSIObj_lst.length(); i++)
  {
    num_of_pixels[i] = (as<NumericMatrix>(as<List>(rMSIObj_lst[i])["pos"])).nrow();
   
    Normalizations.push_back(std::vector<PixelNorms>(num_of_pixels[i]));
    for( int j = 0; j < num_of_pixels[i]; j++)
    {
      Normalizations[i][j].MAX = 0.0;
      Normalizations[i][j].TIC = 0.0;
      Normalizations[i][j].RMS = 0.0;
    }
   
    averageSpectrum[i] = NumericVector(massAxis.length());
    baseSpectrum[i] = NumericVector(massAxis.length());
  }
}

MTNormalizationMeanSpectra::~MTNormalizationMeanSpectra()
{
  //Empty destructor
}

List MTNormalizationMeanSpectra::Run()
{
  Rcpp::Rcout<<"Calculating normalizations...\n";
  
  //Run in multi-threading
  runMSIProcessingCpp();
  
  //Add averages and norms to each image and return the rMSIobj
  for( int i = 0; i < rMSIObj_lst.length(); i++)
  {
    NumericVector vTic(num_of_pixels[i]);
    NumericVector vRms(num_of_pixels[i]);
    NumericVector vMax(num_of_pixels[i]);
    for( int j = 0; j < num_of_pixels[i]; j++)
    {
      vTic[j] = Normalizations[i][j].TIC;
      vRms[j] = Normalizations[i][j].RMS;
      vMax[j] = Normalizations[i][j].MAX;
    }
    
    DataFrame NormsDF = DataFrame::create( Named("TIC") = vTic, Named("RMS") = vRms, Named("MAX") = vMax);
    (as<List>(rMSIObj_lst[i])["mean"]) = averageSpectrum[i]; 
    (as<List>(rMSIObj_lst[i])["base"]) = baseSpectrum[i]; 
    (as<List>(rMSIObj_lst[i])["normalizations"]) = NormsDF;
  }
  
  return rMSIObj_lst; 
}


void MTNormalizationMeanSpectra::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  double TIC, RMS, MAX;
  std::vector<std::vector<double>> thread_average(ioObj->get_images_count());
  std::vector<std::vector<double>> thread_base(ioObj->get_images_count());
  for(unsigned int i = 0; i < ioObj->get_images_count(); i++)
  {
    thread_average[i].resize(cubes[threadSlot]->ncols);
    thread_base[i].resize(cubes[threadSlot]->ncols);
  }
  
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    TIC = 0.0;
    RMS = 0.0;
    MAX = 0.0;
    int imgID = ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    int pixelID = ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
    
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      TIC += cubes[threadSlot]->dataInterpolated[j][k];
      RMS += (cubes[threadSlot]->dataInterpolated[j][k] * cubes[threadSlot]->dataInterpolated[j][k]);
      MAX = cubes[threadSlot]->dataInterpolated[j][k] > MAX ? cubes[threadSlot]->dataInterpolated[j][k] : MAX;
      
      thread_average[imgID][k] += cubes[threadSlot]->dataInterpolated[j][k] / ((double)(num_of_pixels[imgID]));
      thread_base[imgID][k] = cubes[threadSlot]->dataInterpolated[j][k] > thread_base[imgID][k] ? cubes[threadSlot]->dataInterpolated[j][k] : thread_base[imgID][k];
    }
    RMS = sqrt(RMS);
    
    Normalizations[imgID][pixelID].TIC = TIC;
    Normalizations[imgID][pixelID].RMS = RMS;
    Normalizations[imgID][pixelID].MAX = MAX;
  }
  
  mutex_copyData.lock();
  for(unsigned int i = 0; i < ioObj->get_images_count(); i++)
  {
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      averageSpectrum[i][k] += thread_average[i][k];
      baseSpectrum[i][k] = thread_base[i][k] > baseSpectrum[i][k] ? thread_base[i][k] : baseSpectrum[i][k];
    }
  }
  mutex_copyData.unlock();
  
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
List CNormalizationsAndMeans(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB,
                               Rcpp::NumericVector commonMassAxis)
{
  List out;
  
  try
    {
    MTNormalizationMeanSpectra myNorms(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB,
                      commonMassAxis);
    out = myNorms.Run();
    }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}
