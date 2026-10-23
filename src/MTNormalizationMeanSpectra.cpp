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
#include <algorithm>    // std::fill
#include <future>
#include "MTNormalizationMeanSpectra.h"
using namespace Rcpp;

MTNormalizationMeanSpectra::MTNormalizationMeanSpectra(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis), 
  rMSIObj_lst(rMSIObj_list)
{
  num_of_pixels.resize(rMSIObj_lst.length());
  threadBuffers.resize(numOfThreadsDouble);
  
  for( int threadSlot = 0; threadSlot < numOfThreadsDouble; threadSlot++)
  {
    threadBuffers[threadSlot].thread_average.resize(ioObj->get_images_count());
    threadBuffers[threadSlot].thread_base.resize(ioObj->get_images_count());
    
    for(unsigned int i = 0; i < ioObj->get_images_count(); i++)
    {
      threadBuffers[threadSlot].thread_average[i].resize(massAxis.length());
      threadBuffers[threadSlot].thread_base[i].resize(massAxis.length());
      
      std::fill(threadBuffers[threadSlot].thread_average[i].begin(), threadBuffers[threadSlot].thread_average[i].end(), 0.0);
      std::fill(threadBuffers[threadSlot].thread_base[i].begin(), threadBuffers[threadSlot].thread_base[i].end(), 0.0);
    }
  }
  
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
   
  }
  
  averageSpectrum.assign(rMSIObj_lst.length(), std::vector<double>(massAxis.length(), 0.0));
  baseSpectrum.assign(rMSIObj_lst.length(), std::vector<double>(massAxis.length(), 0.0));
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
  
  //Final acumulator of mean and base spectra
  unsigned int numImages = ioObj->get_images_count();
  int numMasses = massAxis.length();
  size_t batchSize = static_cast<size_t>(numOfThreadsDouble);
  
  for (size_t iStart = 0; iStart < numImages; iStart += batchSize)
  {
    size_t iEnd = std::min(iStart + batchSize, static_cast<size_t>(numImages));
    std::vector<std::future<void>> futures;
    
    // Launch a batch
    for (size_t i = iStart; i < iEnd; i++)
    {
      futures.push_back(std::async(std::launch::async, [this, i, numMasses]() {
        
        for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
        {
          for(int k = 0; k < numMasses; k++)
          {
            // Mean
            averageSpectrum[i][k] += threadBuffers[iThread].thread_average[i][k];
            
            // Base
            if(threadBuffers[iThread].thread_base[i][k] > baseSpectrum[i][k])
            {
              baseSpectrum[i][k] = threadBuffers[iThread].thread_base[i][k];
            }
          }
        }
        
      }));
    }
    
    // Esperar a que termine este lote antes de pasar al siguiente
    for(auto& f : futures)
    {
      f.get();
    }
  }
  
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
    (as<List>(rMSIObj_lst[i])["mean"]) = Rcpp::wrap(averageSpectrum[i]);
    (as<List>(rMSIObj_lst[i])["base"]) = Rcpp::wrap(baseSpectrum[i]); 
    (as<List>(rMSIObj_lst[i])["normalizations"]) = NormsDF;
  }
  
  return rMSIObj_lst; 
}


void MTNormalizationMeanSpectra::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  double TIC, RMS, MAX;
  
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    TIC = 0.0;
    RMS = 0.0;
    MAX = 0.0;
    int imgID = ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    int pixelID = ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
    

    auto& avg_buf = threadBuffers[threadSlot].thread_average[imgID];
    auto& base_buf = threadBuffers[threadSlot].thread_base[imgID];
    auto& data_row = cubes[threadSlot]->dataInterpolated[j];
    double inv_num_pixels = 1.0 / ((double)(num_of_pixels[imgID]));
    
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      double val = data_row[k];
      TIC += val;
      RMS += (val * val);
      if (val > MAX) MAX = val;
      
      avg_buf[k] += val * inv_num_pixels;
      if (val > base_buf[k])
      {
        base_buf[k] = val;
      }
    }
    RMS = sqrt(RMS);
    
    Normalizations[imgID][pixelID].TIC = TIC;
    Normalizations[imgID][pixelID].RMS = RMS;
    Normalizations[imgID][pixelID].MAX = MAX;
  }
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
