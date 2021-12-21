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
#include "mtaverage.h"
using namespace Rcpp;

MTAverage::MTAverage(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis, double minTIC, double maxTIC) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis),
  TICmin(minTIC), 
  TICmax(maxTIC)
{
  AverageSpectrum = NumericVector(ioObj->getMassAxisLength());
  for(int i = 0; i < AverageSpectrum.length(); i++)
  {
    AverageSpectrum[i] = 0.0;
  }

  validPixelCount = new unsigned int[ioObj->getNumberOfCubes()];  
  for(int i = 0; i < ioObj->getNumberOfCubes() ; i++)
  {
    validPixelCount[i] = 0;
  }
  
  //Get normalizations
  TicNormalizations = new double*[rMSIObj_list.size()];
  for(int i = 0; i < rMSIObj_list.size(); i++)
  {
    NumericVector tics = Rcpp::as<NumericVector>(Rcpp::as<Rcpp::DataFrame>((Rcpp::as<Rcpp::List>(rMSIObj_list[i]))["normalizations"])["TIC"]);
    TicNormalizations[i] = new double[tics.length()];
    memcpy(TicNormalizations[i], tics.begin(), tics.length() * sizeof(double));
  }
}

MTAverage::~MTAverage()
{
  for (int i = 0; i < ioObj->get_images_count(); i++)
  {
    delete[] TicNormalizations[i];
  }
  delete[] TicNormalizations;
  
  delete[] validPixelCount;
}

NumericVector MTAverage::Run()
{
  Rcpp::Rcout<<"Calculating overall average spectrum...\n";
  
  //Run in multi-threading
  runMSIProcessingCpp();
  
  //Sum the total number of pixels used to calculate the average
  unsigned int pixelCount = 0;
  for(int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
    pixelCount += validPixelCount[i];
  }
  
  if(pixelCount > 0)
  {
    for(int i = 0; i < ioObj->getMassAxisLength(); i++)
    {
      AverageSpectrum[i] /= (double)pixelCount;
    }
  }
  
  return AverageSpectrum;
}


void MTAverage::ProcessingFunction(int threadSlot)
{
  
  double *partialAverage = new double[cubes[threadSlot]->ncols];
  for(int i = 0; i < cubes[threadSlot]->ncols; i++)
  {
    partialAverage[i] = 0.0;
  }
  
  //Perform the average value of each mass channel in the current loaded cube
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    int imageIndex = ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    int pixelIndex = ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
    
    double TICval = TicNormalizations[imageIndex][pixelIndex];
    
    if(TICval >= TICmin && TICval <= TICmax)
    {
      for (int k= 0; k < cubes[threadSlot]->ncols; k++)
      {
        partialAverage[k] += (cubes[threadSlot]->dataInterpolated[j][k])/TICval; //Average with TIC Normalization
        validPixelCount[cubes[threadSlot]->cubeID]++;
      }
    }
  }
  
  //Partial average
  if(validPixelCount[cubes[threadSlot]->cubeID] > 0)
  {
    averageMutex.lock();
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      AverageSpectrum[k] += partialAverage[k];
    }
    averageMutex.unlock();
  }
  
  delete[] partialAverage;
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
NumericVector COverallAverageSpectrum(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB,
                               Rcpp::NumericVector commonMassAxis,
                               double minTIC, double maxTic)
{
  NumericVector out;
  try
  {
    MTAverage myAverage(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB, 
                      commonMassAxis,
                      minTIC,
                      maxTic);
 
    out = myAverage.Run();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}
