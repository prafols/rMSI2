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
#include <stdexcept>
#include "mtinternalref.h"
using namespace Rcpp;

MTInternalRef::MTInternalRef(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis, Rcpp::NumericVector reference) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis), ref(reference)
{
  if(commonMassAxis.length() != reference.length())
  {
    throw std::runtime_error("ERROR: mass axis and the reference spectrum have a different length.\n");
  }
  
  //Compute the mean and the sdev of the reference
  ref_mean = computeMean(ref.begin());
  
  ref_sdev = 0.0;
  for(unsigned int i = 0; i < ref.length(); i++)
  {
    ref_sdev += (ref[i] - ref_mean) * (ref[i] - ref_mean);
  }
  ref_sdev = sqrt(ref_sdev);
  
  maxCorGlobal.cor = 0.0;
  maxCorGlobal.pixelID = -1;
  maxCorGlobal.imageID = -1;
}

MTInternalRef::~MTInternalRef()
{
  
}

List MTInternalRef::Run()
{
  Rcpp::Rcout<<"Calculating internal reference spectrum...\n";
  
  //Run in multi-threading
  runMSIProcessingCpp();
  
  return List::create(Named("score") = maxCorGlobal.cor, Named("imgIndex") = maxCorGlobal.imageID + 1, Named("ID") = maxCorGlobal.pixelID + 1 ); //Transform to R indexing here
}


void MTInternalRef::ProcessingFunction(int threadSlot)
{
  double x;
  double current_mean;
  double current_correlation;
  double covariance;
  double sdev;
  MaxCor threadMaxCor;
  threadMaxCor.cor = 0.0;
  threadMaxCor.imageID = -1;
  threadMaxCor.pixelID = -1;
  
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    current_mean = computeMean(cubes[threadSlot]->dataInterpolated[j]);
    covariance = 0.0;
    sdev = 0.0;
   
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      x = (cubes[threadSlot]->dataInterpolated[j][k]); //Extract current spectrum value
      covariance += (x - current_mean) * (ref[k] - ref_mean); 
      sdev += (x - current_mean)*(x - current_mean);
    }
    sdev = sqrt(sdev);
    
    current_correlation = covariance/(sdev * ref_sdev);
    if(current_correlation > threadMaxCor.cor)
    {
      threadMaxCor.cor = current_correlation;
      threadMaxCor.pixelID =  ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
      threadMaxCor.imageID =  ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    }
  
  }
  
  maxCorMutex.lock();
  if(threadMaxCor.cor > maxCorGlobal.cor)
  {
    maxCorGlobal = threadMaxCor;
  }
  maxCorMutex.unlock();
}

double MTInternalRef::computeMean(double *data)
{
  double mean_x = 0.0;
  for(unsigned int i = 0; i < ref.length(); i++)
  {
    mean_x += data[i];
  }
  return mean_x/((double)ref.length());
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
List CInternalReferenceSpectrum(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB,
                               Rcpp::NumericVector referenceSpectrum,  
                               Rcpp::NumericVector commonMassAxis)
{
  List out;
  try
  {
    MTInternalRef  myRef(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB, 
                      commonMassAxis,
                      referenceSpectrum);
 
    out = myRef.Run();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}
