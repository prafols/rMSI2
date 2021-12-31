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
#include "peakbinning.h"
#include "progressbar.h"
#include <limits>       // std::numeric_limits
#include <stdexcept>
using namespace Rcpp;

PeakBinning::PeakBinning(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, 
                         Rcpp::Reference preProcessingParams):
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, Rcpp::NumericVector(), DataCubeIOMode::PEAKLIST_READ)
{
  //Calculate the total number of pixels
  totalNumOfPixels = 0;
  for(int i=0; i < ioObj->getNumberOfCubes(); i++ )
  {
    totalNumOfPixels += ioObj->getNumberOfPixelsInCube(i);
  }
  
  //Get the parameters
  Rcpp::Reference binningParams = preProcessingParams.field("peakbinning");
  tolerance = binningParams.field("tolerance");
  tolerance_in_ppm = binningParams.field("tolerance_in_ppm");
  binFilter = binningParams.field("binFilter"); 
}

PeakBinning::~PeakBinning()
{

}

void PeakBinning::AppendMassChannel2MassBins(MassBin newMassBin, std::vector<MassBin> &targetMassBins)
{
  //Search the closest mass bin
  double minMassDistance = std::numeric_limits<double>::max();
  int minDistanceIndex = -1;
  double currentMassDistance;
  for(int icol = 0; icol < targetMassBins.size(); icol++)
  {
    currentMassDistance = newMassBin.mass - targetMassBins[icol].mass;
    
    if(fabs(currentMassDistance) < fabs(minMassDistance))
    {
      minMassDistance = currentMassDistance;
      minDistanceIndex = icol;
    }
    
    if(currentMassDistance < 0)
    {
      break;
    }
  }
  
  //Select the kind of binning tolerance
  double compTolerance;
  if(tolerance_in_ppm)
  {
    minMassDistance = 1e6*(fabs(minMassDistance)/newMassBin.mass); //Compute distance in ppm
    compTolerance = tolerance;
  }
  else
  {
    compTolerance = tolerance * newMassBin.binSize;
  } 
  
  if(minDistanceIndex >= 0)
  {
    if( (fabs(minMassDistance) < compTolerance))
    {
      //The peak must be binned with the minDistanceIndex column of the peak matrix
      targetMassBins[minDistanceIndex].mass = 
        ((double)targetMassBins[minDistanceIndex].counts)/((double)(targetMassBins[minDistanceIndex].counts + 1)) * targetMassBins[minDistanceIndex].mass 
        + newMassBin.mass/((double)(targetMassBins[minDistanceIndex].counts + 1));
      
      targetMassBins[minDistanceIndex].binSize = 
      ((double)targetMassBins[minDistanceIndex].counts)/((double)(targetMassBins[minDistanceIndex].counts + 1)) * targetMassBins[minDistanceIndex].binSize 
        + newMassBin.binSize/((double)(targetMassBins[minDistanceIndex].counts + 1));
      
      targetMassBins[minDistanceIndex].counts+=newMassBin.counts;
    }
    else
    {
      if(minMassDistance > 0)
      {
        //Insert the new mass channel after the minDistanceIndex
        targetMassBins.insert(targetMassBins.begin() + minDistanceIndex + 1, newMassBin);
      }
      else
      {
        //Insert the new mass channel before the minDistanceIndex
        targetMassBins.insert(targetMassBins.begin() + minDistanceIndex, newMassBin);
      }
    }
  }
  else
  {
    //targetMassBins is empty, so add the first element to it
    targetMassBins.push_back(newMassBin);
  }
}

void PeakBinning::ProcessingFunction(int threadSlot)
{
  std::vector<MassBin> thread_binMass; //The mass name for each matrix column (local thread space)
  MassBin current_bin;

  //Thread local worker
  PeakPicking::Peaks *mpeaks; //Pointer to the current peaklist
  for( int irow = 0; irow < cubes[threadSlot]->nrows; irow++)
  {
    mpeaks = cubes[threadSlot]->peakLists[irow];
    for(int ipeak = 0; ipeak < mpeaks->mass.size(); ipeak++)
    {
      current_bin.mass = mpeaks->mass[ipeak];
      current_bin.binSize = mpeaks->binSize[ipeak];
      current_bin.counts = 1;
      AppendMassChannel2MassBins(current_bin, thread_binMass);
    }
  }
  
  //Append local thread result to the main bins
  if(!thread_binMass.empty())
  {
    mainBinsMutex.lock();
    for( auto it = thread_binMass.begin(); it != thread_binMass.end(); ++it)
    {
      AppendMassChannel2MassBins(*it, mainMassBins);
    }
    mainBinsMutex.unlock();
  }
}


List PeakBinning::Run()
{
  //Run the mass binning in multithreading
  mainMassBins.clear();
  Rcout<<"Binning peaks...\n";
  runMSIProcessingCpp(); //After the multithreaded binning the mainMassBins object contains the sorted mass channels and the counts on each
  
  //Apply binFilter
  NumericVector massR;
  NumericVector binSizeR;
  for( auto it = mainMassBins.begin(); it != mainMassBins.end(); ++it)
  {
    if( ((double)(it->counts)) > binFilter*(double)totalNumOfPixels)
    {
      massR.push_back(it->mass);
      binSizeR.push_back(it->binSize);
    }
  }
  if(massR.length() == 0)
  {
    throw std::runtime_error("ERROR: all peaks were removed by the bin filter. Consider decresin the bin filter or SNR.\n");
  }
  Rcout<<"Bining complete with a total number of "<<massR.length()<<" bins\n"; 
  
  //Prepare the R bin matrix
  NumericMatrix binMatIntensity(totalNumOfPixels, massR.length());
  NumericMatrix binMatSNR(totalNumOfPixels, massR.length());
  NumericMatrix binMatArea(totalNumOfPixels, massR.length());
  
  return List::create( Named("mass") = massR, Named("binSize") = binSizeR, Named("intensity") = binMatIntensity, Named("SNR") = binMatSNR, Named("area") = binMatArea );
}

// [[Rcpp::export]]
List CRunPeakBinning(Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                     Rcpp::Reference preProcessingParams)
{
  List out;
  try
  {
    PeakBinning myPeakBinning(rMSIObj_list, numOfThreads, memoryPerThreadMB, 
                              preProcessingParams);
  
    out = myPeakBinning.Run(); 
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}

