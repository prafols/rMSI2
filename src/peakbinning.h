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

#ifndef PEAKBINNING_H
#define PEAKBINNING_H
#include <Rcpp.h>
#include <mutex>
#include "threadingmsiproc.h"
#include "peakpicking.h"

class PeakBinning : public ThreadingMsiProc 
{
public:
  
  //Constructor arguments:
  // rMSIObj_list: A list of rMSI objects to process
  // numberOfThreads: Total number of threads to use during processing
  // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
  // preProcessingParams: An R reference class with the pre-processing parameters.
  PeakBinning(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, 
              Rcpp::Reference preProcessingParams);
  
  ~PeakBinning();
  
  //Perform the peak binning, returns a binned peak matrix with all the mass channels in place ready to be filled by the peak-fill algorithm
  Rcpp::List Run(); 

private:
  double binFilter;
  double tolerance;
  bool tolerance_in_ppm; //If true the binning tolerance is specified in  ppm, if false then the number of datapoints per peak is used instead
  int totalNumOfPixels;
  
  //Mass bin data structure
  typedef struct
  {
    double mass; //The mass channel for the mass bin
    double binSize; //The bin size for thes mass bin
    unsigned int counts; //Number of counts (pixels) in the mass bin
  }MassBin;
  
  std::vector<MassBin> mainMassBins; //The main mass bins object
  std::mutex mainBinsMutex; //Mutex to lock the main mass bins object for safe multithreaded editing
  
  //Thread Processing function definition
  void ProcessingFunction(int threadSlot);
  
  //newMassBin: the new mass channel to add in the peak matrix.
  void AppendMassChannel2MassBins(MassBin newMassBin, std::vector<MassBin> &targetMassBins);
};
#endif
