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

#ifndef MT_NORMALIZATIONMEANSPECTRA_H
  #define MT_NORMALIZATIONMEANSPECTRA_H
#include <Rcpp.h>
#include <vector>
#include "threadingmsiproc.h"

class MTNormalizationMeanSpectra : public ThreadingMsiProc 
{
  public:

    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    // commonMassAxis: The common mass axis used to process and interpolate multiple datasets.
    MTNormalizationMeanSpectra(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis);
    ~MTNormalizationMeanSpectra();
    
    //Execute a full imatge processing using threaded methods
    Rcpp::List Run();
    
  private:
    
    typedef struct
    {
      double TIC;
      double RMS;
      double MAX;
    }PixelNorms;
    
    std::vector<std::vector<PixelNorms>> Normalizations;
    
    std::mutex mutex_copyData; //Mutex to avoid coping data to lsNorms (R object) from various threads simultaneously
    
    std::vector<Rcpp::NumericVector> averageSpectrum;
    std::vector<Rcpp::NumericVector> baseSpectrum;
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
    Rcpp::List rMSIObj_lst; //Copy of the rMSI object
    std::vector<unsigned int> num_of_pixels; //Number of pixel in each rMSI object
};
#endif
