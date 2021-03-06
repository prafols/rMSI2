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

#ifndef MT_AVERAGE_H
  #define MT_AVERAGE_H
#include <Rcpp.h>
#include <mutex>
#include "threadingmsiproc.h"

class MTAverage : public ThreadingMsiProc 
{
  public:

    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    // commonMassAxis: The common mass axis used to process and interpolate multiple datasets.
    // minTIC and maxTIC: spectra with a TIC value outside this range will not be used for average calculation
    MTAverage(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis, double minTIC, double maxTIC);
    ~MTAverage();
    
    //Execute a full imatge processing using threaded methods
    Rcpp::NumericVector Run();
    
  private:
    Rcpp::NumericVector AverageSpectrum;
    unsigned int *validPixelCount; //Count the number of pixels that meet the TIC condition in each datacube
    double TICmin, TICmax;
    double **TicNormalizations; //Copy of TIC normalization values for each image in rMSIObj_list. Need to allow a safe axes to TIC values from working threads
    std::mutex averageMutex;
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
};
#endif
