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

#ifndef MT_INTERNAL_REF_H
  #define MT_INTERNAL_REF_H
#include <Rcpp.h>
#include <mutex>
#include "threadingmsiproc.h"

class MTInternalRef : public ThreadingMsiProc 
{
  public:

    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    // commonMassAxis: The common mass axis used to process and interpolate multiple datasets.
    // reference: The reference spectrum
    MTInternalRef(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis, Rcpp::NumericVector reference);
    ~MTInternalRef();
    
    //Execute a full imatge processing using threaded methods
    Rcpp::List Run();
    
  private:
    
    Rcpp::NumericVector ref; //Local copy of the reference
    double ref_sdev; //The standard deviation of the reference spectrum
    double ref_mean; //The mean of the reference spectrum
    
    typedef struct{
      int pixelID;
      int imageID;
      double cor;
    }MaxCor;
    
    std::mutex maxCorMutex;
    MaxCor maxCorGlobal; //Will contain the finally selected pixel
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
    //Methods to compute the mean, number of elements is the same as the global mass axis
    double computeMean(double *data);
    
};
#endif
