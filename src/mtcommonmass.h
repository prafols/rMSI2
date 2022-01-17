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

#ifndef MT_COMMONMASS_H
  #define MT_COMMONMASS_H
#include <Rcpp.h>
#include <vector>
#include <mutex>
#include "threadingmsiproc.h"

class MTCommonMass : public ThreadingMsiProc 
{
  public:

    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // dummy_mass_axis_length: an integer with the mean number of mass channels in all the peak lists. It is used to calculate the memory data cube for each thread.
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    MTCommonMass(Rcpp::List rMSIObj_list, int dummy_mass_axis_length, int numberOfThreads, double memoryPerThreadMB);
    ~MTCommonMass();
    
    //Execute a full imatge processing using threaded methods
    //Returns a List with the common mass axis and a boolean idicating if resampling is needed (should be true for FTICR data and false for Orbitrap)
    Rcpp::List Run();
    
  private:
    std::mutex commonMutex;
    bool bNoNeed2Resample;
    
    typedef struct
    {
      unsigned int level;
      bool bMerge;
      std::vector<double> mass;
      std::vector<double> bins;
    }MergeTree;
    
    //Object containging the global mass axis
    MergeTree globalMergedSpc;
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
    //' CalcMassAxisBinSize.
    //' 
    //' Calc the bin size of a mass axis at each mass channels using simple peak-picking information.
    //' 
    //' @param mass the mass axis.
    //' @param intensity the intensity of a given spectrum.
    //' 
    //' @return the bin size of each m/z channel.
    std::vector<double> CalcMassAxisBinSize(std::vector<double> &mass, std::vector<double> &intensity);
    
    //' MergeMassAxis.
    //' 
    //' Merges two mass axis in a single one using an apropiate bin size.
    //' The resulting mass axis will display a bin size equal to the minimum of two supplied vectors. 
    //' The bin size must be supplied along each input mass axis.
    //' The first mass axis (mz1) can be a zero-length vector.
    //' 
    //' @param mz1 the first mass axis to merge.
    //' @param mz2 the second mass axis to merge.
    //' 
    //' @return the common mass axis that represents mz1 and mz1 accurately.
    //' 
    MergeTree MergeMassAxis(MergeTree &mz1, MergeTree &mz2);
    
};
#endif
