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
#include "mtcommonmass.h"
using namespace Rcpp;

MTCommonMass::MTCommonMass(Rcpp::List rMSIObj_list, int dummy_mass_axis_length, int numberOfThreads, double memoryPerThreadMB):
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, Rcpp::NumericVector(dummy_mass_axis_length), DataCubeIOMode::PEAKLIST_READ),
  bNoNeed2Resample(false)
{

}

MTCommonMass::~MTCommonMass()
{

}

List MTCommonMass::Run()
{
  Rcpp::Rcout<<"Calculating the new mass axi...\n";
  globalMergedSpc.mass.clear();
  globalMergedSpc.bins.clear();
  globalMergedSpc.bMerge = false;
  globalMergedSpc.level = 0;
  bNoNeed2Resample = true;
  
  //Run in multi-threading
  runMSIProcessingCpp();
  
  //Copy to the common mass in R object format
  NumericVector RcommonMass(globalMergedSpc.mass.size()); 
  memcpy(RcommonMass.begin(), globalMergedSpc.mass.data(), sizeof(double)*globalMergedSpc.mass.size());
  
  return List::create(Named("mass") = RcommonMass, Named("NoNeed2Resample") = bNoNeed2Resample);
}

std::vector<double> MTCommonMass::CalcMassAxisBinSize(std::vector<double> &mass, std::vector<double> &intensity)
{
  //Simple peak detector
  double slope_ant = 0, slope;
  std::vector<double> pkMass;
  std::vector<double> pkIntensity;
  std::vector<int> pkIndex;
  std::vector<double> pkBinSize;
  std::vector<int> pkIndexLeft;
  std::vector<int> pkIndexRight;
  for( int i = 0; i < (mass.size()-1); i++ )
  {
    if( i < (mass.size()-1) )
    {
      slope = intensity[i+1] - intensity[i];
    }
    if( (slope_ant >  0 ) && (slope < 0 ) && i > 0 ) 
    {
      pkMass.insert(pkMass.end(), mass[i] );
      pkIntensity.insert(pkIntensity.end(), intensity[i] );
      pkIndex.insert(pkIndex.end(), i ); //inserting C index
    }
    slope_ant = slope;
  }
  
  //Calculate the mass bin size for each peak
  for( int i = 0; i < pkIndex.size(); i++)
  {
    double localBinSize = 0;
    int count = 0;
    //Left region
    if(pkIndex[i] > 0)
    {
      for( int j = pkIndex[i]; j >= 0; j--)
      {
        if( (intensity[j] <= intensity[j-1])  )
        {
          pkIndexLeft.insert(pkIndexLeft.end(), j ); //inserting C index
          break; //End of peak
        }
        localBinSize += (mass[j] - mass[j-1]);
        count++;
      }
    }
    
    //Right region
    if(pkIndex[i] < (mass.size()-1))
    {
      for( int j = pkIndex[i]; j < (mass.size()-1); j++)
      {
        if( (intensity[j] <= intensity[j+1])  )
        {
          pkIndexRight.insert(pkIndexRight.end(), j ); //inserting C index
          break; //End of peak
        }
        localBinSize += (mass[j+1] - mass[j]);
        count++;
      }
    }
    if(count >0 )
    {
      localBinSize /= (double)count;
      pkBinSize.insert(pkBinSize.end(), localBinSize);
    }
  }
  //End of simple peak detector
  
  //Calc bin size at each input vector position
  std::vector<double> bins2(mass.size());
  int currInd = 0;
  

  for(int i=0; i < pkIndexRight.size(); i++ )
  {
    for( int  j = currInd; j <= pkIndexRight[i]; j++)
    {
      bins2[j] = pkBinSize[i];
    }
    currInd = pkIndexRight[i]; 
  }
  while(currInd < bins2.size() )
  {
    bins2[currInd] = pkBinSize[pkBinSize.size()-1];
    currInd++;
  }
  
  return bins2;
}

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
MTCommonMass::MergeTree MTCommonMass::MergeMassAxis(MTCommonMass::MergeTree &mz1, MTCommonMass::MergeTree &mz2)
{
  //Error check
  if(mz2.mass.size() == 0)
  {
    throw std::runtime_error("Error: mz2 does not contain any element.\n");
  }
  
  //Vectors concatenation and sorting
  double *newMz = new double[mz1.mass.size() + mz2.mass.size()];
  double *newBins = new double[mz1.bins.size() + mz2.bins.size()];
  int i1 = 0;
  int i2 = 0;
  int inew = 0;
  if( mz1.mass.size() == 0 )
  {
    //mz1 is empty, so just copy mz2
    for(int i = 0; i < mz2.mass.size(); i++)
    {
      newMz[inew] = mz2.mass[i];
      newBins[inew] = mz2.bins[i];
      inew++;
    }
  }
  else
  {
    //Merge mz1 and mz2
    while(true)
    {
      if(i1 >= mz1.mass.size() && i2 >= mz2.mass.size())
      {
        //We run out of elements on both mass axes
        break;
      }
      else if(i1 >= mz1.mass.size())
      {
        //We run out of elements in mz1
        if(inew == 0)
        {
          newMz[inew] = mz2.mass[i2];
          inew++;
        }
        else if( (mz2.mass[i2] -  newMz[inew - 1]) >= mz2.bins[i2] )
        {
          newMz[inew] = mz2.mass[i2];
          newBins[inew] = mz2.bins[i2];
          inew++;
        }
        i2++;
      }
      else if(i2 >= mz2.mass.size())
      {
        //We run out of elements in mz2
        if(inew == 0)
        {
          newMz[inew] = mz1.mass[i1];
          inew++;
        }
        else if( (mz1.mass[i1] -  newMz[inew - 1]) >= mz1.bins[i1] )
        {
          newMz[inew] = mz1.mass[i1];
          newBins[inew] = mz1.bins[i1];
          inew++;
        }
        i1++;
      }
      else
      {
        //There are remaining elements in both mass axes
        if(mz1.mass[i1] < mz2.mass[i2])
        {
          if(inew == 0)
          {
            newMz[inew] = mz1.mass[i1];
            inew++;
          }
          else if( (mz1.mass[i1] -  newMz[inew - 1]) >= mz1.bins[i1] )
          {
            newMz[inew] = mz1.mass[i1];
            newBins[inew] = mz1.bins[i1];
            inew++;
          }
          i1++;
        }
        else
        {
          if(inew == 0)
          {
            newMz[inew] = mz2.mass[i2];
            inew++;
          }
          else if( (mz2.mass[i2] -  newMz[inew - 1]) >= mz2.bins[i2] )
          {
            newMz[inew] = mz2.mass[i2];
            newBins[inew] = mz2.bins[i2];
            inew++;
          }
          i2++;
        }
      }
    }
  }
  
  //copy used elements to the end result
  MergeTree resultMass;
  resultMass.mass.resize(inew);
  resultMass.bins.resize(inew);
  resultMass.bMerge = mz2.bMerge;
  resultMass.level = mz1.level + 1;
  
  memcpy(resultMass.mass.data(), newMz,  sizeof(double)*inew);
  memcpy(resultMass.bins.data(), newBins,  sizeof(double)*inew);
  
  delete[] newMz;
  delete[] newBins;
  
  return resultMass;
}

void MTCommonMass::ProcessingFunction(int threadSlot)
{
  bool thread_bNoNeed2Resample = true;
  unsigned int icurr = 0;
  bool bLoad = false;
  bool bMassMerge;
  std::vector<MergeTree> threadMergedSpc;
  PeakPicking::Peaks *mpeaks; //Pointer to the current peaklist
  
  //Read only the first mass axis to compare if others are identical, this is the case for Bruker FTICR
  mpeaks = cubes[threadSlot]->peakLists[0];
  
  //TODO I implemented this for Bruker's FTICR data which has some duplicates in the mass axis, now it is commented out... think about it and implement in the imzMLbinreader if needed.
  //Fix duplicates and zero drops
  //FirstSpectrumFixed <- fixImzMLDuplicates(mzdd, dd)
          
  //Calculate first mass axis bin size to avoid having to calculate it at each iteration
  PeakPicking::Peaks firstSpectrum;
  firstSpectrum.mass = mpeaks->mass;
  firstSpectrum.intensity = mpeaks->intensity;
  firstSpectrum.binSize = CalcMassAxisBinSize( firstSpectrum.mass, firstSpectrum.intensity );
            
  while(true)
  {
    if(bLoad)
    {
      //Read the current spectrum 
      mpeaks = cubes[threadSlot]->peakLists[icurr];
          
      //Get Bin size at peaks
      if(firstSpectrum.mass == mpeaks->mass)
      {
        bMassMerge = false;
        mpeaks->binSize = firstSpectrum.binSize;
      }
      else
      {
        bMassMerge = true;
        mpeaks->binSize = CalcMassAxisBinSize( mpeaks->mass, mpeaks->intensity);
      }
      
      thread_bNoNeed2Resample &= (!bMassMerge);
              
      //Shift register and fill the fist element
      threadMergedSpc.emplace(threadMergedSpc.begin()); //Add the new element at the begining so previous elements are moved to the rigth
      threadMergedSpc[0].level = 0;
      threadMergedSpc[0].mass = mpeaks->mass;
      threadMergedSpc[0].bins = mpeaks->binSize;
      threadMergedSpc[0].bMerge = bMassMerge;
      
      icurr++;
      bLoad = false;
    } //end if bload
      
      if( threadMergedSpc.size( ) > 1 )
      {
        if(threadMergedSpc[0].level == threadMergedSpc[1].level || icurr >= cubes[threadSlot]->nrows )
        { 
          if( threadMergedSpc[1].bMerge)
          {
            //Merge!
            threadMergedSpc[0] = MergeMassAxis(threadMergedSpc[0], threadMergedSpc[1]); 
          }
          else
          {
            //Both mass axes are identical so there is no need to merge them, this is the case for Bruker FTICR data
            threadMergedSpc[0].bMerge = threadMergedSpc[1].bMerge;
            threadMergedSpc[0].level += 1;
          }
          
          //Shift register
          threadMergedSpc.erase(threadMergedSpc.begin() + 1); //Just remove the second element
          
          //End Condition
          if(icurr >= cubes[threadSlot]->nrows && threadMergedSpc.size() == 1)
          {
            break;
          }
        }
        else
        {
          bLoad = true;
        }
      }
      else
      {
        bLoad = true;
      }
  } //while end
    
  //Update the global mass axis using the threadMergedSpc[0].mass which contains the cube merged mass axis
  commonMutex.lock();
  globalMergedSpc = MergeMassAxis(globalMergedSpc, threadMergedSpc[0]);
  bNoNeed2Resample &= thread_bNoNeed2Resample; 
  commonMutex.unlock();
}

// Calculate the common mass axis for imzML data in processed mode.
// [[Rcpp::export]]
List CcommonMassAxis(List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB)
{
  List out;
  unsigned long long mean_number_mass_channels = 0;
  unsigned int total_num_of_pixels = 0;
  List data;
  List imzML;
  DataFrame imzMLrun;
  IntegerVector massLengths;
  
  try
  {
    
    for( auto rMSIobj : rMSIObj_list) //For each rMSI object in the list
    {
      data = (as<List>(rMSIobj))["data"];
      imzML = data["peaklist"]; //Data must be set as a peaklist even it is not a peak list!
      
      //Check if data is in processed mode
      if(as<bool>(imzML["continuous_mode"]))
      {
        throw std::runtime_error("Error: imzML must be in processed mode to calculate a common mass axis.\n");
      }
      
      imzMLrun = as<DataFrame>(imzML["run_data"]);
      massLengths = imzMLrun["mzLength"];
      for( int item : massLengths)
      {
        mean_number_mass_channels += item;
        total_num_of_pixels++;
      }
    }
  
    mean_number_mass_channels /= total_num_of_pixels;
    
    //Apply a lower limit, assume at least 100 mass channelks
    mean_number_mass_channels = mean_number_mass_channels < 100 ? 100 : mean_number_mass_channels;
    
    
    MTCommonMass myCommMass(rMSIObj_list, 
                              mean_number_mass_channels,
                              numOfThreads, 
                              memoryPerThreadMB);
   
    out = myCommMass.Run();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}
