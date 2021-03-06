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

#include "threadingmsiproc.h" 
#include "progressbar.h"
#include <stdexcept>

//#define __DEBUG__ // comment out in the final release!

ThreadingMsiProc::ThreadingMsiProc(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis,
                                   DataCubeIOMode storeDataModeimzml, Rcpp::StringVector uuid, Rcpp::String outputImzMLPath, Rcpp::StringVector outputImzMLfnames):
  dataStoreMode(storeDataModeimzml), massAxis(commonMassAxis)
{
  if( (dataStoreMode == DataCubeIOMode::DATA_STORE) || (dataStoreMode == DataCubeIOMode::PEAKLIST_STORE) )
  {
    if(uuid.size() != rMSIObj_list.length())
    {
      throw std::runtime_error("Error: uuid vector must have the same length as rMSIObj_list\n");
    }
    
    if(outputImzMLPath == "")
    {
      throw std::runtime_error("Error: outputImzMLPath is empty\n");
    }
      
    if(outputImzMLfnames.size() != rMSIObj_list.length())
    {
      throw std::runtime_error("Error: outputImzMLfname vector must have the same length as rMSIObj_list\n");
    }
  }
  
  //Force dataStoreMode to peak list read mode if no spectral data is available
  if( dataStoreMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ && massAxis.length() == 0)
  {
    //common mass axis has a length of zero, so it is assumed no spectral data available
    dataStoreMode = DataCubeIOMode::PEAKLIST_READ;
  }
  
  numOfThreadsDouble = 2*numberOfThreads;
  ioObj = new CrMSIDataCubeIO( massAxis, memoryPerThreadMB, dataStoreMode, outputImzMLPath);
  
  //Call the append method for each image in the list
  for(int i = 0; i < rMSIObj_list.length(); i++)
  {
    if( (dataStoreMode == DataCubeIOMode::DATA_STORE) || (dataStoreMode == DataCubeIOMode::PEAKLIST_STORE) )
    {
      Rcpp::String RcppStr_uuid = uuid[i];
      Rcpp::String RcppStr_outImzmls = outputImzMLfnames[i];
      ioObj->appedImageData(rMSIObj_list[i], RcppStr_uuid.get_cstring(), RcppStr_outImzmls.get_cstring());
    }
    else
    {
      ioObj->appedImageData(rMSIObj_list[i]);
    }
  }
  
  cubes = new CrMSIDataCubeIO::DataCube*[numOfThreadsDouble];
  iCube = new int[numOfThreadsDouble];
  bDataReady = new bool[numOfThreadsDouble];
  bRunningThread = new bool[numOfThreadsDouble];
  tworkers = new std::thread[numOfThreadsDouble]; //There will be double of thread objects than the actually running threads
  
  numPixels = 0;
  for (int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
    numPixels += ioObj->getNumberOfPixelsInCube(i);
  }
  
}

ThreadingMsiProc::~ThreadingMsiProc()
{
  delete[] cubes;
  delete[] iCube;
  delete[] bDataReady;
  delete[] bRunningThread;
  delete[] tworkers; 
  delete ioObj;
}

void ThreadingMsiProc::runMSIProcessingCpp()
{
  //Initialize processing data cube
  life_end = false;
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    iCube[i] = -1; //-1 means that there is no any cube assigned to worker thread
    bDataReady[i] = false;
    bRunningThread[i] = false;
  }
  
  int nextCubeLoad = 0; //Point to the next datacube to load
  int nextCubeStore = 0; //Point to the next datacube to store
  int runningThreads = 0; //Total number of running threads
  bool end_of_program = false;
  
  while( !end_of_program ) 
  {
    //Load data for future threads and start working threads
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if(iCube[iThread] == -1 && nextCubeLoad < ioObj->getNumberOfCubes()) //No cube assigned then no thread running in this slot
      {
        progressBar(nextCubeLoad, ioObj->getNumberOfCubes(), "=", " ");
        iCube[iThread] = nextCubeLoad;
        nextCubeLoad++;
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Loading data for cube = " << iCube[iThread] << " on thread " << iThread << "\n";
#endif
        cubes[iThread] = ioObj->loadDataCube(iCube[iThread]);
        if(2*runningThreads < numOfThreadsDouble)
        {
          bRunningThread[iThread] = true;
          tworkers[iThread] = std::thread(std::bind(&ThreadingMsiProc::ProcessingThread, this, iThread)); //Start Thread 
          runningThreads++;
        }
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Started work for cube = " << iCube[iThread] << " on thread " << iThread << "\n";
#endif
      }
    }

    //Wait for thread ends
#ifdef __DEBUG__
    Rcpp::Rcout << "DBG: waiting for some thread...\n"; 
#endif
    if(runningThreads > 0 )
    {
      WaitForSomeThreadEnd();
    }
#ifdef __DEBUG__
    Rcpp::Rcout << "DBG: waiting for data mutex...\n";
#endif
    mtx.lock(); //Any thread locked to this mutex is actually waiting to save data
#ifdef __DEBUG__
    Rcpp::Rcout << "DBG: got the data mutex\n";
#endif
    
    //Update number of running threads
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      runningThreads = 0;
      if(bRunningThread[iThread])
      {
        runningThreads++;
      }
    }

    //Start thread that have data loaded ready to be processed
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if(iCube[iThread] != -1 && !bRunningThread[iThread] && !bDataReady[iThread] && 2*runningThreads < numOfThreadsDouble)
      {
        bRunningThread[iThread] = true;
        tworkers[iThread] = std::thread(std::bind(&ThreadingMsiProc::ProcessingThread, this, iThread)); //Start Thread 
        runningThreads++;
      }      
    }

    //Save data and free thread slots
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      if(bDataReady[iThread] && 
          ((nextCubeStore == iCube[iThread]) ||
          (dataStoreMode == DataCubeIOMode::DATA_READ) ||
          (dataStoreMode == DataCubeIOMode::PEAKLIST_READ) || 
          (dataStoreMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ) ) ) 
      {
        //If destination imzML is set then store the data
        if( (dataStoreMode == DataCubeIOMode::DATA_STORE) || (dataStoreMode == DataCubeIOMode::PEAKLIST_STORE) )
        {
#ifdef __DEBUG__
          Rcpp::Rcout << "DBG: Storing cube = " << cubes[iThread]->cubeID << " on thread " << iThread << "\n";
#endif
          ioObj->storeDataCube(cubes[iThread]);
          nextCubeStore++;
          
#ifdef __DEBUG__
          Rcpp::Rcout << "DBG: Storeing completed\n";
#endif
        }
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Frree cube = " << cubes[iThread]->cubeID << " on thread " << iThread << "\n";
#endif
        ioObj->freeDataCube(cubes[iThread]);
        iCube[iThread] = -1; //Mark thread as stopped
        bDataReady[iThread] = false; //Reset data ready state;
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Free completed, joning thread...\n";
#endif
        tworkers[iThread].join();
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Join completed\n";
#endif
      } 
    }

#ifdef __DEBUG__
    Rcpp::Rcout <<"\nDBG Thread report:\n"; 
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      Rcpp::Rcout <<"Thread "<< iThread <<" : iCube = " << iCube[iThread] << " : running = " << bRunningThread[iThread] << " : dataReady = " << bDataReady[iThread]<<"\n";
    }
    Rcpp::Rcout <<"\n\n\n";
#endif
    
    //Check end condition
    if( nextCubeLoad >= ioObj->getNumberOfCubes() )
    {
      end_of_program = true;
      for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
      {
        end_of_program &= (iCube[iThread] == -1);
      }
    }
    mtx.unlock();
  }
#ifdef __DEBUG__
  Rcpp::Rcout << "DBG: MT proc END\n";
#endif
  Rcpp::Rcout<<"\n";
}

void ThreadingMsiProc::ProcessingThread( int threadSlot )
{
  ioObj->interpolateDataCube(cubes[threadSlot]); 
  
  //Call the processing function for this thread
  ProcessingFunction(threadSlot);
  
  //Save data to R Session and Store the new state of total processed cubes
  mtx.lock();
  bDataReady[threadSlot] = true;
  bRunningThread[threadSlot] = false;
  mtx.unlock();
  
  //Notify Main thread
  std::unique_lock<std::mutex> lock(life_end_mutex);
  if(!life_end) //Notify only if it haven been done already by another thread
  {
    life_end = true;
    life_end_cond.notify_all();
  }
}

void ThreadingMsiProc::WaitForSomeThreadEnd()
{
  std::unique_lock<std::mutex> lock(life_end_mutex);
  while (!life_end)
  {
    life_end_cond.wait(lock);
  }
  life_end =false; //Reset it for the next thread
}

void ThreadingMsiProc::ProcessingFunction(int threadSlot)
{
  Rcpp::Rcout<<"Function ThreadingMsiProc::ProcessingFunction(int threadSlot) has been called from thread slot: "<<threadSlot<<"\n";
  Rcpp::Rcout<<"The base class ThreadingMsiProc can not be used directly and must be derived reimplementing the ProcessingFunction()\n";
}
