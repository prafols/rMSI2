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
  slot_mutexes = new std::mutex[numOfThreadsDouble];

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
  delete[] slot_mutexes;
  delete ioObj;
}

void ThreadingMsiProc::runMSIProcessingCpp()
{
  // Initialize processing data cube
  {
    std::lock_guard<std::mutex> lock(life_end_mutex);
    life_end = false;
  }
  
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    std::lock_guard<std::mutex> lock(slot_mutexes[i]);
    iCube[i] = -1; //-1 means that there is no any cube assigned to worker thread
    bDataReady[i] = false;
    bRunningThread[i] = false;
  }
  
  int nextCubeLoad = 0; // Point to the next datacube to load
  int nextCubeStore = 0; // Point to the next datacube to store
  int runningThreads = 0; // Total number of running threads (managed exclusively by main thread)
  int totalCubes = ioObj->getNumberOfCubes();
  bool end_of_program = false;
  
  while(!end_of_program) 
  {
    // 1. Load data and start working threads strictly respecting the concurrency limit
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      slot_mutexes[iThread].lock();
      bool can_load_and_run = (iCube[iThread] == -1 && nextCubeLoad < totalCubes && (2 * runningThreads < numOfThreadsDouble));
      int currentCubeToLoad = -1;
      
      if(can_load_and_run) 
      {
        currentCubeToLoad = nextCubeLoad;
        iCube[iThread] = nextCubeLoad;
        nextCubeLoad++;
        bRunningThread[iThread] = true;
      }
      slot_mutexes[iThread].unlock();
      
      if(can_load_and_run)
      {
        progressBar(currentCubeToLoad, totalCubes, "=", " ");
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Loading data for cube = " << currentCubeToLoad << " on thread slot " << iThread << "\n";
#endif
        cubes[iThread] = ioObj->loadDataCube(currentCubeToLoad);
        
        runningThreads++;
        tworkers[iThread] = std::thread(&ThreadingMsiProc::ProcessingThread, this, iThread); // Start Thread
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Started work for cube = " << currentCubeToLoad << " on thread slot " << iThread << "\n";
#endif
      }
    }
    
    // Check if there are actually any threads still running before waiting
    bool any_thread_running = false;
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      slot_mutexes[iThread].lock();
      if(bRunningThread[iThread]) 
      {
        any_thread_running = true;
      }
      slot_mutexes[iThread].unlock();
      if(any_thread_running) break;
    }
    
    //2. Wait for thread ends only if there is actively something running
#ifdef __DEBUG__
    Rcpp::Rcout << "DBG: waiting for some thread... any_running = " << any_thread_running << "\n"; 
#endif
    
    if(runningThreads > 0 && any_thread_running)
    {
      WaitForSomeThreadEnd();
    }
    
#ifdef __DEBUG__
    Rcpp::Rcout << "DBG: thread notification received, processing slots...\n";
#endif
    
    // 3. Save data and free thread slots
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      slot_mutexes[iThread].lock();
      bool readyToStore = bDataReady[iThread];
      int cubeIdToCheck = iCube[iThread];
      slot_mutexes[iThread].unlock();
      
      if(readyToStore && 
         ((nextCubeStore == cubeIdToCheck) ||
         (dataStoreMode == DataCubeIOMode::DATA_READ) ||
         (dataStoreMode == DataCubeIOMode::PEAKLIST_READ) || 
         (dataStoreMode == DataCubeIOMode::DATA_AND_PEAKLIST_READ) ) ) 
      {
        // If destination imzML is set then store the data
        if((dataStoreMode == DataCubeIOMode::DATA_STORE) || (dataStoreMode == DataCubeIOMode::PEAKLIST_STORE))
        {
#ifdef __DEBUG__
          Rcpp::Rcout << "DBG: Storing cube = " << cubes[iThread]->cubeID << " on thread slot " << iThread << "\n";
#endif
          ioObj->storeDataCube(cubes[iThread]);
          nextCubeStore++;
          
#ifdef __DEBUG__
          Rcpp::Rcout << "DBG: Storing completed\n";
#endif
        }
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Freeing cube = " << cubes[iThread]->cubeID << " on thread slot " << iThread << "\n";
#endif
        ioObj->freeDataCube(cubes[iThread]);
        
        slot_mutexes[iThread].lock();
        iCube[iThread] = -1; // Mark thread as stopped
        bDataReady[iThread] = false; // Reset data ready state
        bRunningThread[iThread] = false; // Ensure running flag is cleared
        slot_mutexes[iThread].unlock();
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Free completed, joining thread...\n";
#endif
        
        if(tworkers[iThread].joinable())
        {
          tworkers[iThread].join();
        }
        
        runningThreads--; // Decrement active running threads count safely
        
        // Signal thread completion
        ThreadCompleteCallback(iThread);
        
#ifdef __DEBUG__
        Rcpp::Rcout << "DBG: Join completed\n";
#endif
      } 
    }
    
#ifdef __DEBUG__
    Rcpp::Rcout << "\nDBG Thread report:\n"; 
    for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
    {
      slot_mutexes[iThread].lock();
      Rcpp::Rcout << "Thread " << iThread << " : iCube = " << iCube[iThread] << " : running = " << bRunningThread[iThread] << " : dataReady = " << bDataReady[iThread] << "\n";
      slot_mutexes[iThread].unlock();
    }
    Rcpp::Rcout << "\n\n\n";
#endif
    
    // Check end condition
    if(nextCubeLoad >= totalCubes)
    {
      end_of_program = true;
      for(int iThread = 0; iThread < numOfThreadsDouble; iThread++)
      {
        slot_mutexes[iThread].lock();
        if(iCube[iThread] != -1) 
        {
          end_of_program = false;
        }
        slot_mutexes[iThread].unlock();
      }
    }
  }
#ifdef __DEBUG__
  Rcpp::Rcout << "DBG: MT proc END\n";
#endif
  Rcpp::Rcout << "\n";
}

void ThreadingMsiProc::ProcessingThread(int threadSlot)
{
  
  ioObj->interpolateDataCube(cubes[threadSlot]); 
  
  // Execute the actual data processing task for this cube
  ProcessingFunction(threadSlot);
  
  // Safely mark data as ready and clear the running flag for this specific slot
  {
    std::lock_guard<std::mutex> lock(slot_mutexes[threadSlot]);
    bDataReady[threadSlot] = true;
    bRunningThread[threadSlot] = false;
  } 
  
  // Signal the main thread that work has been completed
  {
    std::lock_guard<std::mutex> lock(life_end_mutex);
    life_end = true;
  }
  life_end_cond.notify_one();
}

void ThreadingMsiProc::WaitForSomeThreadEnd()
{
  std::unique_lock<std::mutex> lock(life_end_mutex);
  life_end_cond.wait(lock, [this] { return life_end; });
  life_end = false; // Reseteamos la bandera para el siguiente evento
}

void ThreadingMsiProc::ProcessingFunction(int threadSlot)
{
  Rcpp::Rcout<<"Function ThreadingMsiProc::ProcessingFunction(int threadSlot) has been called from thread slot: "<<threadSlot<<"\n";
  Rcpp::Rcout<<"The base class ThreadingMsiProc can not be used directly and must be derived reimplementing the ProcessingFunction()\n";
}

void ThreadingMsiProc::ThreadCompleteCallback(int threadSlot)
{
  //Do nothing by default. Derived classes may overrride this to execute code post-thread completion
}
