/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2019 Pere Rafols Soler
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

#include <future>

#include "fastbinning.h"
#include "progressbar.h"

using namespace Rcpp;

RapidBinning::RapidBinning(List rMSIobject, NumericVector peak_masses, NumericVector peak_widths, NumericMatrix PeakMatrix, unsigned int nThreads, unsigned int PeakMatrixRowOffset):
  img(rMSIobject),
  pmasses(peak_masses),
  pwidths(peak_widths),
  pMat(PeakMatrix),
  n_threads(nThreads),
  img_row_offset(PeakMatrixRowOffset)
{
  n_peaks = pmasses.length();
  if(n_peaks != pwidths.length())
  {
    Rcpp::stop("Dimension mismatch: 'pmasses' (length %d) and 'peak_widths' (length %d) must have the same length.", 
               n_peaks, peak_widths.length());
  }
  
  //Init rMSIXBin objects
  vec_XBin.reserve(n_threads);
  for (unsigned int t = 0; t < n_threads; ++t)
  {
    vec_XBin.push_back(std::make_unique<rMSIXBin>(rMSIobject, 1));
  }
  
  //Get the mass axis
  massAxis = as<NumericVector>(rMSIobject["mass"]);
  
  // Load coordinate matrix converting from R 1-based to C++ 0-based indexing
  NumericMatrix XYCoords = img["pos"];
  size_t num_pixels = XYCoords.nrow();
  
  vec_XYCoords.reserve(num_pixels);
  
  for (size_t i = 0; i < num_pixels; ++i)
    {
    vec_XYCoords.emplace_back(XYPoint{
      static_cast<unsigned int>(XYCoords(i, 0) - 1),
      static_cast<unsigned int>(XYCoords(i, 1) - 1)
    });
  }

}

RapidBinning::~RapidBinning()
{
  
}

void RapidBinning::run()
{
  Rcout << "Rapid Peak Binning..." << std::endl;
  
  std::vector<TaskHandle> active_tasks;
  active_tasks.reserve(n_threads);
  
  std::vector<unsigned int> free_thread_ids;
  free_thread_ids.reserve(n_threads);
  for (unsigned int t = 0; t < n_threads; ++t)
  {
    free_thread_ids.push_back(t);
  }
  
  //Mass search constants
  size_t last_mass_index = 0;
  
  for (unsigned int j = 0; j < n_peaks; ++j)
  {
    // If we hit our thread budget, poll/wait until AT LEAST ONE slot opens up
    while (active_tasks.size() >= n_threads)
    {
      for (auto it = active_tasks.begin(); it != active_tasks.end(); )
      {
        // Check if thread finished without blocking indefinitely
        if (it->future.wait_for(std::chrono::microseconds(100)) == std::future_status::ready)
        {
          it->future.get(); // Propagate exceptions if any occurred
          
          // Return freed thread ID back to the available pool
          free_thread_ids.push_back(it->thread_id);
          
          // Erase completed task handle
          it = active_tasks.erase(it);
          break; // Free slot found! Exit check loop
        }
        else
        {
          ++it;
        }
      }
      
      // Prevent excessive CPU polling if all workers are currently busy
      if (active_tasks.size() >= n_threads) {
        std::this_thread::sleep_for(std::chrono::microseconds(100));
      }
      
    }
    
    // Pop an available thread ID from the free pool
    unsigned int assigned_thread_id = free_thread_ids.back();
    free_thread_ids.pop_back();
    
    //Refresh progress...
    progressBar(j, n_peaks, "=", " ");
    
    //Peak mass range search
    PeakColumnRange current_colRange = calculateSinglePeakRange( pmasses[j], pwidths[j], last_mass_index);
    
    // Emplace next peak column task immediately into the freed slot
    active_tasks.push_back({
      std::async(
        std::launch::async, 
        &RapidBinning::calculatePeakMatColumn, 
        this, 
        assigned_thread_id, 
        current_colRange.start_col, 
        current_colRange.count, 
        j 
        ),
      assigned_thread_id
    });
  }
  
  // Drain remaining active threads
  for (auto& task : active_tasks) {
    task.future.get();
  }
  active_tasks.clear();
  
  Rcout << std::endl << "Binning complete!" << std::endl;
}


// startingIonIndex: starting ion/mass channel index in the binary stream axis.
// ionCount: Total number of contiguous mass channels to decode and integrate for this peak
// target_peak_index: column index in pMat where decoded pixel intensities for this mass peak will be written.
void RapidBinning::calculatePeakMatColumn(unsigned int thread_id, 
                                          unsigned int startingIonIndex, unsigned int ionCount, 
                                          unsigned int target_peak_index)
{
  if (ionCount == 0 )
  {
    return; 
  }  
  
  // 1. Decode stream slice into native C++ vector (0 R memory allocations)
  std::vector<double> decodedImg = vec_XBin[thread_id]->decodeImgStream2Buffer(startingIonIndex, ionCount);
  
  // 2. Direct raw memory pointer write to pre-allocated R matrix column
  size_t total_mat_rows = pMat.nrow();
  double* col_write_ptr = REAL(pMat) 
    + (static_cast<size_t>(target_peak_index) * total_mat_rows) 
    + img_row_offset;
  
  // 3. Iterate over pre-converted coordinates and extract pixel intensities
  size_t n_pixels = vec_XYCoords.size();
  unsigned int width = vec_XBin[thread_id]->getImgWidth();
  
  for (size_t i = 0; i < n_pixels; ++i)
  {
    const auto& point = vec_XYCoords[i];
    
    // Extract intensity from decoded raster matrix using 0-based (X, Y) pixel lookup
    // and store directly into the designated row offset in pMat
    col_write_ptr[i] = decodedImg[point.x + (point.y * width)];
  }
}


RapidBinning::PeakColumnRange RapidBinning::calculateSinglePeakRange( double target_mass, double target_width, size_t& last_mass_index)
{
  double half_w = target_width / 2.0;
  double min_m = target_mass - half_w;
  double max_m = target_mass + half_w;
  size_t n_axis = massAxis.length();
  
  // Advance cursor to min_m
  while (last_mass_index < n_axis && massAxis[last_mass_index] < min_m)
  {
    ++last_mass_index;
  }
  
  // EDGE CASE: Entire target peak is beyond the upper bound of massAxis
  if (last_mass_index >= n_axis) {
    return PeakColumnRange{ 0, 0 }; // count = 0 signals peak not found
  }
  
  size_t start_idx = last_mass_index;
  
  // 2. Adjust lower bound across gaps (nearest point before gap)
  if (start_idx > 0)
  {
    double dist_curr = std::abs(massAxis[start_idx] - min_m);
    double dist_prev = std::abs(min_m - massAxis[start_idx - 1]);
    
    if (dist_prev < dist_curr)
    {
      start_idx = start_idx - 1;
    }
  }

  // 3. Scan forward to upper bound (max_m)
  size_t end_idx = start_idx;
  while (end_idx < n_axis && massAxis[end_idx] <= max_m)
  {
    ++end_idx;
  }
  unsigned int count = static_cast<unsigned int>(end_idx - start_idx);
  
  // If start_idx lands past max_m (e.g. inside a massive gap where no points exist in [min_m, max_m])
  if (massAxis[start_idx] > max_m) {
    count = 0;
  }
  
  // 4. Update last_mass_index state cursor for the NEXT peak search
  last_mass_index = start_idx;
  
  // 5. Construct and return range struct
  return PeakColumnRange{
    static_cast<unsigned int>(start_idx),
    count
  };
}


//' Cload_RapidBinning
 //' 
 //' Extracts and bins peak images across an MSI dataset in parallel into a pre-allocated R matrix.
 //' 
 //' @param rMSIobj An rMSI object prefilled with parsed imzML/binary metadata.
 //' @param pmasses A numeric vector of target peak mass centers (sorted).
 //' @param pwidths A numeric vector of peak window widths (same length as pmasses).
 //' @param pMat A pre-allocated NumericMatrix where extracted peak columns will be written.
 //' @param img_row_offset Integer specifying the starting row offset in pMat (1-based from R).
 //' @param number_of_threads Number of worker threads for parallel decoding.
 //' 
 // [[Rcpp::export]]
 void C_RapidBinning(
     List rMSIobj,
     NumericVector pmasses,
     NumericVector pwidths,
     NumericMatrix pMat,
     unsigned int img_row_offset,
     int number_of_threads)
 {
   // 1. Input Validation
   if (pmasses.length() != pwidths.length()) {
     Rcpp::stop("ERROR in C_RapidBinning: 'pmasses' and 'pwidths' must have the same length.");
   }
   
   if (img_row_offset < 1) {
     Rcpp::stop("ERROR in C_RapidBinning: 'img_row_offset' must be >= 1 (1-based R indexing).");
   }
   
   if (number_of_threads < 1) {
     Rcpp::stop("ERROR in C_RapidBinning: 'number_of_threads' must be at least 1.");
   }
   
   // 2. Convert R 1-based row offset to C++ 0-based offset
   unsigned int c_img_row_offset = img_row_offset - 1;
   
   // 3. Thread-safe execution wrapper
   try
   {
     // Instantiate the RapidBinning pipeline object
     RapidBinning binning_engine(
         rMSIobj,
         pmasses,
         pwidths,
         pMat,
         static_cast<unsigned int>(number_of_threads),
         c_img_row_offset );
     
     // Execute parallel extraction and direct matrix population
     binning_engine.run();
   }
   catch (const std::exception &e)
   {
     // Safely propagate C++ exceptions to R session without crashing the R process
     Rcpp::stop(e.what());
   }
   catch (...)
   {
     Rcpp::stop("ERROR in C_RapidBinning: An unknown exception occurred during processing.");
   }
 }
