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


#ifndef RAPID_BINNING_H
#define RAPID_BINNING_H

#include <vector>
#include <Rcpp.h>
#include "rMSIXBin.h"

class RapidBinning
{
  public:
  
    // rMSIobject: rMSI image object
    // peak_masses: the detected peak masses
    // peak_widths: the width of eacj peak detected (same length as the peak_masses)
    // PeakMatrix: An R matrix pre-allocated that will be modified during the peak-binning
    // nThreads: number of processing threads
    // PeakMatrixRowOffset: peak matrix offset if multiple images are merged
    RapidBinning(Rcpp::List rMSIobject, Rcpp::NumericVector peak_masses, Rcpp::NumericVector peak_widths, Rcpp::NumericMatrix PeakMatrix, unsigned int nThreads, unsigned int PeakMatrixRowOffset);

    ~RapidBinning();
    
    void run();
  
  private:
    Rcpp::List img; //Points to the original rMSIobject
    Rcpp::NumericVector pmasses; //Points to the original peak_masses
    Rcpp::NumericVector pwidths; //Points to the original peak_widths 
    Rcpp::NumericMatrix pMat; //Points to the original PeakMatrix 
    Rcpp::NumericVector massAxis; //The MSI data mass axis
    unsigned int n_peaks; //length of pmasses and pwidths
    unsigned int n_threads;
    unsigned int img_row_offset; //starting row index offset in pMat. Accounts for multi-dataset concatenation where pixel data from image $i$ is appended below preceding image datasets.
    std::vector<std::unique_ptr<rMSIXBin>> vec_XBin; //Each thread needs and independent rMSIXBin object.
  
    //The method to run in parallel
    void calculatePeakMatColumn(unsigned int thread_id, 
                                unsigned int startingIonIndex, unsigned int ionCount, 
                                unsigned int target_peak_index);
    
    //Struct to store each worker with its thread
    struct TaskHandle {
      std::future<void> future;
      unsigned int thread_id;
    };
    
    //Struct to hold X, Y coords pairs
    struct XYPoint {
      unsigned int x;
      unsigned int y;
    };
    std::vector<XYPoint> vec_XYCoords;
    
    //Struct for peak mass searches
    struct PeakColumnRange {
      unsigned int start_col;
      unsigned int count;
    };
    PeakColumnRange calculateSinglePeakRange( double target_mass, double target_width, size_t& last_mass_index);
    
};

#endif
