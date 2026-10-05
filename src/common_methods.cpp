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

#include "common_methods.h"

std::string parse_xml_uuid(std::string uuid)
{
  std::size_t ipos = uuid.find('{');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  do
  {
    ipos = uuid.find('-'); 
    if( ipos != std::string::npos )
    {
      uuid.erase(ipos, 1);
    } 
  } while ( ipos != std::string::npos);
  ipos = uuid.find('}');
  if( ipos != std::string::npos )
  {
    uuid.erase(ipos, 1);
  }
  for( unsigned int i=0; i < uuid.length(); i++)
  {
    uuid[i] = toupper(uuid[i]);
  }
  return uuid;
}


//' Accumulate Overall Mean and Base Peak Spectra Across Mass Spectrometry Images
 //' 
 //' Interpolates an individual image's mean and base peak spectra onto a unified 
 //' target mass axis using single-precision 1D linear interpolation and updates 
 //' the running overall maximum vectors in-place.
 //' 
 //' @param common_mass A numeric vector containing the target/reference mass axis (m/z values).
 //' @param img_mass A numeric vector containing the mass axis (m/z values) of the current image.
 //' @param img_mean A numeric vector containing the mean intensity spectrum of the current image.
 //' @param img_base A numeric vector containing the base peak spectrum of the current image.
 //' @param overall_mean A numeric vector storing the accumulated maximum mean intensities (modified in-place).
 //' @param overall_base A numeric vector storing the accumulated maximum base peak intensities (modified in-place).
 //' 
 //' @return Void. The `overall_mean` and `overall_base` vectors are updated directly in memory by reference.
 //' 
 //' @keywords internal
 // [[Rcpp::export]]
void CaccumulateOverallSpectra(Rcpp::NumericVector common_mass,
                                 Rcpp::NumericVector img_mass,
                                 Rcpp::NumericVector img_mean,
                                 Rcpp::NumericVector img_base,
                                 Rcpp::NumericVector overall_mean,
                                 Rcpp::NumericVector overall_base)
{
  
  int n_common = common_mass.size();
  int n_src = img_mass.size();
  
  if (n_src == 0 || n_common == 0) return;
  
  // Temporary scratch vectors for 1D interpolated results
  std::vector<double> interp_mean(n_common);
  std::vector<double> interp_base(n_common);
  
  // 1. Interpolate mean spectrum
  mlinterp::interp(
    &n_src, n_common,
    img_mean.begin(), interp_mean.data(),
    img_mass.begin(), common_mass.begin()
  );
  
  // 2. Interpolate base spectrum
  mlinterp::interp(
    &n_src, n_common,
    img_base.begin(), interp_base.data(),
    img_mass.begin(), common_mass.begin()
  );
  
  // 3. In-place parallel max update (NO R-level memory allocation)
  for (int i = 0; i < n_common; ++i)
  {
    if (interp_mean[i] > overall_mean[i]) overall_mean[i] = interp_mean[i];
    if (interp_base[i] > overall_base[i]) overall_base[i] = interp_base[i];
  }
}
