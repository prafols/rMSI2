/*************************************************************************
 *     rMSI2 - R package for MSI data handling and visualization
 *     Copyright (C) 2022 Lluc Sementé Fernández
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

#include "annotation_rMSIannotation.h"
using namespace Rcpp;

//Constructor
rMSIannotation::rMSIannotation(MSIparameters* myMSIparams, NumericMatrix myPeakIntensityMatrix, 
                               NumericVector myPeakMassAxis, NumericVector myAverageSpectrumMassAxis, 
                               NumericVector myAverageSpectrumIntensity, DataFrame myAdductElementsDF)
{
  MSIparams = myMSIparams;
  //Rcpp::Rcout << "  MSIparams ... \n";
  
  peakIntensityMatrix = new double*[MSIparams->numberPixels];
  for(int i = 0; i < MSIparams->numberPixels; i++)
  {
    peakIntensityMatrix[i] = new double[MSIparams->numberPeaks];
    for(int j = 0; j < MSIparams->numberPeaks; j++)
    {
      peakIntensityMatrix[i][j] = myPeakIntensityMatrix(i,j);
    }
  }
  //Rcpp::Rcout << "  peakIntensityMatrix ... \n";
  
  peakMassAxis = new double[MSIparams->numberPeaks]; 
  for(int i = 0; i < MSIparams->numberPeaks; i++)
  {
    peakMassAxis[i] = myPeakMassAxis(i);
  }  
  //Rcpp::Rcout << "  peakMassAxis ... \n";
  
  averageSpectrumMassAxis = new double[MSIparams->numberSpectrumDataPoints]; 
  for(int i = 0; i < MSIparams->numberSpectrumDataPoints; i++)
  {
    averageSpectrumMassAxis[i] = myAverageSpectrumMassAxis(i);
  }  
  //Rcpp::Rcout << "  averageSpectrumMassAxis ... \n";
  
  averageSpectrumIntensity = new double[MSIparams->numberSpectrumDataPoints]; 
  for(int i = 0; i < MSIparams->numberSpectrumDataPoints; i++)
  {
    averageSpectrumIntensity[i] = myAverageSpectrumIntensity(i);
  } 
  //Rcpp::Rcout << "  averageSpectrumIntensity ... \n";
  
  adductElementsDF = myAdductElementsDF;
  Rcpp::Rcout << "0.DONE\n";
}

//Deconstructor
rMSIannotation::~rMSIannotation()
{
  //Deallocate memory from Heap
  Rcpp::Rcout << "4.Deallocate rMSIannotation ... \n";
  delete[] averageSpectrumIntensity;
  delete[] averageSpectrumMassAxis;
  delete[] peakMassAxis;
  for(int i = 0; i < MSIparams->numberPixels; i++)
  {
    delete[] peakIntensityMatrix[i];
  }
  delete[] peakIntensityMatrix;
  Rcpp::Rcout << "4.DONE" << "\n";
  Rcpp::Rcout << "*** Leaving C++ ***" << "\n";
}

//Generate PeakInformation instances for each peak in the peak matrix
void rMSIannotation::generatePeakInformation()
{
  int tmp = 0;
  peakInformationVector.reserve(MSIparams->numberPeaks);
  for(int i = 0; i < MSIparams->numberPeaks; i++)
  {
    peakInformationVector.emplace_back(peakIntensityMatrix, peakMassAxis, MSIparams->numberPixels, tmp);
    tmp++;
  }
}

// Starting function of rMSIannotation
void rMSIannotation::run()
{
  Rcpp::Rcout << "1.Generate PeakInformation instances ... \n";
  generatePeakInformation();
  Rcpp::Rcout << "1.DONE" << "\n";
  
  
  Rcpp::Rcout << "2.Generate IsotopeAnnotation instance ... \n";
  IsotopeAnnotation myIsotopeAnnotation(MSIparams, &peakInformationVector);
  Rcpp::Rcout << "2.DONE" << "\n";
  
  Rcpp::Rcout << "3.Starting isotope annotation ... \n";
  myIsotopeAnnotation.run();
  myIsotopeAnnotation.reportIsotopicPatterns();
  Rcpp::Rcout << "3.DONE" << "\n";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
void rMSIannotation_C(int numberPixels, 
                      int numberPeaks,
                      int numberSpectrumDataPoints, 
                      int tolerance,
                      double scoreThreshold,
                      bool toleranceInScans,
                      DataFrame adductElementsDF,
                      NumericMatrix peakIntensityMatrix, 
                      NumericVector peakMassAxis,
                      NumericVector averageSpectrumMassAxis, 
                      NumericVector averageSpectrumIntensity)
{
  Rcpp::Rcout << "*** Entering C++ ***" << "\n";
  MSIparameters myMSIparams;
  myMSIparams.numberPixels = numberPixels;
  myMSIparams.numberPeaks = numberPeaks;
  myMSIparams.numberSpectrumDataPoints = numberSpectrumDataPoints;
  myMSIparams.tolerance = tolerance;
  myMSIparams.scoreThreshold = scoreThreshold;
  myMSIparams.toleranceInScans = toleranceInScans;
  
  Rcpp::Rcout << "0.Instantiate rMSIannotation" << "\n";
  rMSIannotation myAnnotation(&myMSIparams, peakIntensityMatrix, peakMassAxis,
                              averageSpectrumMassAxis, averageSpectrumIntensity, adductElementsDF); 
  
  Rcpp::Rcout << "Run rMSIannotation" << "\n";
  myAnnotation.run();
}







