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
PeakInformation::PeakInformation(double** peakIntensityMatrix, double* peakMassAxis, int numberPixels, int myColumnIndex)
{
  isotopicPosition = -1;
  mass = peakMassAxis[myColumnIndex];
  columnIndex = myColumnIndex;
  intensityVector = new double[numberPixels];
  averageIntensity = 0;
  for(int pixel = 0; pixel < numberPixels; pixel++)
  {
    intensityVector[pixel] = peakIntensityMatrix[pixel][columnIndex];
    averageIntensity += intensityVector[pixel];
  }
  averageIntensity = averageIntensity/(numberPixels);

  // if(columnIndex < 10)
  // {
  //   Rcpp::Rcout << "  Starting instance " << columnIndex << "\n";
  //   Rcpp::Rcout << "  ...mass "<< mass <<"\n";
  //   Rcpp::Rcout << "  ...columnIndex " << columnIndex << "\n";
  //   Rcpp::Rcout << "  ...averageIntensity "<< averageIntensity <<" \n";
  //   Rcpp::Rcout << "  ...instance " << columnIndex << " DONE \n";
  // }
  
}

//Destructor
PeakInformation::~PeakInformation()
{
  //Rcpp::Rcout << "Deallocate PeakInformation" << "\n";
  delete[] intensityVector;
}

//Check if a peak is the M+0 peak of a pattern
bool PeakInformation::isMonoisotope()
{
  if(isotopicPosition == 0)
  {
    return true;
  }
  return false;
}

//Check if a peak is part of an isotopic pattern
bool PeakInformation::inPattern()
{
  if(isotopicPosition != -1)
  {
    return true;
  }
  return false;
}

//Intensity correlation between two peaks
double PeakInformation::peakCorrelation(PeakInformation* targetPeak, int numberPixels)
{
  double sumA = 0, sumB = 0, sumAA = 0, sumAB = 0,  sumBB = 0;
  double numerator = 0, denominator = 0, correlation = 0;
    
  for(int i = 0; i < numberPixels ; i++)
  {
      sumAB  += intensityVector[i]*targetPeak->intensityVector[i];
      sumA   += intensityVector[i];
      sumB   += targetPeak->intensityVector[i];
      sumAA  += intensityVector[i]*intensityVector[i];
      sumBB  += targetPeak->intensityVector[i]*targetPeak->intensityVector[i];
  }
  numerator = (numberPixels * sumAB - (sumA * sumB));
  denominator = (numberPixels * (sumAA)-(sumA * sumA)) * (numberPixels * (sumBB)-(sumB * sumB));
  denominator = sqrt(denominator);
  correlation = numerator/denominator;
  return(correlation);
}

//Set the isotopic position of a peak in the pattern
void PeakInformation::setIsotopicPosition(int position)
{
  isotopicPosition = position;
}

double PeakInformation::getMass()
{
  return(mass);
}

double PeakInformation::getAverageIntensity()
{
  return(averageIntensity);
}

double PeakInformation::getPixelIntensity(int pixel)
{
  return(intensityVector[pixel]);
}

int PeakInformation::getColumnIndex()
{
  return(columnIndex);
}

int PeakInformation::getIsotopicPosition()
{
  return(isotopicPosition);
}

void PeakInformation::writeIsotopeScoreTable(double ILS, double morphologyScore, double intensityScore, double massScore, double massError, double intensityRatio, int carbonAtoms)
{
  isotopeScoreTable.ILS = ILS;
  isotopeScoreTable.morphologyScore = morphologyScore;
  isotopeScoreTable.intensityScore = intensityScore;
  isotopeScoreTable.massScore = massScore;
  isotopeScoreTable.massError = massError;
  isotopeScoreTable.intensityRatio = intensityRatio;
  isotopeScoreTable.carbonAtoms = carbonAtoms;
}

isotopeScoreTableStr* PeakInformation::getIsotopeScoreTable()
{
  return(&isotopeScoreTable);
}