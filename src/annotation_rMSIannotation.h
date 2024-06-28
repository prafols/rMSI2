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

#ifndef RMSIANNOTATION_H
#define RMSIANNOTATION_H

#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

typedef struct
{
  int numberPixels; //Number of pixels in the peak matrix
  int numberPeaks; //Number of peaks in the peak matrix
  int numberSpectrumDataPoints; //Number of data points in the average spectrum
  int tolerance;             //Mass tolerance for peak candidates
  double scoreThreshold;     //Final score needed to pass the isotope test
  bool toleranceInScans;     //If true, tolerance is applied in scans, else is in ppm
}MSIparameters;

typedef struct
{
  double ILS;                //Isotopic likelihood scores for each isotope.
  double morphologyScore; 
  double intensityScore; 
  double massScore; 
  double massError;
  double intensityRatio;
  int carbonAtoms;
}isotopeScoreTableStr;

typedef struct
{
  std::string element;
  double mass;
}adductElement;

//*PeakInformation
//*
//* Class to store the information of mass peak of the peak matrix
//*
//*
//**//

class PeakInformation
{
public:
  //Methods
  PeakInformation(double** peakIntensityMatrix,
                  double* peakMassAxis,
                  int numberPixels,
                  int myColumnIndex); //Constructor
  ~PeakInformation(); //Destructor
  bool isMonoisotope(); //Check if a peak is the M+0 peak of a pattern
  bool inPattern(); //Check if a peak is part of an isotopic pattern
  double peakCorrelation(PeakInformation* targetPeak, int numberPixels); //Intensity correlation between two peaks
  void setIsotopicPosition(int position); //Set the isotopic position of a peak in the pattern
  double getMass(); 
  double getPixelIntensity(int pixel);
  double getAverageIntensity();
  int getColumnIndex();
  int getIsotopicPosition();
  void writeIsotopeScoreTable(double ILS, double morphologyScore, double intensityScore, double massScore, double massError, double intensityRatio, int carbonAtoms);
  isotopeScoreTableStr* getIsotopeScoreTable();
  
private:
  //Variables
  isotopeScoreTableStr isotopeScoreTable;
  int isotopicPosition; //Numeric position of the peak in an isotopic pattern. 0 -> M+0, 1 -> M+1, etc... -1 -> Not part of a pattern. -2 ->
  int columnIndex; //Index of the peak in the peak matrix
  double mass; //Mass of the peak
  double averageIntensity; //Average intensity of the peak
  double* intensityVector; //Intensity values of the peak at each pixel
  adductElement* element;
};

//IsotopicPattern
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class IsotopicPattern
{
public:
  //variables
  
  //methods
  IsotopicPattern(PeakInformation* monoisotope, PeakInformation* isotope);
  int getNumberOfPeaksInPattern();
  void addIsotope(PeakInformation* isotope);
  void reportToConsole();
  bool inNetworkAsAdduct(adductElement* element);
  int getAdductNetworkId();
  
private:
  //variables
  int numberOfPeaks; //Total number of peaks in the pattern.
  int adductNetworkID;
  adductElement* adduct;
  std::vector<PeakInformation*> isotopicPatternPeaks; //Pointers to PeakInformation instances belonging to the pattern 
  
  //methods  
};

//IsotopeAnnotation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class IsotopeAnnotation
{
public:
  IsotopeAnnotation(MSIparameters* MSIparams,
                    std::vector<PeakInformation>* myPeakInformationVector); //Constructor.
  void run();   //Program execution function.
  void reportIsotopicPatterns();

  
private:
  MSIparameters* MSIparams;
  PeakInformation* bestCandidate; //pointer to the peakInformation corresponding to the best candidate.
  std::vector<PeakInformation>* peakInformationVector; //vector with pointers to all PeakInformation instances.
  std::vector<PeakInformation*> isotopeCandidatesVector; //vector with pointers to the PeakInformation instances to check for an isotopic pattern.
  std::vector<IsotopicPattern> isotopicPatternsVector; 
  
  //methods
  bool findCandidates(PeakInformation* peakOfInterest, int isotopeNumber); // Search candidates at mass distance for isotopes.
  bool scoreCandidates(PeakInformation* peakOfInterest, int isotopeNumber); //Scores all the candidates and retrieves the best one in bestCandidate.
};

//AdductNetwork
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AdductNetwork
{
public:
  //variables
  
  //methods
  AdductNetwork(IsotopicPattern* adductA, adductElement* elementA, IsotopicPattern* adductB, adductElement* elementB);
  void addAdduct(IsotopicPattern* adduct, adductElement* element);
  
private:
  //variables

  //methods  
};

//AdductAnnotation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AdductAnnotation
{
public:
  AdductAnnotation(MSIparameters* myMSIparams,
                   std::vector<PeakInformation>* myPeakInformationVector, 
                   std::vector<IsotopicPattern>* myIsotopicPatternsVector,
                   DataFrame myAdductElementsDF); //Constructor.
  void run();   //Program execution function.
  void reportAdductNetworks();
  
private:
  MSIparameters* MSIparams;
  IsotopicPattern* bestAdductCandidate;
  adductElement* bestAdductElement;
  std::vector<PeakInformation>* peakInformationVector; //vector with pointers to all PeakInformation instances.
  std::vector<IsotopicPattern>* isotopicPatternsVector; 
  std::vector<IsotopicPattern*> adductCandidatesVector; 
  std::vector<adductElement> adductElementVector; // TODO: ensure they are ordered by mass.
  std::vector<AdductNetwork> adductNetworkVector;
  
  //methods
  void fillAdductElementVector(DataFrame adductElementDataframe);
  bool findCandidatesAsAdduct(IsotopicPattern* isotopicPatternOfInterest, adductElement* candidateElement);
  bool findCandidates(IsotopicPattern* isotopicPatternOfInterest); // Search candidates at mass distance for isotopes.
  bool scoreCandidatesAsAdduct(IsotopicPattern* isotopicPatternOfInterest, adductElement* candidateElement); //Scores all the candidates and retrieves the best one in bestCandidate.
};


//rMSIannotation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class rMSIannotation
{
public:

  //methods
  rMSIannotation(MSIparameters* myMSIparams,
                 NumericMatrix myPeakIntensityMatrix,
                 NumericVector myPeakMassAxis,
                 NumericVector myAverageSpectrumMassAxis,
                 NumericVector myAverageSpectrumIntensity,
                 DataFrame myAdductElementsDF);
  ~rMSIannotation();
  void run();

private:
  //variables
  MSIparameters* MSIparams;
  double** peakIntensityMatrix; //Peak intensity data of the peak matrix
  double* averageSpectrumMassAxis; //Mass axis of the average spectrum
  double* peakMassAxis; //Mass axis of the peak matrix
  double* averageSpectrumIntensity; //Intensity vector of the average spectrum
  std::vector<PeakInformation> peakInformationVector; //vector with all PeakInformation instances
  DataFrame adductElementsDF;
  
  //methods
  void generatePeakInformation();
};


#endif