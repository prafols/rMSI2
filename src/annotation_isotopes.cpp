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
using namespace std;



int factorialIterativeCustom(int n)
{
  long f = 1;
  for (int i=1; i<=n; ++i)
    f *= i;
  return f;
}


/*************************************************************************
 * 
 *                            IsotopeAnnotation
 *     
 **************************************************************************/


IsotopeAnnotation::IsotopeAnnotation(MSIparameters* myMSIparams, std::vector<PeakInformation>* myPeakInformationVector)
{
  MSIparams = myMSIparams;
  peakInformationVector = myPeakInformationVector;
}


void IsotopeAnnotation::run()
{
  for(int i = 0; i < peakInformationVector->size(); i++) // Iterate for each peakInformation.
  {
    PeakInformation* myPeakOfInterest = &peakInformationVector->operator[](i);

     if( !myPeakOfInterest->inPattern() ) // Search if is not assigned in an isotopic pattern.
     {
       if( findCandidates(myPeakOfInterest, 1) )  // Check if it has possible isotope candidates.
       {
         if( scoreCandidates(myPeakOfInterest, 1) ) // Check if a candidate pass the ILS test.
         {
          isotopicPatternsVector.emplace_back(myPeakOfInterest, bestCandidate); // Instantiate an IsotopicPattern
          while( findCandidates(myPeakOfInterest, isotopicPatternsVector.back().getNumberOfPeaksInPattern()) ) // Search for more isotopes for the isotopic pattern.
          {
           if( scoreCandidates(myPeakOfInterest, isotopicPatternsVector.back().getNumberOfPeaksInPattern()) ) // Check if a candidate pass the ILS test.
           {
             isotopicPatternsVector.back().addIsotope(bestCandidate); // Add isotope to the isotopic pattern
           }
           else break;
          }
        }
      }
    }
  }
  Rcpp::Rcout << "\nIsotopeAnnotation run successful! " << isotopicPatternsVector.size() << " isotopic patterns detected." << "\n";
}


bool IsotopeAnnotation::findCandidates(PeakInformation* peakOfInterest, int isotopeNumber)
{
  bool candidatesFlag = false;
  isotopeCandidatesVector.clear();

  if( ((peakOfInterest->getColumnIndex())+1) <= peakInformationVector->size() )
  {
    for( int i = ((peakOfInterest->getColumnIndex())+1); i < peakInformationVector->size(); i++ )
    {
      PeakInformation* myCandidatePeak = &peakInformationVector->operator[](i);
      double massDifference = myCandidatePeak->getMass() - (peakOfInterest->getMass() + 1.003355*isotopeNumber);
      double massWindow = MSIparams->tolerance*peakOfInterest->getMass()/1000000;
      
      if( (massDifference >= -massWindow) && (massDifference<= massWindow) && !myCandidatePeak->inPattern())
      {
        isotopeCandidatesVector.emplace_back(myCandidatePeak);
        candidatesFlag = true;
      }
      
      if( (massDifference > massWindow) )
      {
        return candidatesFlag;
      }
    }
  }
  return candidatesFlag;
}


bool IsotopeAnnotation::scoreCandidates(PeakInformation* peakOfInterest, int isotopeNumber)
{
  bool testPass = false; // flag to indicate if a candidate pass the ILS test 
  double tmpMaximumILS = 0; // slot to save the previous maximum ILS between the candidates
  double scoreMass = 0, scoreMorphology = 0, isotopeRatio = 0, ILS = 0, ppm = 0, scoreIntensity = 0;
  int carbonAtoms = 0;
  // std::vector<int> nullPixelsPOI; //vector containing the pixels with a 0 in intensity for the peak of interest.
  // for(int i = 0; i < MSIparams->numberPixels; i++)
  // {
  //   if( peakOfInterest->getPixelIntensity(i) == 0 )
  //   {
  //     nullPixelsPOI.emplace_back(i);
  //   }
  // }
  for(int i = 0; i < isotopeCandidatesVector.size(); i++) //Starting the evaluation of all candidates.
  {
    PeakInformation* candidatePeak = isotopeCandidatesVector[i]; // pointer to the candidate to test.

    // std::vector<int> nullPixelsCandidate; //vector containing the pixels with a 0 in intensity for a candidate peak.
    // std::vector<int> valuePixels; //vector containing the pixels with intensity values at the same tame between the POI and a candidate.
    // for(int j = 0; j <= MSIparams->numberPixels; j++)
    // {
    //   if( candidatePeak->getPixelIntensity(i) == 0 )
    //   {
    //     nullPixelsCandidate.emplace_back(i);
    //   }
    // }
    //
    // for(int j = 0; j <= MSIparams->numberPixels; j++)
    // {
    //   if( !(find(nullPixelsPOI.begin(), nullPixelsPOI.end(), j) || find(nullPixelsCandidate.begin(), nullPixelsCandidate.end(), j)) )
    //   {
    //     valuePixels.emplace_back(j);
    //   }
    // }

    //*********************************** Mass score
    scoreMass = 0;
    ppm = fabs(1000000*(candidatePeak->getMass()-(peakOfInterest->getMass()+(1.003355*isotopeNumber)))/(peakOfInterest->getMass()+(1.003355*isotopeNumber)));

    scoreMass  =  ppm/(sqrt((1/log(2)))*MSIparams->tolerance); //Adjusting the maximum tolerance to score 0.5
    scoreMass *= -scoreMass;
    scoreMass  = exp(scoreMass);

    //*********************************** Morphology score
    scoreMorphology = peakOfInterest->peakCorrelation(candidatePeak, MSIparams->numberPixels);
    scoreMorphology = (scoreMorphology < 0) ?  0 : scoreMorphology;
    
    //*********************************** Intensity score
    double scoreIntHMDB = 0, scoreIntPeptideAtlas = 0, A = 0, B = 0;
    
    //Linear model interpolation
    for(int j = 0; j < MSIparams->numberPixels; j++)
    {
      A += (peakOfInterest->getPixelIntensity(j) - peakOfInterest->getAverageIntensity())*(candidatePeak->getPixelIntensity(j) - candidatePeak->getAverageIntensity());
      B += (peakOfInterest->getPixelIntensity(j) - peakOfInterest->getAverageIntensity())*(peakOfInterest->getPixelIntensity(j) - peakOfInterest->getAverageIntensity());
    }
    isotopeRatio = (A/B < 0) ?  0.00005 : A/B;  //Isotope ratio from the data

    //intercept = candidatePeak->getAverageIntensity() - ratio_slope*peakOfInterest->getAverageIntensity();

    double ratioHMDB  = 0.000702 * peakOfInterest->getMass() - 0.03851; //HMDB
    double ratioPeptideAtlas  = 0.0004712976 * peakOfInterest->getMass() - 0.0080417319; //PeptideAtlas

    double adjustedIsotopeRatio = (isotopeNumber == 1) ? isotopeRatio : (pow(factorialIterativeCustom(isotopeNumber)*isotopeRatio, 1 / double(isotopeNumber)) + 0.01081573 * ((isotopeNumber - 1) / 2));
    carbonAtoms = (int)std::nearbyint(adjustedIsotopeRatio * (0.9893 / 0.0107));

    if( peakOfInterest->getMass() <= 1300 ) //Condition for evaluation with the HMDB model
    {
      scoreIntHMDB  = fabs(ratioHMDB - adjustedIsotopeRatio)/(sqrt(1/(-log(0.7)))*0.2);  //Adjusting 0.2 difference to score 0.7
      scoreIntHMDB *= -scoreIntHMDB;
      scoreIntHMDB  = exp(scoreIntHMDB);
      scoreIntensity = scoreIntHMDB;
    }

    if( peakOfInterest->getMass() >= 600 ) //Condition for evaluation with the peptideAtlas model
    {
      scoreIntPeptideAtlas  = fabs(ratioPeptideAtlas - adjustedIsotopeRatio)/(sqrt(1/(-log(0.7)))*0.2258983);  //Adjusting 4*S  to score 0.7, this includes 99% of the peptides in the library
      scoreIntPeptideAtlas *= -scoreIntPeptideAtlas;
      scoreIntPeptideAtlas  = exp(scoreIntPeptideAtlas);
    }

    if(scoreIntPeptideAtlas >= scoreIntHMDB)
    {
      scoreIntensity = scoreIntPeptideAtlas;
    }
    
    //ILS
    ILS = scoreIntensity * scoreMorphology * scoreMass;
    
    if ( ILS > tmpMaximumILS )
    {
      tmpMaximumILS = ILS;
    }

    if( tmpMaximumILS >= MSIparams->scoreThreshold )
    {
      testPass = true;
      bestCandidate = candidatePeak;
    }
  }
  
  if(testPass)
  {
    bestCandidate->writeIsotopeScoreTable(ILS, scoreMorphology, scoreIntensity, scoreMass, ppm, isotopeRatio, carbonAtoms);
  }

  return testPass;
}


void IsotopeAnnotation::reportIsotopicPatterns()
{
  for(int i = 0; i<isotopicPatternsVector.size() ;i++)
  {
    Rcpp::Rcout << "\n--- IsotopicPattern instance "<< i <<" ---\n";
    isotopicPatternsVector[i].reportToConsole();
  }
}


/*************************************************************************
 * 
 *                            IsotopicPattern
 *     
 **************************************************************************/


IsotopicPattern::IsotopicPattern(PeakInformation* monoisotope, PeakInformation* isotope)
{
  numberOfPeaks = 2;
  monoisotope->setIsotopicPosition(0);
  isotope->setIsotopicPosition(1);
  isotopicPatternPeaks.emplace_back(monoisotope);
  isotopicPatternPeaks.emplace_back(isotope);
  adductNetworkID = -1;
  adductElement = "";
}


void IsotopicPattern::addIsotope(PeakInformation* isotope)
{
  isotope->setIsotopicPosition(numberOfPeaks);
  numberOfPeaks = numberOfPeaks + 1;
  isotopicPatternPeaks.emplace_back(isotope);
}


int IsotopicPattern::getNumberOfPeaksInPattern()
{
  return(numberOfPeaks);
}


void IsotopicPattern::reportToConsole()
{
  for(int i = 0; i < numberOfPeaks; i++)
  {
    if(i == 0)
    {
      Rcpp::Rcout << "M+"<< isotopicPatternPeaks[i]->getIsotopicPosition() <<
        ": Peak number " << isotopicPatternPeaks[i]->getColumnIndex() <<
          ", m/z " << isotopicPatternPeaks[i]->getMass() << " \n";
    }
    else
    {
     isotopeScoreTableStr* scores = isotopicPatternPeaks[i]->getIsotopeScoreTable();
      Rcpp::Rcout << "M+"<< isotopicPatternPeaks[i]->getIsotopicPosition() <<
        ": Peak number " << isotopicPatternPeaks[i]->getColumnIndex() <<
          ", m/z " << isotopicPatternPeaks[i]->getMass() <<
            ", ILS "<< scores->ILS <<
              ", carbonAtoms "<< scores->carbonAtoms <<
                ", massError "<< scores->massError << " \n";
    }
  }

}


bool IsotopicPattern::inNetworkAsAdduct(adductElement* element)
{
  if(adductNetworkID == -1)
  {
    return false;
  }
  
  if( adduct->mass == element->mass )
  {
    return true;
  }
    
  return false;
}

int getAdductNetworkId()
{
  return adductNetworkID;
}
