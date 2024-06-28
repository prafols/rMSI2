/*************************************************************************
 *     rMSI2 - R package for MSI data handling and visualization
 *     Copyright (C) 2023 Lluc Sementé Fernández
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




/*************************************************************************
 * 
 *                            AdductAnnotation
 *     
 **************************************************************************/



AdductAnnotation::AdductAnnotation(MSIparameters* myMSIparams, std::vector<PeakInformation>* myPeakInformationVector, 
                                   std::vector<IsotopicPattern>* myIsotopicPatternsVector, DataFrame myAdductElementsDF)
{
  MSIparams = myMSIparams;
  peakInformationVector = myPeakInformationVector;
  isotopicPatternsVector= myIsotopicPatternsVector;
  adductElementVector.clear();
  fillAdductElementVector(myAdductElementsDF);
}



void AdductAnnotation::fillAdductElementVector(DataFrame adductElementDataframe)
{
  NumericVector massVector = adductElementDataframe[1];
  CharacterVector nameVector = adductElementDataframe[0];
  for(int i = 0; i<adductElementDataframe.nrows() ;i++)
  {
    adductElement myElement;
    myElement.mass = massVector[i];
    myElement.element = nameVector[i];
    adductElementVector.emplace_back(myElement);
  }
}

void AdductAnnotation::run()
{
  for(int i = 0; i < isotopicPatternsVector->size(); i++) // Iterate for each IsotopicPattern.
  {
    IsotopicPattern* myIsotopicPatternOfInterest = &isotopicPatternsVector->operator[](i);
    for( int j = 0; j < adductElementVector.size(); j++)
    {
      if( !myIsotopicPatternOfInterest->inNetworkAsAdduct(&adductElementVector[j]) ) // Search if is not assigned in an AdductNetwork.
      {
        if (findCandidatesAsAdduct( myIsotopicPatternOfInterest,  &adductElementVector[j]))
        {
          if (scoreCandidatesAsAdduct( myIsotopicPatternOfInterest, &adductElementVector[j] ))
          {
            if(bestAdductCandidate->inNetworkAsAdduct(bestAdductElement)) // If the best candidate is already in an AdductNetwork 
            {
              adductNetworkVector[bestAdductCandidate->getAdductNetworkId()].addAdduct(myIsotopicPatternOfInterest, &adductElementVector[j]);
            }
            else
            {
              adductNetworkVector.emplace_back( myIsotopicPatternOfInterest, &adductElementVector[j],
                                                bestAdductCandidate, bestAdductElement ); //New AdductNetwork instance
            }
          }
        }
      }
    }
  }
  Rcpp::Rcout << "\nAdductAnnotation run successful! " << adductNetworkVector.size() << " adduct networks created." << "\n";
}




