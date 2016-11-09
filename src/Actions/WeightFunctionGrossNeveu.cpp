/*
 * WeightFunctionGrossNeveu.cpp
 *
 *  Created on: 09.11.2016
 *      Author: dschmidt
 */

#include "WeightFunctionGrossNeveu.h"

namespace FermiOwn {

WeightFunctionGrossNeveu::WeightFunctionGrossNeveu( const GrossNeveuKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling ) :
		BasicWeightFunctionTemplate( boolField, timeSize, spatialSize, dim, numFlavours, coupling ),
		change( dslash, kfield ),
		locWeights( numFlavours + 1 )
{
	// generate list of local weights
	for( size_t kx = 0; kx <= 2*numFlavours; kx += 2 ) {
		double val = (1.+kx)/2.;
		locWeights[kx/2] = std::tgamma( val ) * std::pow( kappa, val );
	}
//	std::cout << "Weights:" << std::endl;
//	for( double w : locWeights ) std::cout << w << std::endl;
}

WeightFunctionGrossNeveu::~WeightFunctionGrossNeveu() {
	// TODO Auto-generated destructor stub
}

Complex WeightFunctionGrossNeveu::calculateWeight() {
	if( needsKeepDecision ) {
		std::cerr << "Error in WeightFunctionGrossNeveu::calculateWeight(): WeightFunction is not updated properly! Call keep/reset before calculateWeight()." << std::endl;
		exit(1);
	}
	double factor = 1.;
	idxMap KxNxMap = kfield.tally();
	for( idxPair kxnx : KxNxMap ) {
		if( kxnx.first%2 != 0 ) {
			std::cerr << "In WeightFunctionGrossNeveu::calculateWeight(): Constraint kx even violated! kx = " << kxnx.first << std::endl;
		}
		factor *= std::pow( locWeights[kxnx.first/2], kxnx.second );
	}

	auto chlist = change.calculateDifference( initialField );
//	initialField.Print();
//	change.Print();
//	kfield.Print();
	dslash.calculateUpdateMatrices( chlist );
	Complex det = dslash.getDet();

	needsKeepDecision = true;
//	std::cout << "Factor: " << factor << " det: " << det << std::endl;
	return det*factor;
}

Complex WeightFunctionGrossNeveu::updateWeight( const std::set<size_t>& changedAt ) {
	if( needsKeepDecision ) {
		std::cerr << "Error in WeightFunctionGrossNeveu::updateWeight(): WeightFunction is not updated properly! Call keep/reset before calculateWeight()." << std::endl;
		exit(1);
	}
}

} /* namespace FermiOwn */
