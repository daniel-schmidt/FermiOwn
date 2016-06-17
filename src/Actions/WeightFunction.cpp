/*
 * WeightFunction.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include "WeightFunction.h"

namespace FermiOwn {

WeightFunction::WeightFunction( FieldBoolean& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling ) :
		kxiab( boolField ),
//		slac( timeSize, spatialSize, dim, numFlavours ),
		V( timeSize*pow(spatialSize, dim) ),
		kappa( 1./coupling )
{}

WeightFunction::~WeightFunction() {}

Complex WeightFunction::calculateWeight() {
	size_t k = kxiab.sumAll();
	size_t ntilde = kxiab.countOffdiagonal2();
	double factor = 1.;
	for( size_t x = 0; x< V; x++ ) {
		auto nx = kxiab.countSummedSpin( x );
		factor *= getHypergeometricFactor( nx(1), nx(2) );
	}
	factor *= std::pow( 2., ntilde );
	Complex dw = std::pow( Complex(-kappa), double(k)/2. );

//	slac.erase( kxiab );

//	Complex det = slac.det();
	Complex det = 0.; //TODO: change this
	return factor*dw*det;
}

Complex WeightFunction::updateWeight() {
//	size_t k = kxiab.sumAll();
//	size_t ntilde = kxiab.countOffdiagonal2();
//	double factor = 1.;
//	for( size_t x = 0; x< V; x++ ) {
//		auto nx = kxiab.countSummedSpin( x );
//		factor *= getHypergeometricFactor( nx(1), nx(2) );
//	}
//	factor *= std::pow( 2., ntilde );
//	Complex dw = std::pow( Complex(-kappa), double(k)/2. );
//
//	FieldBoolean changed = kxiab.different(oldField);
//
////	slac.erase( kxiab );
//	slac.update( kxiab, changed, SlacOperatorMatrix::separateUpdate );
//	Complex det = slac.det();
////	slac.setFull();
//
//	return factor*dw*det;
}

/*=====================================================================
 * Private Functions
 *=====================================================================*/


double WeightFunction::getHypergeometricFactor( int flavour ) {
	switch( flavour ){
	case 0: return 1.;
	case 1: return 1.5;
	case 2: return 2.75;
	case 3: return 53./8.;
	case 4: return 345./16.;
	case 5: return 2947./32.;
	case 6: return 31411./64.;
	default: std::cout << "Value of the confluent hypergeometric function for a=" << flavour << " not implemented!" << std::endl;
	return -1.;
	}
}

double WeightFunction::getHypergeometricFactor( int n1, int n2 ) {
	if( n1%2==1 ) return 0.;
	switch( n1 ) {
	case 0: return getHypergeometricFactor( n2 );
	case 2:
		switch( n2 ) {
		case 0:	return 0.5;
		case 1:	return 5./4.;
		case 2:	return 31./8.;
		case 3: return 239./16.;
		case 4: return 2257./32.;
		default: std::cout << "Value of the confluent hypergeometric function for n2=" << n2 << " not implemented!" << std::endl;
		return -1.;
		}
		break;
		case 4:
			switch( n2 ) {
			case 0: return 3./4.;
			case 1: return 21./8.;
			case 2: return 177./16.;
			case 3: return 1779./32.;
			case 4: return 21003./64.;
			default: std::cout << "Value of the confluent hypergeometric function for n2=" << n2 << " not implemented!" << std::endl;
			return -1.;
			}
			break;
			default: std::cout << "Value of the confluent hypergeometric function for n1=" << n1 << " and " << n2 << " not implemented!" << std::endl;
			return -1.;
	}
}

} /* namespace FermiOwn */
