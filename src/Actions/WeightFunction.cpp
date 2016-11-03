/*
 * WeightFunction.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include "WeightFunction.h"

namespace FermiOwn {

WeightFunction::WeightFunction( const ThirringKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling ) :
				kxiab( boolField ),
				initialField( kxiab ),
				savedState( kxiab ),
				dslash( timeSize, spatialSize, dim, numFlavours ),
				V( timeSize*pow(spatialSize, dim) ),
				kappa( 1./coupling ),
				needsKeepDecision( false )
{}

WeightFunction::~WeightFunction() {}

Complex WeightFunction::calculateWeight() {
	if( needsKeepDecision ) {
		std::cerr << "Error in calculateWeight(): WeightFunction is not updated properly! Call keep/reset before calculateWeight()." << std::endl;
		exit(1);
	}

	size_t k = kxiab.sumAll();
	size_t ntilde = kxiab.countOffdiagonal2();
	double factor = 1.;
	for( size_t x = 0; x< V; x++ ) {
		auto nx = kxiab.countSummedSpin( x );
		factor *= getHypergeometricFactor( nx(1), nx(2) );
	}
	factor *= std::pow( 2., ntilde );
	Complex dw = std::pow( Complex(-kappa), double(k)/2. );

	dslash.calculateUpdateMatrices( kxiab, kxiab.different( initialField ) );
	Complex det = dslash.getDet();
	//	dslash.reset();

	needsKeepDecision = true;

	return factor*dw*det;
}

Complex WeightFunction::updateWeight( const std::set< size_t > & changedAt ) {

	if( needsKeepDecision ) {
		std::cerr << "Error in updateWeight(): WeightFunction is not updated properly! Call keep/reset before calculateWeight()." << std::endl;
		exit(1);
	}

	int dk = kxiab.sumAll() - savedState.sumAll(); //TODO: should be able to obtain from kxiab.different();
	int dntilde = kxiab.countOffdiagonal2() - savedState.countOffdiagonal2();
	double factor = 1.;

	for( size_t x : changedAt ) {
//		std::cout << x << std::endl;
		Eigen::ArrayXi newNx = kxiab.countSummedSpin( x );
		Eigen::ArrayXi oldNx = savedState.countSummedSpin( x );
		factor *= getHypergeometricFactor( newNx(1), newNx(2) ) / getHypergeometricFactor( oldNx(1), oldNx(2) );
	}


	factor *= std::pow( 2., dntilde );
	Complex dw = std::pow( Complex(-kappa), double(dk)/2. );

	dslash.calculateUpdateMatrices( kxiab, kxiab.different( savedState ) );
	Complex det = dslash.updateDet();

//	kxiab.different( savedState ).Print();
//	std::cout << "dk=" << dk << " dntilde=" << dntilde << " det=" << det << " dw=" << dw << " factor=" << factor << std::endl;

	needsKeepDecision = true;
	return factor*dw*det;
}

//				std::cout << "Matrix before erase: " << std::endl << slac.getMatrix() << std::endl << std::endl;
//	slac.erase( kxiab );
//	//				std::cout << slac.getMatrix() << std::endl << std::endl;
//	det = slac.det();
//	slac.setFull();
//
//	double factor = getHypergeometricFactor( nxNew(1), nxNew(2) )  * getHypergeometricFactor( nyNew(1), nyNew(2) );
//	factor /= getHypergeometricFactor( nxOld(1), nxOld(2) ) * getHypergeometricFactor( nyOld(1), nyOld(2) );
//
//	//					std::cout << "\tfactor1=" << factor;
//
//	//	if( dntilde != 0 ) std::cerr << "\tdntilde = " << dntilde << std::endl;
//	factor *= std::pow( 2., dntilde );
//	Complex dw = std::pow( Complex(-kappa), double(dk)/2. );
//	Complex detRatio = det/detOld;
//	fWeight << std::real(detRatio) << "\t" << std::imag(detRatio) << "\t" << dw << "\t" << std::real(factor*dw*(detRatio)) << "\t" << std::imag(factor*dw*(detRatio)) << std::endl;
////	std::cout << "reDet=" << std::real(detRatio) << "\t imDet=" << std::imag(detRatio) << "\t dk=" << dk << "\t kappa=" << kappa << "\t dw=" << dw << "\t factor=" << factor << "\t dnt=" << dntilde << "\t" << std::fabs(factor*dw*(detRatio)) << std::endl;
//	return factor*dw*(det/detOld);

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
