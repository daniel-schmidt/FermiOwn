/*
 * FermiBoolMetropolis.cpp
 *
 *  Created on: 31.03.2016
 *      Author: dschmidt
 */

#include "FermiBoolMetropolis.h"

namespace FermiOwn {

FermiBoolMetropolis::FermiBoolMetropolis( FieldBoolean& boolField, const Lattice & lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator ) :
							kappa( 2./lambda ),
							Nf( numFlavours ),
							kxiab( boolField ),
							lat(lattice),
							rndGen( randomGenerator ),
							confGen( 2, numFlavours, randomGenerator ),
							weightFun( boolField, lat.getTimeSize(), lat.getSpaceSize(), lat.getDim(), numFlavours, lambda ),
							uni_real_dist(),
							intV_dist( 0, lat.getVol()-1 ),
							int2mu_dist( 0, 2*lat.getDim()-1 ),
							intNf_dist( 0, numFlavours-1 ),
							intSpin_dist( 0, 1 ),
							oldField(kxiab),
							acceptanceCounter(0),
							fWeight( "weight" + std::to_string(lambda) + ".dat" )
{
	confGen.generateAllowedConfs();
//	generateAllowedConfs( numFlavours );
//	intConfIndex = std::uniform_int_distribution<int>( 0, allowedConfs.rows()-1 );
//
//	slac.erase( kxiab );
//	detOld = slac.det();
//	slac.setFull();
	//	std::string filename = "weight";
	//	filename = filename + std::to_string(lambda) + ".dat";
}

FermiBoolMetropolis::~FermiBoolMetropolis() {
}

bool FermiBoolMetropolis::updateField() {

	oldField = kxiab;

//	Complex weight = 1./weightFun.calculateWeight();

	weightFun.saveState();
	// draw random points, which get a new configuration

	std::set<size_t> changedPoints;

	for( int cnt=0; cnt < 4; cnt++ ) {
		int x = intV_dist(*rndGen);
		RowVectorXb newConf = confGen.getRandomConf();
		kxiab.setRow( newConf, x );
		changedPoints.insert( size_t( x ) );
	}

	Complex weight = weightFun.updateWeight( changedPoints );

//	std::cout << weight << std::endl;
	bool accepted = accept( weight );

	return accepted;
}

bool FermiBoolMetropolis::accept( Complex weight ) {
	double r = uni_real_dist(*rndGen);
	bool accepted = false;
	if( std::fabs(weight) > r ) {
		accepted = true;
		acceptanceCounter++;
		detOld = det;
		weightFun.keep();
	} else {
		kxiab = oldField;
		weightFun.reset();
	}
//	std::cout << "accepted: " << accepted << std::endl;
	return accepted;
}


//void FermiBoolMetropolis::updateNaive( size_t x ) {
//	if( intSpin_dist(*rndGen) ) {
//
//		int spin = intSpin_dist(*rndGen);
//		int a = intNf_dist(*rndGen);
//		int b = intNf_dist(*rndGen);
//
//		kxiab.invert( x, spin, a, b );
//	}
//}

//bool FermiBoolMetropolis::updateField( size_t x, size_t spin, size_t a, size_t b ) {
//
//	kxiab.invert( x, spin, a, b );
//	kxiab.enforceConstraint( x, spin, a, b );
//
////	std::cout << std::endl << "First change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
//	if( a == b ) {
////		 updating with same flavour,
////		 choose, if we update another spin (un-/setting a kxaa=2) or another flavour (setting two kxaa=1, keeping n1 even)
//		size_t c = intNf_dist(*rndGen);
//		size_t spin2 = 0;
//		if( a==c ) {
//			// same flavour drawn again, we have to update the second spin (otherwise we would revert our initial update)
//			spin2 = 1-spin;
////			std::cout << "New c drawn." << std::endl;
//		} else {
//			// different flavour, spin can be chosen randomly
//			spin2 = intSpin_dist(*rndGen);
////			std::cout << "New spin and b drawn." << std::endl;
//		}
//		kxiab.invert( x, spin2, c, c );
//		kxiab.enforceConstraint( x, spin2, c, c );
//	} else {
//		// updating with two different flavours, enforce kxab = kxba
//		spin = intSpin_dist(*rndGen);
////		std::cout << "New spin drawn." << std::endl;
//		kxiab.invert( x, spin, b, a );
//		kxiab.enforceConstraint( x, spin, b, a );
//	}
////	std::cout << "Second change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
////	kxiab.Print();
//	bool success = true;
//	if( kxiab.constraintViolated( x ) ) {
////		std::cout << "x violated! Reset and continue loop..." << std::endl;
//		kxiab = oldField;
//		success = false;
//	}
//
//	return success;
//}

//Complex FermiBoolMetropolis::calculateWeight( int dk, int dntilde ) {
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
//}


} /* namespace FermiOwn */
