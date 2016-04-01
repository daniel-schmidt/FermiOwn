/*
 * FermiBoolMetropolis.cpp
 *
 *  Created on: 31.03.2016
 *      Author: dschmidt
 */

#include "FermiBoolMetropolis.h"

namespace FermiOwn {

FermiBoolMetropolis::FermiBoolMetropolis( FieldBoolean& boolField, SlacOperatorMatrix& slacMat, const Lattice & lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator ) :
			kappa( 2./lambda ),
			kxiab( boolField ),
			oldField(kxiab),
			slac(slacMat),
			lat(lattice),
			rndGen( randomGenerator ),
			uni_real_dist(),
			intV_dist( 0, lat.getVol()-1 ),
			int2mu_dist( 0, 2*lat.getDim()-1 ),
			intNf_dist( 0, numFlavours-1 ),
			intSpin_dist( 0, 1 ),
			acceptanceCounter(0)
{
	slac.erase( kxiab );
	detOld = slac.det();
	slac.setFull();
}

FermiBoolMetropolis::~FermiBoolMetropolis() {
}

bool FermiBoolMetropolis::updateField() {

	oldField = kxiab;

	// draw random point, spin and 2 flavours
	int x = intV_dist(*rndGen);
	int spin = intSpin_dist(*rndGen);
	int a = intNf_dist(*rndGen);
	int b = intNf_dist(*rndGen);

	int dk = -kxiab.sumAll();
	int dntilde = -kxiab.countOffdiagonal2();
	nxOld = kxiab.countSummedSpin( x );
	updateField( x, spin, a, b );
	nxNew = kxiab.countSummedSpin( x );

	int mu = int2mu_dist(*rndGen);
	int y = lat.getNeighbours( x )[mu];
	nyOld = kxiab.countSummedSpin( y );
	updateField( y, spin, a, b );
	nyNew = kxiab.countSummedSpin( y );

	dk += kxiab.sumAll();
	dntilde += kxiab.countOffdiagonal2();

	Complex weight = calculateWeight( dk, dntilde );

	if( imag(weight) > ZERO_TOL || ( std::fabs(weight) > ZERO_TOL && std::real(weight) < 0 ) ) {
		std::cerr << "Warning, non-positive weight: " << weight << std::endl;
	}

	bool accepted = accept( weight );
	return accepted;
}

bool FermiBoolMetropolis::updateField( size_t x, size_t spin, size_t a, size_t b ) {
	kxiab.invert( x, spin, a, b );
	if( a == b ) {
		// updating with same flavour,
		// choose, if we update another spin (setting a kxaa=2) or another point (setting two kxaa=1, keeping n1 even)
		//					double r = uni_real_dist(gen);
		//					if( 0.5 > r ) {
		kxiab.invert( x, 1-spin, a, a );
		//					} else {
		//					}
	} else {
		// updating with two different flavours, enforce kxab = kxba

//		int spin2 = intSpin_dist(*rndGen);

		//		if( dk2 < 0 ) {
		//			// we have to delete another spin in kxba to keep the sum constraint, choose one randomly and check,
		//			// if it is set, otherwise the other one must be set, since we have only two spins.
		//			if( kxiab.getValue( x, spin2, b, a ) ) {
		//				kxiab.invert( x, spin2, b, a );
		//			} else {
		//				kxiab.invert( x, 1-spin2, b, a );
		//			}
		//		} else {
		//			// we have to set another spin in kxba as before
		//			if( !kxiab.getValue( x, spin2, b, a ) ) {
		//				kxiab.invert( x, spin2, b, a );
		//			} else {
		//				kxiab.invert( x, 1-spin2, b, a );
		//			}
		//		}
	}
	bool success = true;
	if( kxiab.constraintViolated( x ) ) {
		std::cout << "x violated! Reset and continue loop..." << std::endl;
		kxiab = oldField;
		success = false;
	}
	return success;
}

Complex FermiBoolMetropolis::calculateWeight( int dk, int dntilde ) {
	//				std::cout << "Matrix before erase: " << std::endl << slac.getMatrix() << std::endl << std::endl;
	slac.erase( kxiab );
	//				std::cout << slac.getMatrix() << std::endl << std::endl;
	det = slac.det();
	slac.setFull();

	double factor = getHypergeometricFactor( nxNew(1), nxNew(2) ) * getHypergeometricFactor( nyNew(1), nyNew(2) );
	factor /= getHypergeometricFactor( nxOld(1), nxOld(2) ) * getHypergeometricFactor( nyOld(1), nyOld(2) );
//					std::cout << "\tfactor1=" << factor;

//	std::cout << "\tdntilde = " << dntilde << std::endl;
	factor *= pow( 2, -dntilde );
	double dw = std::pow(kappa, double(dk)/2.);
	return factor*dw*(det/detOld);
}

bool FermiBoolMetropolis::accept( Complex weight ) {
	double r = uni_real_dist(*rndGen);
	bool accepted = false;
	if( std::fabs(weight) > r ) {
		accepted = true;
		acceptanceCounter++;
		detOld = det;
	} else {
		kxiab = oldField;
	}
	return accepted;
}

double FermiBoolMetropolis::getHypergeometricFactor( int flavour ) {
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

double FermiBoolMetropolis::getHypergeometricFactor( int n1, int n2 ) {
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
