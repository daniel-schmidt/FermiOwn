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
							Nf( numFlavours ),
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
							acceptanceCounter(0),
							fWeight( "weight" + std::to_string(lambda) + ".dat" )
{
	generateAllowedConfs( numFlavours );
	intConfIndex = std::uniform_int_distribution<int>( 0, allowedConfs.rows()-1 );

	slac.erase( kxiab );
	detOld = slac.det();
	slac.setFull();
	//	std::string filename = "weight";
	//	filename = filename + std::to_string(lambda) + ".dat";
}

FermiBoolMetropolis::~FermiBoolMetropolis() {
}

void  FermiBoolMetropolis::generateAllowedConfs( size_t numFlavours ) {
	size_t numCols = (numFlavours*numFlavours)*2;
	size_t numConfigs = 1;
	numConfigs = numConfigs << numCols;
//	std::cout << "Num cols " << numCols << std::endl;
	allowedConfs = Eigen::MatrixXi::Zero( numConfigs, numCols );
//	std::cout << numConfigs << std::endl;
//	std::cout << allowedConfs << std::endl;
	size_t numAllowedConfs = 0;
	for( size_t conf = 0; conf < numConfigs; conf++ ) {
		size_t bits = conf;
		size_t bits2 = conf;
		FieldBoolean fConf( 1, 2, numFlavours, NULL, zeroInit );

		for( size_t spin = 0; spin < 2; spin++ ) {
			for( size_t a = 0; a < numFlavours; a++ ) {
				for( size_t b = 0; b < numFlavours; b++ ) {
					if( bits2 % 2 == 1 ) {
						fConf.setValue( true, 0, spin, a, b );
					}
					bits2 = bits2 >> 1;
				}
			}
		}

		if( !fConf.constraintViolated(0) ) {
			for( size_t bit = 0; bit < numCols; bit++ ) {
				if( bits % 2 == 1 ) {
					allowedConfs( numAllowedConfs, numCols-1-bit ) = 1;
				}
				bits = bits >> 1;
			}
			numAllowedConfs ++;
		}
	}
	allowedConfs.conservativeResize( numAllowedConfs, Eigen::NoChange );
//	std::cout << "All configs:" << std::endl << allowedConfs << std::endl;
}

void FermiBoolMetropolis::sumAllConfs() {
	size_t confsPerX = allowedConfs.rows();
	size_t numConfigs = std::pow( confsPerX, lat.getVol() );
//	std::cout << "numConfigs: " << numConfigs << std::endl;

	Complex sum = 0.;

	for( size_t conf = 0; conf < numConfigs; conf++ ) {
		size_t confNum = conf;
//		std::cout << "Conf " << conf << " confNum: ";

		Eigen::VectorXi confList = Eigen::VectorXi::Zero( lat.getVol() );
		size_t count = 0;
		while( confNum >= confsPerX ) {
			confList( count ) = confNum % confsPerX;
			confNum /= confsPerX;
			count ++;
		}

		confList( count ) = confNum;
//		std::cout << std::endl<< confList <<std::endl;

		for( size_t x = 0; x < lat.getVol(); x++ ) {
			size_t newConfIndex = confList( x );
			for( size_t spin = 0; spin < 2; spin++ ) {
				for( size_t a = 0; a < Nf; a++ ) {
					for( size_t b = 0; b < Nf; b++ ) {
						kxiab.setValue( bool(allowedConfs( newConfIndex, spin*Nf*Nf + b*Nf + a)), x, spin, a, b );
					}
				}
			}
		}
//		kxiab.Print();

		sum += calculateWeight();
	}
	std::cout << 2./kappa << "\t" << sum << std::endl;
}

Complex FermiBoolMetropolis::calculateWeight() {
	size_t k = kxiab.sumAll();
	size_t ntilde = kxiab.countOffdiagonal2();
	double factor = 1.;
	for( size_t x = 0; x< lat.getVol(); x++ ) {
		auto nx = kxiab.countSummedSpin( x );
		factor *= getHypergeometricFactor( nx(1), nx(2) );
	}
	factor *= std::pow( 2., ntilde );
	Complex dw = std::pow( Complex(-kappa), double(k)/2. );

	return factor*dw;
}

Complex FermiBoolMetropolis::calculateWeightChange() {
	size_t k = kxiab.sumAll();
	size_t ntilde = kxiab.countOffdiagonal2();
	double factor = 1.;
	for( size_t x = 0; x< lat.getVol(); x++ ) {
		auto nx = kxiab.countSummedSpin( x );
		factor *= getHypergeometricFactor( nx(1), nx(2) );
	}
	factor *= std::pow( 2., ntilde );
	Complex dw = std::pow( Complex(-kappa), double(k)/2. );

	FieldBoolean changed = kxiab.different(oldField);

//	slac.erase( kxiab );
	slac.update( kxiab, changed, SlacOperatorMatrix::separateUpdate );
	det = slac.det();
//	slac.setFull();

	return factor*dw*det;
}

bool FermiBoolMetropolis::updateField() {

//	std::cout << "Field before update " << std::endl;
//	kxiab.Print();
	oldField = kxiab;
	SlacOperatorMatrix oldSlac = slac;
//	int dk = -kxiab.sumAll();
//	int dntilde = -kxiab.countOffdiagonal2();

	Complex weight = 1./calculateWeight();

//	slac.erase( kxiab );
	det = slac.det();
//	std::cout << "before: " << det << " ";
//	slac.setFull();
	weight /= det;

	// draw random point, spin and 2 flavours


	for( int cnt=0; cnt < 4; cnt++ ) {
//	nxOld = kxiab.countSummedSpin( x );
		int x = intV_dist(*rndGen);
		int newConfIndex = intConfIndex( *rndGen );
	//	std::cout << "Setting x=" << x << " to conf " << newConfIndex << std::endl;
		for( size_t spin = 0; spin < 2; spin++ ) {
			for( size_t a = 0; a < Nf; a++ ) {
				for( size_t b = 0; b < Nf; b++ ) {
					kxiab.setValue( bool(allowedConfs( newConfIndex, spin*Nf*Nf + b*Nf + a)), x, spin, a, b );
				}
			}
		}
	}

	weight *= calculateWeight();
//	slac.erase( kxiab );

	FieldBoolean changed = kxiab.different(oldField);
	slac.update( kxiab, changed, SlacOperatorMatrix::separateUpdate );
	det = slac.det();
//	Eigen::MatrixXcd upSlac = slac.getMatrix();
//	slac.setFull();
//	slac.erase( kxiab );
//	det = slac.det();
//	if( !upSlac.isApprox( slac.getMatrix() ) ) {
//		std::cout << "Matrix falsch!" << std::endl << upSlac << std::endl<< std::endl << slac.getMatrix() << std::endl;
//	}
//	std::cout << abs( updateDet-det ) << std::endl;
//	if( abs( updateDet-det ) > 1e-10 ) {
//		std::cout << "update: " << updateDet << " ";
//		std::cout << "full: " << det << std::endl;
//	}
//	slac.setFull();
	weight *= det;
	bool accepted = accept( weight );

	if( !accepted ) slac = oldSlac;
//	std::cout << "Accepted: " << accepted << ", field after update:" << std::endl << std::endl;
//	kxiab.Print();
	return accepted;
}

void FermiBoolMetropolis::updateNaive( size_t x ) {
	if( intSpin_dist(*rndGen) ) {

		int spin = intSpin_dist(*rndGen);
		int a = intNf_dist(*rndGen);
		int b = intNf_dist(*rndGen);

		kxiab.invert( x, spin, a, b );
	}
}

bool FermiBoolMetropolis::updateField( size_t x, size_t spin, size_t a, size_t b ) {

	kxiab.invert( x, spin, a, b );
	kxiab.enforceConstraint( x, spin, a, b );

//	std::cout << std::endl << "First change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
	if( a == b ) {
//		 updating with same flavour,
//		 choose, if we update another spin (un-/setting a kxaa=2) or another flavour (setting two kxaa=1, keeping n1 even)
		size_t c = intNf_dist(*rndGen);
		size_t spin2 = 0;
		if( a==c ) {
			// same flavour drawn again, we have to update the second spin (otherwise we would revert our initial update)
			spin2 = 1-spin;
//			std::cout << "New c drawn." << std::endl;
		} else {
			// different flavour, spin can be chosen randomly
			spin2 = intSpin_dist(*rndGen);
//			std::cout << "New spin and b drawn." << std::endl;
		}
		kxiab.invert( x, spin2, c, c );
		kxiab.enforceConstraint( x, spin2, c, c );
	} else {
		// updating with two different flavours, enforce kxab = kxba
		spin = intSpin_dist(*rndGen);
//		std::cout << "New spin drawn." << std::endl;
		kxiab.invert( x, spin, b, a );
		kxiab.enforceConstraint( x, spin, b, a );
	}
//	std::cout << "Second change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
//	kxiab.Print();
	bool success = true;
	if( kxiab.constraintViolated( x ) ) {
//		std::cout << "x violated! Reset and continue loop..." << std::endl;
		kxiab = oldField;
		success = false;
	}

	return success;
}

Complex FermiBoolMetropolis::calculateWeight( int dk, int dntilde ) {
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
//	std::cout << "accepted: " << accepted << std::endl;
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
