/*
 * FermiOwnGraph.cpp
 *
 *  Created on: 15.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include "Lattice.h"
#include "FieldBoolean.h"
#include "SlacOperatorMatrix.h"

namespace FermiOwn {

double getHypergeometricFactor( int flavour ) {
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

double getHypergeometricFactor( int n1, int n2 ) {
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

}

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 8 ) {
		std::cerr << "Need 7 arguments: Nt, Ns, numThermal, numMeasures, upPerMeasure, coupling, Nf" << std::endl;
		exit(1);
	}

	// read parameters from console
	size_t dim = 3;
	size_t Nt = atoi( argv[1] );
	size_t Ns = atoi( argv[2] );
	size_t V = Nt*Ns*Ns;
	int numThermal = atoi( argv[3] );
	int numMeasures = atoi( argv[4] );
	int upPerMeasure = atoi( argv[5] );
	double lambda = atof( argv[6] );
	int Nf = atof( argv[7] );
	int numSpin = 2;
	// initialize classes
	Lattice lat( Nt, Ns, dim );
	SlacOperatorMatrix slac( Nt, Ns, dim, Nf );
//	std::cout << slac.getMatrix() << std::endl << std::endl;

	// initialize random number generator and distributions
	std::ranlux48 gen;
	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist( 0, V-1 );
	std::uniform_int_distribution<int> int2mu_dist( 0, 2*dim-1 );
	std::uniform_int_distribution<int> intNf_dist( 0, Nf-1 );
	std::uniform_int_distribution<int> intSpin_dist( 0, 1 );

	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
		double kappa = 2./lambda;

		// initialize single-flavour part to true, the rest to false
		FieldBoolean kxab( lat.getVol(), numSpin, Nf, &gen, zeroInit );

		for( int a = 0; a < Nf; a++ ) {
			for( int spin = 0; spin < numSpin; spin++ ) {
				for( size_t x = 0; x < lat.getVol(); x++ ) {
					kxab.setValue( true, x, spin, a, a );
				}
			}
		}

		double av_k = 0.;
		double accrate = 0;
		slac.erase( kxab );
//		std::cout << slac.getMatrix() << std::endl << std::endl;
		Complex detOld = slac.det();
		slac.setFull();

		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				// draw random point, spin and 2 flavours
				int x = intV_dist(gen);
				int spin = intSpin_dist(gen);
				int a = intNf_dist(gen);
				int b = intNf_dist(gen);

//				 calculate old value of na
				Eigen::ArrayXi nxOld = kxab.countSummedSpin( x );
				Eigen::ArrayXi nxNew;
				Eigen::ArrayXi nyOld = Eigen::ArrayXi::Zero( numSpin+1 );
				Eigen::ArrayXi nyNew = Eigen::ArrayXi::Zero( numSpin+1 );

				FieldBoolean kxabOld(kxab);

//				std::cout << "Field before update:" << std::endl;
//				kxab.Print();

				// update at the drawn index
				int dk = -kxab.sumAll();
				kxab.invert( x, spin, a, b );
				int dk2 = dk + kxab.sumAll();
//				std::cout << "Field after first update:" << std::endl;
//				kxab.Print();

				if( a == b ) {
					// updating with same flavour,
					// choose, if we update another spin (setting a kxaa=2) or another point (setting two kxaa=1, keeping n1 even)
//					double r = uni_real_dist(gen);
//					if( 0.5 > r ) {
						kxab.invert( x, 1-spin, a, a );
//					} else {
						int mu = int2mu_dist(gen);
						int y = lat.getNeighbours( x )[mu];
						nyOld = kxab.countSummedSpin( y );
						kxab.invert( y, spin, a, a );
						kxab.invert( y, 1-spin, a, a );

						if( kxab.constraintViolated( y ) ) {
							kxab = kxabOld;
							std::cout << "violated! Reset and continue loop..." << std::endl;
							continue;
						}
						nyNew = kxab.countSummedSpin( y );
//					}
				} else {
					// updating with two different flavours, enforce kxab = kxba
					int spin2 = intSpin_dist(gen);
					if( dk2 < 0 ) {
						// we have to delete another spin in kxba, choose one randomly and check,
						// if it is set, otherwise the other one must be set, since we have only two spins.
						if( kxab.getValue( x, spin2, b, a ) ) {
							kxab.invert( x, spin2, b, a );
						} else {
							kxab.invert( x, 1-spin2, b, a );
						}
					} else {
						// we have to set another spin in kxba as before
						if( !kxab.getValue( x, spin2, b, a ) ) {
							kxab.invert( x, spin2, b, a );
						} else {
							kxab.invert( x, 1-spin2, b, a );
						}
					}
				}
				dk += kxab.sumAll();
				nxNew = kxab.countSummedSpin( x );
//				std::cout << "Field after second update:" << std::endl;
//				kxab.Print();
				if( kxab.constraintViolated( x ) ) {
					kxab = kxabOld;
					std::cout << "violated! Reset and continue loop..." << std::endl;
					continue;
				}
//				std::cout << "Matrix before erase: " << std::endl << slac.getMatrix() << std::endl << std::endl;
				slac.erase( kxab );
//				std::cout << slac.getMatrix() << std::endl << std::endl;
				Complex det = slac.det();
				slac.setFull();
//				std::cout << "nxOld=" << std::endl << nxOld << std::endl << "nxNew=" << std::endl << nxNew << std::endl
//						<< "nyOld=" << std::endl << nyOld << std::endl << "nyNew=" << std::endl << nyNew << std::endl;
//				std::cout << "dk: " << dk << "\tdetOld: " << detOld << "\tdet: " << det;
				double factor = 1.;
				factor /= getHypergeometricFactor( nxOld(1), nxOld(2) ) * getHypergeometricFactor( nyOld(1), nyOld(2) );
//				std::cout << "\tfactor1=" << factor;
				factor *= getHypergeometricFactor( nxNew(1), nxNew(2) ) * getHypergeometricFactor( nyNew(1), nyNew(2) );
				double dw = std::pow(kappa, double(dk)/2.);
				Complex weight =  factor*dw*(det/detOld);

				double r = uni_real_dist(gen);
//				bool accepted = false;
				if( std::fabs(weight) > r ) {
//					accepted = true;
					accrate++;
					detOld = det;
				} else {
					kxab = kxabOld;
				}
//				std::cout << "\tfactor: " << factor << "\tdw: " << dw << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;
			}
			if( measure >= numThermal ) av_k += kxab.sumAll()/double(2*V);
		}
		av_k/=double(numMeasures);
		accrate /= double((numMeasures+numThermal)*upPerMeasure);
		std::cerr << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
