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
#include "SlacOperatorMatrix.h"
#include "Lattice.h"

namespace FermiOwn {

Complex calcDet( SlacOperatorMatrix& dslac, const VectorXb& kx ) {
	for( int x = 0; x < kx.size(); x++ ) {
		if( kx(x) ) {
			dslac.deletePoint(x);
		}
	}
	return pow(I,dslac.getMatrix().rows())*dslac.det();
}
int linIndex(int N, int i, int j) {
	if( i == j) {
		return i;
	} else {
		if( i > j ) {
			int tmp = j;
			j = i;
			i = tmp;
		}
		return N+N*(N-1)/2 - (N-i)*(N-i-1)/2+j-i-1;
	}
}

Complex calcDet( std::vector<SlacOperatorMatrix> dslac, const MatrixXb& kxab0, const MatrixXb& kxab1, int Nf ) {
	Complex det = 1.;
	for( size_t b = 0; b < Nf; b++ ) {
		for( size_t a = 0; a < b; a++ ) {
			VectorXb kx0 = kxab0.col( linIndex( Nf, a, b ) );
			VectorXb kx1 = kxab1.col( linIndex( Nf, a, b ) );
			dslac[a].eraseCols( kx0, kx1 );
			dslac[b].eraseRows( kx0, kx1 );
		}
	}
	for( size_t a = 0; a < Nf; a++ ) {
		VectorXb kx0 = kxab0.col( a );
		VectorXb kx1 = kxab1.col( a );
		dslac[a].eraseCols( kx0, kx1 );
		dslac[a].eraseRows( kx0, kx1 );
		det *= pow(I,dslac[a].getMatrix().rows())*dslac[a].det();
	}
	return det;
}


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

bool constraintViolated( const MatrixXb& kxab0, const MatrixXb& kxab1, int Nf, int x ) {
	bool violated = false;
	for( int m = 0; m < Nf; m++ ) {
		int constrSum = 0;
		for( int n = 0; n < Nf; n++ )
		{
			constrSum += kxab0( x, linIndex( Nf, m, n ) ) + kxab1( x, linIndex( Nf, m, n ) );
		}
		//		std::cout << "constraint for a=" << m << " is " << constrSum << std::endl;
		if( constrSum > 1 ) violated = true;
	}
	return violated;
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

	// initialize classes
	Lattice lat( Nt, Ns, dim );
	std::vector<SlacOperatorMatrix> slacNf( Nf, SlacOperatorMatrix( Nt, Ns, dim ) );

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

		// Data layout: first Nf columns hold kxaa for a=1,..,Nf, then kxab with a<b orderd (1,2), (1,3), (2,3), ...
		// initialize single-flavour part to true, the rest to false
		std::vector<MatrixXb> kxab;
		for( int spin = 0; spin <=1 ; spin++ ) {
			kxab.push_back(MatrixXb::Constant( V, Nf*(Nf+1)/2, false ));
			kxab.back().leftCols(Nf) = MatrixXb::Constant( V, Nf, true);
		}

		double av_k = 0.;
		double accrate = 0;
		Complex detOld = calcDet(slacNf, kxab[0], kxab[1], Nf);

		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				// draw random point and 2 flavours
				int x = intV_dist(gen);
				//				int mu = int2mu_dist(gen);
				int a = intNf_dist(gen);
				int b = intNf_dist(gen);
				int spin = intSpin_dist(gen);
				//				int y = lat.getNeighbours(x)[mu];

				//				double factor = 1.;
				int dk = 0;		// change in kxab
				//				int dkt = 0;	// change in tilde kxab
				// calculate old value of na

				Eigen::VectorXi kxSpinSum = kxab[0].leftCols(Nf).cast<int>() + kxab[1].leftCols(Nf).cast<int>();
				int n1 = 0;
				int n2 = 0;

				n1 = -(kxSpinSum.array() == 1).count();
				n2 = -(kxSpinSum.array() == 2).count();

//				for( int n = 0; n <= Nf; n++ ) {
//					na(n) = -(kx.array() == n).count();
//				}

				//			std::cout << "Constraint before update: " << std::endl;
				//			constraintViolated( kxab, Nf, x );
				//			constraintViolated( kxab, Nf, y );
				//			std::cout << "updating..." << std::endl;

				std::vector<MatrixXb> kxabOld(kxab);

				// update at the drawn index
				dk = -kxab[0].count() - kxab[1].count();
				kxab[spin]( x, linIndex( Nf, a, b ) ) = !kxab[spin]( x, linIndex( Nf, a, b ) );
				dk += kxab[spin].count();

				if( a == b ) {
					// updating with same flavour,
					// choose, if we update another spin (setting a kxaa=2) or another point (setting two kxaa=1, keeping n1 even)

					if( 0.5 > uni_real_dist(gen) ) {
						kxab[ 1-spin ]( x, a ) = !kxab[ 1-spin ]( x, a );
					} else {
						int mu = int2mu_dist(gen);
						int y = lat.getNeighbours( x )[mu];
						kxab[spin]( y, a ) = !kxab[spin]( y, a );
						if( constraintViolated( kxab[0], kxab[1], Nf, x ) ) {
							kxab = kxabOld;
//							std::cout << "violated! Reset and continue loop..." << std::endl;
							continue;
						}
					}
				} else {
					// updating with two different flavours, enforce kxab = kxba
					if( dk < 0 ) {
						// we have to delete another spin in kxba, choose one randomly and check,
						// if it is set, otherwise the other one must be set, since we have only two spins.
						int spin2 = intSpin_dist(gen);
						if( kxab[spin2]( x, linIndex( Nf, b, a ) ) ) {
							kxab[spin2]( x, linIndex( Nf, b, a ) ) = ! kxab[spin2]( x, linIndex( Nf, b, a ) );
						} else {
							kxab[1-spin2]( x, linIndex( Nf, b, a ) ) = ! kxab[1-spin2]( x, linIndex( Nf, b, a ) );
						}
					} else {
						// we have to set another spin in kxba as before
						int spin2 = intSpin_dist(gen);
						if( !kxab[spin2]( x, linIndex( Nf, b, a ) ) ) {
							kxab[spin2]( x, linIndex( Nf, b, a ) ) = ! kxab[spin2]( x, linIndex( Nf, b, a ) );
						} else {
							kxab[1-spin2]( x, linIndex( Nf, b, a ) ) = ! kxab[1-spin2]( x, linIndex( Nf, b, a ) );
						}
					}

					//					if( 0.5 > uni_real_dist(gen) ) {
					//						kxab(x,linIndex(Nf,a, b)) = !kxab(x,linIndex(Nf,a, b));
					//						kxab(y,linIndex(Nf,a, b)) = !kxab(y,linIndex(Nf,a, b));
					//					} else {
					//						kxab(x,linIndex(Nf,a, b)) = !kxab(x,linIndex(Nf,a, b));
					//						kxab(x,a) = !kxab(x,a);
					//						kxab(x,b) = !kxab(x,b);
					//					}
					//					factor *= pow( 2., double(dk) );
				}
				dk += kxab[0].count() + kxab[1].count();

				if( constraintViolated( kxab[0], kxab[1], Nf, x ) ) {
					kxab = kxabOld;
					//				std::cout << "violated! Reset and continue loop..." << std::endl;
					continue;
				}
				Complex det = calcDet(slacNf, kxab[0], kxab[1], Nf);

				// TODO: calculate change in na

//				kxSpinSum = kxab[0].leftCols(Nf).cast<int>() + kxab[1].leftCols(Nf).cast<int>();
//				n1 += (kxSpinSum.array() == 1).count();
//				n2 += (kxSpinSum.array() == 2).count();

//				kx0 = kxab[0].leftCols(Nf).rowwise().count().cast<int>();
//				for( int n = 0; n <= Nf; n++ ) {
//					na(n) += (kx.array() == n).count();
//					factor *= std::pow( getHypergeometricFactor( n ), double(na(n)) );
//					//					std::cout << "n=" << n << "\tfactor=" << factor << "\thyperFactor=" << getHypergeometricFactor(n) << std::endl;
//				}

				//			std::cout << slacNf[0].getMatrix() << std::endl << std::endl;
				//			std::cout << "kxab=" << std::endl << kxab << std::endl << std::endl;
				//			std::cout << "kx=" << std::endl << kx << std::endl;
				//			std::cout << "na=" << std::endl << na << std::endl;
				//			std::cout << "x=" << x << "\ty=" << y << "\ta=" << a << "\tb=" << b;
				//			std::cout << "\tka=" << kxab.leftCols(Nf).count() << "\tkab="<< kxab.rightCols( Nf*(Nf-1)/2 ).count() << "\tdk=" << dk;
				//			std::cout << "\tdetOld: " << detOld << "\tdet: " << det;
				double factor = 1.;	//TODO: implement correct factor for Nf > 1
				double dw = std::pow(kappa, dk);
				Complex weight =  factor*dw*(det/detOld);

				double r = uni_real_dist(gen);
				//			bool accepted = false;
				if( std::fabs(weight) > r ) {
					//				accepted = true;
					accrate++;
					detOld = det;
				} else {
					kxab = kxabOld;
				}
				//			std::cout << "\tdw: " << dw << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;
			}
			if( measure >= numThermal ) av_k += ( kxab[0].count() + kxab[1].count() )/double(V);
		}
		av_k/=double(numMeasures);
		accrate /= double((numMeasures+numThermal)*upPerMeasure);
		std::cout << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
