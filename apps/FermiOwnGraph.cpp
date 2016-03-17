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
typedef Eigen::Matrix< bool, Eigen::Dynamic, 1 > VectorXb;
typedef Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > MatrixXb;

Complex calcDet( SlacOperatorMatrix& dslac, const VectorXb& kx ) {
	for( int x = 0; x < kx.size(); x++ ) {
		if( kx(x) ) {
			dslac.deletePoint(x);
		}
	}
	return pow(I,dslac.getMatrix().rows())*dslac.det();
}

Complex calcDet( std::vector<SlacOperatorMatrix> dslac, const MatrixXb& kxab ) {
	Complex det = 1.;
	for( size_t a = 0; a < dslac.size(); a++ ) {
		//		for( int x = 0; x < kxab.rows(); x++ ) {
		//			if( kxab( x, a ) ) {
		//				dslac[a].deletePoint(x);
		//			}
		//		}
		//		det *= pow(I,dslac[a].getMatrix().rows())*dslac[a].det();
		det *= calcDet( dslac[a], kxab.col(a) );
	}
	return det;
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

	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
		double kappa = 2./lambda;

		// Data layout: first Nf columns hold kxaa for a=1,..,Nf, then kxab with a<b orderd (1,2), (1,3), (2,3), ...

		//		std::cout << "Nf*(Nf+1)/2="<<Nf*(Nf+1)/2 << std::endl;
		MatrixXb kxab = MatrixXb::Constant( V, Nf*(Nf+1)/2, true );

		double av_k = 0.;
		double accrate = 0;
		Complex detOld = calcDet(slacNf, kxab);

		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				// draw random point, direction and 2 flavours and get neighbour
				int x = intV_dist(gen);
				int mu = int2mu_dist(gen);
				int a = intNf_dist(gen);
				int b = intNf_dist(gen);
				int y = lat.getNeighbours(x)[mu];

				double factor = 1.;
				int dk = 0;		// change in kxab

				Eigen::VectorXi kx = kxab.leftCols(Nf).rowwise().count().cast<int>();
				Eigen::VectorXi na(Nf+1);
				for( int n = 0; n <= Nf; n++ ) {
					na(n) = -(kx.array() == n).count();
				}

				if( a == b ) {
					// switching link with same flavour
					dk = -kxab.leftCols(Nf).count();
					kxab(x,linIndex(Nf,a, b)) = !kxab(x,linIndex(Nf,a, b));
					kxab(y,linIndex(Nf,a, b)) = !kxab(y,linIndex(Nf,a, b));
					dk += kxab.leftCols(Nf).count();
				} else {
					// switching flavour at two positions
					dk = -kxab.rightCols( Nf*(Nf-1)/2 ).count();
					kxab(x,linIndex(Nf,a, b)) = !kxab(x,linIndex(Nf,a, b));
					kxab(y,linIndex(Nf,a, b)) = !kxab(y,linIndex(Nf,a, b));
					dk += kxab.rightCols( Nf*(Nf-1)/2 ).count();
					dk *= 2;
					//					factor = sqrt(2);
				}

				Complex det = calcDet(slacNf, kxab);

				kx = kxab.leftCols(Nf).rowwise().count().cast<int>();
				for( int n = 0; n <= Nf; n++ ) {
					na(n) += (kx.array() == n).count();
					factor *= std::pow( getHypergeometricFactor( n ), double(na(n)) );
//					std::cout << "n=" << n << "\tfactor=" << factor << "\thyperFactor=" << getHypergeometricFactor(n) << std::endl;
				}

				//				std::cout << slacNf[0].getMatrix() << std::endl << std::endl;
//								std::cout << "kxab=" << std::endl << kxab << std::endl << std::endl;
//								std::cout << "kx=" << std::endl << kx << std::endl;
//								std::cout << "na=" << std::endl << na << std::endl;
//								std::cout << "x=" << x << "\ty=" << y << "\ta=" << a << "\tb=" << b;
//								std::cout << "\tka=" << kxab.leftCols(Nf).count() << "\tkab="<< kxab.rightCols( Nf*(Nf-1)/2 ).count() << "\tdk=" << dk;
//								std::cout << "\tdetOld: " << detOld << "\tdet: " << det;
				double dw = std::pow(kappa, dk);
				Complex weight =  factor*dw*(det/detOld);

				double r = uni_real_dist(gen);
//				bool accepted = false;
				if( std::fabs(weight) > r ) {
//					accepted = true;
					accrate++;
					detOld = det;
				} else {
					kxab(x,linIndex(Nf,a, b)) = !kxab(x,linIndex(Nf,a, b));
					kxab(y,linIndex(Nf,a, b)) = !kxab(y,linIndex(Nf,a, b));
				}
//				std::cout << "\tdw: " << dw << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;
			}
			if( measure >= numThermal ) av_k += ( kxab.leftCols(Nf).count()+2.*kxab.rightCols( Nf*(Nf-1)/2 ).count() )/double(V);
		}
		av_k/=double(numMeasures);
		accrate /= double((numMeasures+numThermal)*upPerMeasure);
		std::cout << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
