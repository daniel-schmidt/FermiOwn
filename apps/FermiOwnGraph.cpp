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
typedef Eigen::Matrix< bool, 1, Eigen::Dynamic > VectorXb;
typedef Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > MatrixXb;

Complex calcDet( SlacOperatorMatrix& dslac, const VectorXb& kx ) {
	for( int x = 0; x < kx.size(); x++ ) {
		if( kx(x) ) {
			dslac.deletePoint(x);
		}
	}
	return pow(I,dslac.getMatrix().rows())*dslac.det();
}

Complex calcDet( std::vector<SlacOperatorMatrix>& dslac, const MatrixXb& kxab ) {
	Complex det = 1.;
	for( int a = 0; a < kxab.cols(); a++ ) {
		for( int x = 0; x < kxab.rows(); x++ ) {
			if( kxab( a, x ) ) {
				dslac[a].deletePoint(x);
			}
		}
		det *= pow(I,dslac[a].getMatrix().rows())*dslac[a].det();
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

}

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 7 ) {
		std::cerr << "Need 6 arguments: lattice size, numThermal, numMeasures, upPerMeasure, coupling, Nf" << std::endl;
		exit(1);
	}

	// read parameters from console
	size_t N = atoi( argv[1] );
	size_t V = N*(N-1)*(N-1);
	int numThermal = atoi( argv[2] );
	int numMeasures = atoi( argv[3] );
	int upPerMeasure = atoi( argv[4] );
	size_t dim = 3;
	double lambda = atof( argv[5] );
	int Nf = atof( argv[6] );

	// initialize classes
	Lattice lat( N, N-1, dim );
	std::vector<SlacOperatorMatrix> slacNf;
	for( int a = 0; a < Nf; a++) {
		SlacOperatorMatrix dslac3d( N, N-1, dim );
		slacNf.push_back(dslac3d);
	}

	// initialize random number generator and distributions
	std::ranlux48 gen;
	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist( 0, V-1 );
	std::uniform_int_distribution<int> int2mu_dist( 0, 2*dim-1 );
	std::uniform_int_distribution<int> intNf_dist( 0, Nf-1 );

	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
		double kappa = 2./lambda;


//		VectorXb kx = VectorXb::Constant( V, true );
		// Data layout: first Nf columns hold kxaa for a=1,..,Nf, then kxab with a<b orderd (1,2), (1,3), (2,3), ...
		MatrixXb kxab = MatrixXb::Constant( V, Nf*(Nf+1)/2, true );

		double av_k = 0.;
		double accrate = 0;
		Complex detOld = calcDet(slacNf[0], kxab);
		int count = 0;

		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				// draw random point, direction and 2 flavours and get neighbour
				int x = intV_dist(gen);
				int mu = int2mu_dist(gen);
				int a = intNf_dist(gen);
				int b = intNf_dist(gen);
				int y = lat.getNeighbours(x)[mu];

				// calculate ks
				VectorXb kx = kxab.leftCols(Nf).rowwise().count();
				VectorXb ktx = kxab.rightCols( Nf*(Nf-1)/2 ).rowwise().count();

				int dk = 0;		// change in kxaa
				int dkt = 0;	// change in kxab

				if( a == b ) {
					// switching link with same flavour
					kxab(a,x)
					dk = -kx.count();
					kx(x) = !kx(x);
					kx(y) = !kx(y);
					int k = kx.count();
					dk += k;

				} else {
					// switching flavour at
				}

				Complex det = calcDet(slacNf[0], kx);
//				std::cout << "x=" << x << "\ty=" << y << "\tk=" << k << "\tdk=" << dk << "\tdet: " << detOld-det;
				double dw = std::pow(-1.5*kappa, dk);
				Complex weight =  dw*(det/detOld);

				double r = uni_real_dist(gen);

				if( std::fabs(weight) > r ) {
					accrate++;
					detOld = det;
				} else {
					kx(x) = !kx(x);
					kx(y) = !kx(y);
				}
//				std::cout << "\tdw: " << dw << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;

				count ++;
				slacNf[0].setFull();
			}
			if( measure >= numThermal ) av_k += kx.count()/double(V);
		}
		av_k/=double(numMeasures);
		accrate /= double((numMeasures+numThermal)*upPerMeasure);
		std::cout << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
