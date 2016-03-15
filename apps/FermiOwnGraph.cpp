/*
 * FermiOwnGraph.cpp
 *
 *  Created on: 15.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "SlacOperatorMatrix.h"
#include "Lattice.h"

namespace FermiOwn {
typedef Eigen::Matrix< bool, 1, Eigen::Dynamic > VectorXb;

Complex calcDet( SlacOperatorMatrix& dslac, VectorXb kx ) {
	for( int x = 0; x < kx.size(); x++ ) {
		if( kx(x) ) {
			dslac.deletePoint(x);
		}
	}
	return pow(I,dslac.getMatrix().size())*dslac.det();
}

}

int main( int argc, char** argv ) {
	using namespace Eigen;
	using namespace FermiOwn;


	if( argc != 5 ) {
		std::cerr << "Need 4 arguments: lattice size, numMeasures, upPerMeasure, coupling" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1] );
	size_t V = N*(N-1)*(N-1);
	int numMeasures = atoi( argv[2] );
	int upPerMeasure = atoi( argv[3] );
	size_t dim = 3;
	double lambda = atof( argv[4] );

	Lattice lat( N, N-1, dim );
	SlacOperatorMatrix dslac3d( N, N-1, dim );
	std::ranlux48 gen;
	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist( 0, V-1 );
	std::uniform_int_distribution<int> int2mu_dist( 0, 2*dim-1 );


	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
		VectorXb kx = VectorXb::Constant( V, true );
		double av_k = 0.;
		double accrate = 0;
		for( int measure = 0; measure < numMeasures; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				int x = intV_dist(gen);
				int mu = int2mu_dist(gen);
				int y = lat.getNeighbours(x)[mu];
				kx(x) = !kx(x);
				kx(y) = !kx(y);
				int k = kx.count();

//				std::cout << "x=" << x << "\ty=" << y << "\tk=" << k;

				Complex det = calcDet(dslac3d, kx);
				double weight = std::pow(-1.5/lambda, k) * real(det);

				double r = uni_real_dist(gen);
				bool accepted = false;
				if( weight > r ) {
					accepted = true;
					accrate++;
				} else {
					kx(x) = !kx(x);
					kx(y) = !kx(y);
				}

//				std::cout << "\tdet: " << det << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;
				dslac3d.setFull();
			}
			av_k += kx.count();
		}
		av_k/=numMeasures*V;
		accrate /= numMeasures*upPerMeasure;
		std::cout << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
