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
	return pow(I,dslac.getMatrix().rows())*dslac.det();
}

}

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 6 ) {
		std::cerr << "Need 5 arguments: lattice size, numThermal, numMeasures, upPerMeasure, coupling" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1] );
	size_t V = N*(N-1)*(N-1);
	int numThermal = atoi( argv[2] );
	int numMeasures = atoi( argv[3] );
	int upPerMeasure = atoi( argv[4] );
	size_t dim = 3;
	double lambda = atof( argv[5] );

	Lattice lat( N, N-1, dim );
	SlacOperatorMatrix dslac3d( N, N-1, dim );
	std::ranlux48 gen;
	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist( 0, V-1 );
	std::uniform_int_distribution<int> int2mu_dist( 0, 2*dim-1 );


	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
		double kappa = 2./lambda;
		VectorXb kx = VectorXb::Constant( V, true );
		double av_k = 0.;
		double accrate = 0;
		Complex detOld = calcDet(dslac3d, kx);
		int count = 0;
		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				int x = intV_dist(gen);
				int mu = int2mu_dist(gen);
				int y = lat.getNeighbours(x)[mu];
				int dk = -kx.count();
				kx(x) = !kx(x);
				kx(y) = !kx(y);
				int k = kx.count();
				dk += k;

				Complex det = calcDet(dslac3d, kx);
//				std::cout << "x=" << x << "\ty=" << y << "\tk=" << k << "\tdk=" << dk << "\tdet: " << detOld-det;
				double dw = std::pow(-1.5*kappa, dk);
				Complex weight =  dw*(det/detOld);

				double r = uni_real_dist(gen);
				bool accepted = false;
				if( std::fabs(weight) > r ) {
					accepted = true;
					accrate++;
					detOld = det;
				} else {
					kx(x) = !kx(x);
					kx(y) = !kx(y);
				}
//				std::cout << "\tdw: " << dw << "\tweight: " << weight << "\taccepted: " << accepted << std::endl;

				count ++;
				dslac3d.setFull();
			}
			if( measure >= numThermal ) av_k += kx.count()/double(V);
		}
		av_k/=double(numMeasures);
		accrate /= double((numMeasures+numThermal)*upPerMeasure);
		std::cout << lambda << "\t" << av_k << "\t" << accrate << std::endl;
	}
}
