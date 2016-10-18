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
#include "FermiBoolMetropolis.h"
#include "ConfigGenerator.h"

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
	double lambdaInitial = atof( argv[6] );
	double dLambda = 0.1;
	int Nf = atof( argv[7] );
	int numSpin = 2;
	// initialize classes
	Lattice lat( Nt, Ns, dim );

	// initialize random number generator and distributions

#pragma omp parallel for
	for( int step = 0; step < 10; step++ ) {
		double lambda = lambdaInitial + (step+1)*dLambda;
		std::ranlux48 gen;
		std::ofstream kfile( "avk" + std::to_string(lambda) + ".dat" );

		// initialize single-flavour part to true, the rest to false
		FieldBoolean kxiab( lat.getVol(), numSpin, Nf, &gen, zeroInit );
//		for( int a = 0; a < Nf; a++ ) {
//			for( int spin = 0; spin < numSpin; spin++ ) {
//				for( size_t x = 0; x < lat.getVol(); x++ ) {
//					kxiab.setValue( true, x, spin, a, a );
//				}
//			}
//		}

		FermiBoolMetropolis updater( kxiab, lat, lambda/2., Nf, &gen ); //TODO: we are simulating with lambda/2 due to convention in Base...

		double av_k = 0.;
		auto measure = [&](){
			double k = kxiab.sumAll()/double(2*V);
			av_k += k;
			kfile << k << "\t" << std::endl;
		};

		ConfigGenerator confGen( numThermal, numMeasures, upPerMeasure, &updater, measure );
		confGen.run();

//		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
//			for( int i = 0; i < upPerMeasure; i++ ) {
//				updater.step();
//			}
//			if( measure >= numThermal )
//			{
//				double k = kxiab.sumAll()/double(2*V);
//				av_k += k;
//				kfile << k << "\t" << std::endl;
//			}
//		}
		av_k/=double(numMeasures);
		double accrate = updater.getAcceptance();
		Complex expPhase = 0;
		std::cerr << lambda << "\t" << av_k << "\t" << accrate << "\t" << std::real(expPhase) << "\t" << std::imag(expPhase) << std::endl;
		kfile.close();
	}
}
