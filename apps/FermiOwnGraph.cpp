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
#include "FermiBoolMetropolis.h"

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


	for( int step = 0; step < 10; step++ ) {
		lambda += 0.1;
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
		FermiBoolMetropolis updater( kxiab, slac, lat, lambda, Nf, &gen );

		//		updater.sumAllConfs();

		double av_k = 0.;

		for( int measure = 0; measure < numMeasures+numThermal; measure++) {
			for( int i = 0; i < upPerMeasure; i++ ) {
				updater.updateField();
			}
			if( measure >= numThermal )
			{
				double k = kxiab.sumAll()/double(2*V);
				av_k += k;
				kfile << k << "\t" << std::endl;
			}
		}
		av_k/=double(numMeasures);
		double accrate = updater.acceptanceCounter/double((numMeasures+numThermal)*upPerMeasure);
		std::cerr << lambda << "\t" << av_k << "\t" << accrate << std::endl;
		kfile.close();
	}
}
