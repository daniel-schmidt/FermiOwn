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
#include "Timer.h"
#include "Lattice.h"
#include "ThirringKField.h"
#include "GrossNeveuKField.h"
#include "FermiBoolMetropolis.h"
#include "ConfigGenerator.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 10 ) {
		std::cerr << "Need 9 arguments: Nt, Ns, numThermal, numMeasures, upPerMeasure, coupling, dcoupling, Nf, model name (Th or GN)" << std::endl;
		exit(1);
	}

	Timer totalTimer;
	totalTimer.start();

	// read parameters from console
	size_t dim = 3;
	size_t Nt = atoi( argv[1] );
	size_t Ns = atoi( argv[2] );
	size_t V = Nt*Ns*Ns;
	int numThermal = atoi( argv[3] );
	int numMeasures = atoi( argv[4] );
	int upPerMeasure = atoi( argv[5] );
	double lambdaInitial = atof( argv[6] );
	double dLambda = atof(argv[7]);
	size_t Nf = atof( argv[8] );
	size_t numSpin = 2;

	char* model = argv[9];
	// initialize classes
	Lattice lat( Nt, Ns, dim );
	Timer singleRunTimer;
	std::ofstream av_k_file("avk.dat", std::ofstream::app);
	av_k_file << "#lambda" << "\t" << "av_k" << "\t" << "accrate" << "\t" << "real(expPhase)" << "\t" << "imag(expPhase)" << std::endl;
#pragma omp parallel for
	for( int step = 0; step < 10; step++ ) {
		double lambda = lambdaInitial + (step+1)*dLambda;

		std::ranlux48 gen;
		std::ofstream kfile( "avk" + std::to_string(lambda) + ".dat" );

		// initialize single-flavour part to true, the rest to false
		BasicKField* kfield;
		if( model == NULL ) {
			std::cerr << "You did not pass a model: argument Th or GN missing." << std::endl;
			exit(1);
		} if( std::strcmp( model, "Th" ) == 0 ) {
			kfield = new ThirringKField( lat.getVol(), {numSpin, Nf, Nf} );
		} else if( std::strcmp( model, "GN" ) == 0 ) {
			kfield = new GrossNeveuKField( lat.getVol(), {numSpin, Nf} );
		} else {
			std::cerr << "Model not known, specify either Th or GN as option!" << std::endl;
			exit(1);
		}
//		for( int a = 0; a < Nf; a++ ) {
//			for( int spin = 0; spin < numSpin; spin++ ) {
//				for( size_t x = 0; x < lat.getVol(); x++ ) {
//					kxiab.setValue( true, x, spin, a, a );
//				}
//			}
//		}

		FermiBoolMetropolis updater( *kfield, lat, lambda/2., Nf, &gen ); //TODO: we are simulating with lambda/2 due to convention in Base...

		double av_k = 0.;
		auto measure = [&](){
			double k = kfield->sumAll()/double(2*V);
			av_k += k;
			kfile << k << "\t" << std::endl;
		};

		ConfigGenerator confGen( numThermal, numMeasures, upPerMeasure, &updater, measure );

		singleRunTimer.start();
		confGen.run();
		singleRunTimer.stop();
		singleRunTimer.printDuration("Creating configs for current parameter");
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
		Complex expPhase = updater.getAveragePhase();
		kfile.close();
#pragma omp critical
		av_k_file << lambda << "\t" << av_k << "\t" << accrate << "\t" << std::real(expPhase) << "\t" << std::imag(expPhase) << std::endl;
		delete kfield;
	}
	av_k_file.close();
	totalTimer.stop();
	totalTimer.printDuration("Total execution");
}
