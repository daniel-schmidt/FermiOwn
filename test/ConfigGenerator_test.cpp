/*
 * ConfigGenerator_test.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>

#include "ConfigPerPointGeneratorTh.h"
#include "ConfigPerPointGeneratorGN.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 3 ) {
		std::cerr << "Need 2 arguments: spinor size and number of flavours" << std::endl;
		exit(1);
	}

	size_t numSpins = atoi( argv[1] );
	size_t numFlavours = atoi( argv[2] );

	std::ranlux48 rndGen;
	ConfigPerPointGeneratorTh Thgen( numSpins, numFlavours, &rndGen );

	Thgen.generateAllowedConfs();

	MatrixXb allConfs = Thgen.getAllConfs();
	std::cout << "The Thirring model with Nf=" << numFlavours << " has " << allConfs.rows() << " allowed configurations:" << std::endl;
	std::cout << allConfs << std::endl << std::endl;

	std::cout << "Generating 10 random configurations from list of all configs: " << std::endl;
	for( int i = 0; i < 10; i++ ) {
		std::cout << Thgen.getRandomConf() << std::endl;
	}
	std::cout << std::endl;

	ConfigPerPointGeneratorGN GNgen( numSpins, numFlavours, &rndGen );
	GNgen.generateAllowedConfs();

	allConfs = GNgen.getAllConfs();
	std::cout << "All " << allConfs.rows() << " allowed configurations for the Gross-Neveu model:" << std::endl;
	std::cout << allConfs << std::endl << std::endl;

	std::cout << "Generating 10 random configurations from list of all configs: " << std::endl;
	for( int i = 0; i < 10; i++ ) {
		std::cout << GNgen.getRandomConf() << std::endl;
	}
}
