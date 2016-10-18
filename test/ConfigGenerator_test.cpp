/*
 * ConfigGenerator_test.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>

#include "../src/Actions/ConfigPerPointGenerator.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 3 ) {
		std::cerr << "Need 2 arguments: spinor size and number of flavours" << std::endl;
		exit(1);
	}

	size_t numSpins = atoi( argv[1] );
	size_t numFlavours = atoi( argv[2] );

	std::ranlux48 rndGen;
	ConfigPerPointGenerator gen( numSpins, numFlavours, &rndGen );

	gen.generateAllowedConfs();

	MatrixXb allConfs = gen.getAllConfs();
	std::cout << "Nf=" << numFlavours << " has " << allConfs.rows() << " allowed configurations:" << std::endl;
	std::cout << allConfs << std::endl << std::endl;

	std::cout << "Generating 10 random configurations from list of all configs: " << std::endl;
	for( int i = 0; i < 10; i++ ) {
		std::cout << gen.getRandomConf() << std::endl;
	}

}
