/*
 * ConfigGenerator_test.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>
#include "ConfigGenerator.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	if( argc != 3 ) {
		std::cerr << "Need 2 arguments: spinor size and number of flavours" << std::endl;
		exit(1);
	}

	size_t numSpins = atoi( argv[1] );
	size_t numFlavours = atoi( argv[2] );

	ConfigGenerator gen( numSpins, numFlavours );

	gen.generateAllowedConfs();

	MatrixXb allConfs = gen.getAllConfs();
	std::cout << "Nf=" << numFlavours << " has " << allConfs.rows() << " allowed configurations:" << std::endl;
	std::cout << allConfs << std::endl;

}
