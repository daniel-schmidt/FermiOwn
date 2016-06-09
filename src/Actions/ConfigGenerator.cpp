/*
 * ConfigGenerator.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include "ConfigGenerator.h"

namespace FermiOwn {

ConfigGenerator::ConfigGenerator( const size_t numSpins, const size_t numFlavours ) :
	Nf( numFlavours ),
	dimSpinor( numSpins )
{
	if( dimSpinor != 2 ) std::cout << "Warning: ConfigGenerator called with spinor size != 2, this may not be implemented correctly!" << std::endl;
}

ConfigGenerator::~ConfigGenerator() {}

void  ConfigGenerator::generateAllowedConfs() {

	size_t numCols = (Nf*Nf)*dimSpinor;
	size_t numConfigs = 1;
	numConfigs = numConfigs << numCols;

	//TODO: implement something better, which doesn't need to allocate so much memory...
	try {
		allowedConfs = Eigen::MatrixXi::Zero( numConfigs, numCols );
	} catch ( const std::bad_alloc& ) {
		std::cerr << "ERROR: Not enough memory!" << std::endl;
		std::cerr << "ConfigGenerator tries to allocate an array of size " << numConfigs << " x " << numCols << std::endl;
		exit(1);
	}
	size_t numAllowedConfs = 0;
	for( size_t conf = 0; conf < numConfigs; conf++ ) {
		size_t bits = conf;
		size_t bits2 = conf;
		FieldBoolean fConf( 1, dimSpinor, Nf, NULL, zeroInit );

		for( size_t spin = 0; spin < 2; spin++ ) {
			for( size_t a = 0; a < Nf; a++ ) {
				for( size_t b = 0; b < Nf; b++ ) {
					if( bits2 % 2 == 1 ) {
						fConf.setValue( true, 0, spin, a, b );
					}
					bits2 = bits2 >> 1;
				}
			}
		}

		if( !fConf.constraintViolated(0) ) {
			for( size_t bit = 0; bit < numCols; bit++ ) {
				if( bits % 2 == 1 ) {
					allowedConfs( numAllowedConfs, numCols-1-bit ) = 1;
				}
				bits = bits >> 1;
			}
			numAllowedConfs ++;
		}
	}
	allowedConfs.conservativeResize( numAllowedConfs, Eigen::NoChange );
}

} /* namespace FermiOwn */
