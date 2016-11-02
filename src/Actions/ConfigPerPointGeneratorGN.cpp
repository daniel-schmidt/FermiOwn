/*
 * ConfigPerPointGeneratorGN.cpp
 *
 *  Created on: 02.11.2016
 *      Author: dschmidt
 */

#include "ConfigPerPointGeneratorGN.h"

namespace FermiOwn {

ConfigPerPointGeneratorGN::ConfigPerPointGeneratorGN( const size_t numSpins, const size_t numFlavours, std::ranlux48* randomGenerator) :
	ConfigPerPointGenerator( numSpins, numFlavours, randomGenerator)
{}

ConfigPerPointGeneratorGN::~ConfigPerPointGeneratorGN() {}

void  ConfigPerPointGeneratorGN::generateAllowedConfs() {

	size_t numCols = Nf*dimSpinor;
	size_t numConfigs = 1;
	numConfigs = numConfigs << numCols;

	//TODO: implement something better, which doesn't need to allocate so much memory...
	try {
		allowedConfs = MatrixXb::Zero( numConfigs, numCols );
	} catch ( const std::bad_alloc& ) {
		std::cerr << "ERROR: Not enough memory!" << std::endl;
		std::cerr << "ConfigGenerator tries to allocate an array of size " << numConfigs << " x " << numCols << std::endl;
		exit(1);
	}
	size_t numAllowedConfs = 0;
	MatrixXb allConfs = allowedConfs;
	for( size_t conf = 0; conf < numConfigs; conf++ ) {

		// convert conf number to binary and count number of 1s, this must be even.
		std::bitset<sizeof(size_t) * CHAR_BIT> binary( conf );

		if( binary.count()%2 == 0 ) {
			for( size_t bitCount = 0; bitCount < numCols; bitCount++ ) {
				allowedConfs( numAllowedConfs, bitCount ) = binary[ bitCount ];
			}
			numAllowedConfs ++;
		}
//		size_t bits = conf;
//		size_t bits2 = conf;
//		FieldBoolean fConf( 1, dimSpinor, Nf, NULL, zeroInit );
//
//		for( size_t spin = 0; spin < 2; spin++ ) {
//			for( size_t a = 0; a < Nf; a++ ) {
//				if( bits2 % 2 == 1 ) {
//					fConf.setValue( true, 0, spin, a, b );
//				}
//				bits2 = bits2 >> 1;
//			}
//		}
//
//		if( !fConf.constraintViolated(0) ) {
//			for( size_t bit = 0; bit < numCols; bit++ ) {
//				if( bits % 2 == 1 ) {
//					allowedConfs( numAllowedConfs, numCols-1-bit ) = true;
//				}
//				bits = bits >> 1;
//			}
//			numAllowedConfs ++;
//		}
	}
	allowedConfs.conservativeResize( numAllowedConfs, Eigen::NoChange );
	int_dist = std::uniform_int_distribution<int>( 0, numAllowedConfs-1 );
}

} /* namespace FermiOwn */
