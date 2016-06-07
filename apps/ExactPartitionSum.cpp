/*
 * ExactPartitionSum.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include <Eigen/Dense>
#include "FieldBoolean.h"
#include "ConfigGenerator.h"
#include "WeightFunction.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	size_t latVol = 4;
	size_t numSpins = 2;
	size_t numFlavours = 2;

	FieldBoolean kxiab( latVol, numSpins, numFlavours, NULL, zeroInit );
	ConfigGenerator confGen( latVol, numSpins, numFlavours );
	WeightFunction weight( kxiab );

	confGen.generateAllowedConfs();
	Eigen::MatrixXi allowedConfs = confGen.getAllConfs();

	size_t confsPerX = allowedConfs.rows();
	//TODO: this should check somehow, if all configs fit in memory...
	size_t numConfigs = std::pow( confsPerX, latVol );
//	std::cout << "numConfigs: " << numConfigs << std::endl;

	Complex sum = 0.;

	for( size_t conf = 0; conf < numConfigs; conf++ ) {
		size_t confNum = conf;
//		std::cout << "Conf " << conf << " confNum: ";

		Eigen::VectorXi confList = Eigen::VectorXi::Zero( latVol );
		size_t count = 0;
		while( confNum >= confsPerX ) {
			confList( count ) = confNum % confsPerX;
			confNum /= confsPerX;
			count ++;
		}

		confList( count ) = confNum;
//		std::cout << std::endl<< confList <<std::endl;

		for( size_t x = 0; x < latVol; x++ ) {
			size_t newConfIndex = confList( x );
			for( size_t spin = 0; spin < 2; spin++ ) {
				for( size_t a = 0; a < numFlavours; a++ ) {
					for( size_t b = 0; b < numFlavours; b++ ) {
						kxiab.setValue( bool(allowedConfs( newConfIndex, spin*numFlavours*numFlavours + b*numFlavours + a)), x, spin, a, b );
					}
				}
			}
		}
//		kxiab.Print();

		sum += weight.calculateWeight();
	}
	std::cout << 2./kappa << "\t" << sum << std::endl;
}
