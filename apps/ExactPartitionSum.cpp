/*
 * ExactPartitionSum.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>
#include "FieldBoolean.h"
#include "ConfigGenerator.h"
#include "WeightFunction.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	// TODO: even 2x3x3 for Nf=2 is too large to fit in the memory...
	size_t Nt = 2;
	size_t Ns = 3;
	size_t dim = 3;
	size_t latVol = Nt*Ns*Ns;
	size_t numSpins = 2;
	size_t Nf = 1;

	std::ofstream fout( "exactPartitionSum.dat");

	double dlambda = 0.1;
	for( int i = 1; i <= 20; i++ ) {
		double lambda = dlambda*i;
		FieldBoolean kxiab( latVol, numSpins, Nf, NULL, zeroInit );
		ConfigGenerator confGen( numSpins, Nf );
		WeightFunction weight( kxiab, Nt, Ns, dim, Nf, lambda );

		confGen.generateAllowedConfs();
		MatrixXb allowedConfs = confGen.getAllConfs();

		size_t confsPerX = allowedConfs.rows();
		//TODO: this should check somehow, if all configs fit in memory...
		long double numConfigs = std::pow( confsPerX, latVol );
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
					for( size_t a = 0; a < Nf; a++ ) {
						for( size_t b = 0; b < Nf; b++ ) {
							kxiab.setValue( bool(allowedConfs( newConfIndex, spin*Nf*Nf + b*Nf + a)), x, spin, a, b );
						}
					}
				}
			}
			//		kxiab.Print();

			sum += weight.calculateWeight();
			weight.reset();
		}
		fout << lambda << "\t" << sum << std::endl;
	}
	fout.close();
}
