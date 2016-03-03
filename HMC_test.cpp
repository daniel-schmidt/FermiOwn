/*
 * HMC_test.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "FieldScalar.h"
#include "Action.h"
#include "Lattice.h"
#include "HMC.h"

int main() {
	std::cout << "Testing the HMC algorithm..." << std::endl;

	// initialize random number generator
	std::ranlux48 rndGen;

	// initialize lattice
	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	// initialize action
	double kappa = 0.5, lambda = 0.1;
	Action act(lat, kappa, lambda);

	// initialize field
	FieldScalar<Real> phi(lat.getVol(), &rndGen, zeroInit);

	double HMCt = 0.6;
	size_t HMCnt = 8;
	// initialize HMC
	HMC updater( HMCt, HMCnt, act, phi, &rndGen );

	size_t updates = 100;
	size_t stepsPerUpdate=1;

	size_t acceptance = 0.;
	for( size_t u = 0; u<updates; u++ ) {
		for( size_t spu = 0; spu < stepsPerUpdate; spu++ )
		{
			bool accepted = updater.update();
			if(accepted) acceptance++;
		}
	}
	double accRate = acceptance/double(updates*stepsPerUpdate);

	std::cout << "Accepted "<< accRate << std::endl;
	phi.Print();
}
