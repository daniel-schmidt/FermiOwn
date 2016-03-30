/*
 * HMC_test.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "Action.h"
#include "HMC.h"

#include "Field.h"
#include "Lattice.h"

int main() {
	using namespace FermiOwn;
	std::cout << "Testing the HMC algorithm..." << std::endl;

	// initialize random number generator
	std::ranlux48 rndGen;

	// initialize lattice
	size_t Nt = 4, Ns = 4, dim = 3;
	Lattice lat(Nt, Ns, dim);

	// initialize action
	double kappa = 0.10, lambda = 1.145;
	Action act(lat, kappa, lambda);

	// initialize field
	Field<Real> phi(lat.getVol(), 1, &rndGen, gaussianInit);

	double HMCt = 1.;
	size_t HMCnt = 30;
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
