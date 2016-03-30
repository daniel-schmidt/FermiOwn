/*
 * Metropolis_test.cpp
 *
 *  Created on: 10.03.2016
 *      Author: dschmidt
 */

#include "Action.h"
#include "Metropolis.h"

#include "Field.h"
#include "Lattice.h"

int main() {
	using namespace FermiOwn;
	std::cout << "Testing the Metropolis algorithm..." << std::endl;

	// initialize random number generator
	std::ranlux48 rndGen;

	// initialize lattice
	size_t Nt = 8, Ns = 8, dim = 3;
	Lattice lat(Nt, Ns, dim);

	// initialize action
	double lambda = 1.145;

	// initialize field
	Field<Real> phi(lat.getVol(), 1, &rndGen, gaussianInit);


	// initialize Monte Carlo Algorithm
	double delta = 0.7;
	size_t updates = 100;
	size_t stepsPerUpdate=1;

	// file output
	std::ofstream fMagAbs("magAbsMetropolis.dat");

	for( double kappa = 0.10; kappa < 0.30; kappa += 0.005 )
	{
		Action act(lat, kappa, lambda);
		Metropolis updater( delta, act, phi, &rndGen );
		double acceptance = 0.;
		double av_mag = 0.;
		for( size_t u = 0; u<updates; u++ ) {
			for( size_t spu = 0; spu < stepsPerUpdate; spu++ )
			{
				double accepted = updater.update();
				acceptance +=accepted;

			}

			// measuring observables
			double magnetization = 0.;
			for( size_t x = 0; x < phi.getSize(); x++ )
				magnetization += phi(x);
			av_mag += std::fabs(magnetization)/phi.getSize();
		}
		double accRate = acceptance/double(updates*stepsPerUpdate);
		av_mag /= double(updates);
		std::cout << "Accepted "<< accRate << std::endl;
		std::cout << "Average magnetization: " << av_mag << std::endl;
		fMagAbs << kappa << "\t" << av_mag << std::endl;
	}
	fMagAbs.close();
}
