/*
 * Action_test.cpp
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#include <random>
#include "Lattice.h"
#include "Action.h"
#include "Field.h"

// TODO: rewrite as true test, include in makefile

bool Action_test(){
	using namespace FermiOwn;
	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	std::ranlux48 rndGen;

	Field<Real> phi(lat.getVol(), 1, &rndGen, oneInit);

	double lambda = 0.5;
	double kappa = 0.1;
	double Spot = kappa*phi.cwiseMultAndSum(phi) + lambda*(phi*phi).cwiseMultAndSum(phi*phi);
	std::cout << "Potential is: " << Spot  << " and should be 4.8" << std::endl;

	Action act( lat, kappa, lambda );
	std::cout << "Action is: " << act.getAction( phi )  << " and should be 3.2" << std::endl;
	std::cout << "Force is: ";
	act.getForce(phi).Print();
	std::cout<< " and should be vector of 14.8s." << std::endl;
	lambda = 0.9;
	kappa = 0.3;
	Action act2( lat, kappa, lambda );
	std::cout << "Action is: " << act2.getAction(phi)  << " and should be -6.4" << std::endl;
	return true;
}

int main() {
	int ret = 1;
	if( Action_test() ) {
		ret = 0;
	}
	return ret;
}
