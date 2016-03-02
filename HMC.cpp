/*
 * HMC.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "HMC.h"

HMC::HMC() {
	// TODO Auto-generated constructor stub

}

HMC::~HMC() {
	// TODO Auto-generated destructor stub
}

void HMC::update() {

	// integration
	momentum.setGaussian();
	double H_old = Hamiltonian();
	FieldScalar<Real> old = phi;
	integrator.integrate();
	double H_new = Hamiltonian();

	double dH = std::exp(H_old-H_new);

	// accept/reject step
	bool accepted = false;
	double r = uniformDistribution01(randomGenerator);
	if( r < dH )
		accepted = true;
	else
		phi = old;
}

double HMC::Hamiltonian() {
	double H = act.getAction(phi);
	H += 0.5 * momentum.dot(momentum);
	return H;
}
