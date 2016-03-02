/*
 * Integrator.cpp
 *
 *  Created on: 01.03.2016
 *      Author: dschmidt
 */

#include "Integrator.h"

Integrator::Integrator( FieldScalar<Real>& new_phi, const Action& new_act, std::ranlux48& rndGen, const double new_t, const size_t new_nt ) :
	x(new_phi),
	p( FieldScalar<Real>(x.getLattice(), rndGen, randomInit) ),	//TODO: is randomInit the correct init type? uniform or gaussian?
	act(new_act),
	randomGenerator(rndGen),
	t(new_t),
	nt(new_nt),
	dt(t/nt)
{
}

Integrator::~Integrator() {
	// TODO Auto-generated destructor stub
}

void Integrator::integrate() {
	p -= dt/2.*act.getForce(x);	//TODO: test if the sign is correct
	for( size_t n = 0; n < t-1; n++ ) {
		x += p * dt;
		p -= dt*act.getForce(x);
	}
	x += p * dt;
	p -= dt/2.*act.getForce(x);
}
