/*
 * Integrator.cpp
 *
 *  Created on: 01.03.2016
 *      Author: dschmidt
 */

#include "Integrator.h"

namespace FermiOwn {

Integrator::Integrator( Field<Real>& position, Field<Real>& momentum, const BasicAction& new_act, const double new_t, const size_t new_nt ) :
	x(position),
	p(momentum),
	act(new_act),
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
	for( size_t n = 0; n < nt-1; n++ ) {
		x += p * dt;
		p -= dt*act.getForce(x);
	}
	x += p * dt;
	p -= dt/2.*act.getForce(x);
}

} // namespace FermiOwn
