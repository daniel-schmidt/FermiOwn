/*
 * Action.cpp
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#include "Action.h"

namespace FermiOwn {

Action::Action( const Lattice& new_lat, const double new_kappa, const double new_lambda) :
	lat(new_lat),
	kappa(new_kappa),
	lambda(new_lambda)
{}

Action::~Action() {}

double Action::getAction( const Field<Real>& phi) const {
	// potential part of the action
	double S = phi.cwiseMultAndSum(phi) + lambda*(phi*phi-1.).cwiseMultAndSum(phi*phi-1.);

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		double kinetic = 0;
		std::vector<size_t> nnIndex = lat.getNeighbours(x,Lattice::fwd);
		for( size_t nn:nnIndex ) {
			kinetic += phi(nn);
		}
		S -= 2*kappa*phi(x) * kinetic;
	}

	return S;
}

Field<Real> Action::getForce( const Field<Real>& phi) const {
	Field<Real> force = ( 2. + 4.*lambda*(phi.cwiseMultAndSum(phi)-1.) )*phi;

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		double kinetic = 0;
		std::vector<size_t> nnIndex = lat.getNeighbours(x);
		for( size_t nn:nnIndex ) {
			kinetic += phi(nn);
		}
		force(x) -= 2 * kappa * kinetic;
	}

	return force;
}

} // namespace FermiOwn
