/*
 * Action.cpp
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#include "Action.h"

Action::Action( const double new_kappa, const double new_lambda) :
	kappa(new_kappa),
	lambda(new_lambda)
{}

Action::~Action() {}

double Action::getAction( const FieldScalar<Real>& phi) const {
	// potential part of the action
	double S = phi.dot(phi) + lambda*(phi*phi-1.).dot(phi*phi-1.);

	for( size_t x = 0; x < phi.getLattice().getVol(); x++ ) {
		double kinetic = 0;
		std::vector<size_t> nnIndex = phi.getLattice().getNeighbours(x,Lattice::fwd);
		for( size_t nn:nnIndex ) {
			kinetic += phi(nn);
		}
		S -= 2*kappa*phi(x) * kinetic;
	}

	return S;
}

FieldScalar<Real> Action::getForce( const FieldScalar<Real>& phi) const {
	FieldScalar<Real> force = ( 2. + 4.*lambda*(phi.dot(phi)-1.) )*phi;

	for( size_t x = 0; x < phi.getLattice().getVol(); x++ ) {
		double kinetic = 0;
		std::vector<size_t> nnIndex = phi.getLattice().getNeighbours(x);
		for( size_t nn:nnIndex ) {
			kinetic += phi(nn);
		}
		force(x) -= 2 * kappa * kinetic;
	}

	return force;
}
