/*
 * Action.cpp
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#include "Action.h"

Action::Action(FieldScalar<Real>& field, double new_kappa, double new_lambda) :
	phi(field),
	lat( phi.getLattice() ),
	kappa(new_kappa),
	lambda(new_lambda)
{}

Action::~Action() {}

double Action::getAction() {
	// potential part of the action
	double S = phi.dot(phi) + lambda*(phi*phi-1.).dot(phi*phi-1.);

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

double Action::getForce() {
	//TODO implement force
	return 0;
}
