/*
 * Action.cpp
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#include "Action.h"

Action::Action(FieldScalar<Real>& field, double mass, double coupling) :
	phi(field),
	m(mass),
	lambda(coupling)
{}

Action::~Action() {}

double Action::getAction() {
	//TODO kinetic Term
	double S = m*phi.dot(phi) + lambda*(phi*phi).dot(phi*phi);
	return S;
}

double Action::getForce() {
	//TODO implement force
	return 0;
}
