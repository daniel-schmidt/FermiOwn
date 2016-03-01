/*
 * Action.h
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#ifndef ACTION_H_
#define ACTION_H_

#include "FieldScalar.h"
#include "Lattice.h"

class Action {
public:
	Action( FieldScalar<Real>& field, double new_kappa, double new_lambda);
	virtual ~Action();

	double getAction();
	double getForce();

private:
	FieldScalar<Real>& phi;
	const Lattice& lat;
	double kappa;
	double lambda;
};

#endif /* ACTION_H_ */
