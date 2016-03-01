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
	/**
	 * @brief calculates the force for the scalar field
	 * @return a field with the force at every lattice point
	 */
	FieldScalar<Real> getForce();

private:
	FieldScalar<Real>& phi;
	const Lattice& lat;
	double kappa;
	double lambda;
};

#endif /* ACTION_H_ */
