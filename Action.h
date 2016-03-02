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
	Action( const double new_kappa, const double new_lambda);
	virtual ~Action();

	double getAction( const FieldScalar<Real>& phi ) const;

	/**
	 * @brief calculates the force for the scalar field
	 * More exactly, it calculates the gradient of the action, which is more likely the negative force...
	 * @return a field with the force at every lattice point
	 */
	FieldScalar<Real> getForce( const FieldScalar<Real>& phi ) const;

private:
	const double kappa;
	const double lambda;
};

#endif /* ACTION_H_ */
