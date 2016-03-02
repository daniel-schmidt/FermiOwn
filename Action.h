/*
 * Action.h
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#ifndef ACTION_H_
#define ACTION_H_

#include "BasicAction.h"
#include "FieldScalar.h"
#include "Lattice.h"

class Action: public BasicAction {
public:
	Action( const double new_kappa, const double new_lambda);
	virtual ~Action();

	double getAction( const FieldScalar<Real>& phi ) const;

	FieldScalar<Real> getForce( const FieldScalar<Real>& phi ) const;

private:
	const double kappa;
	const double lambda;
};

#endif /* ACTION_H_ */
