/*
 * Action.h
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#ifndef ACTION_H_
#define ACTION_H_

#include "FieldScalar.h"

class Action {
public:
	Action( FieldScalar<Real>& field, double mass, double coupling );
	virtual ~Action();

	double getAction();
	double getForce();

private:
	FieldScalar<Real>& phi;
	double m;
	double lambda;
};

#endif /* ACTION_H_ */
