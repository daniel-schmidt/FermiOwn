/*
 * BasicAction.h
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#ifndef BASICACTION_H_
#define BASICACTION_H_

#include "FieldScalar.h"

class BasicAction {
public:
	BasicAction(){};
	virtual ~BasicAction(){};

	virtual double getAction( const FieldScalar<Real>& phi ) const =0;

	/**
	 * @brief calculates the force for the scalar field
	 * More exactly, it calculates the gradient of the action, which is more likely the negative force...
	 * @return a field with the force at every lattice point
	 */
	virtual FieldScalar<Real> getForce( const FieldScalar<Real>& phi ) const =0;

};

#endif /* BASICACTION_H_ */