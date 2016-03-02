/*
 * Integrator.h
 *
 *  Created on: 01.03.2016
 *      Author: dschmidt
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <random>
#include "FieldScalar.h"
#include "Action.h"

class Integrator {
public:
	Integrator( FieldScalar<Real>& new_phi, const Action& act, std::ranlux48& rndGen, const double new_t, const size_t new_nt );
	virtual ~Integrator();

	void integrate();
	void invertMomentum();
private:
	FieldScalar<Real>& x;
	FieldScalar<Real> p;
	const Action& act;	///< the action function to evaluate during integration
	std::ranlux48& randomGenerator; ///< reference to random number generator engine
	const double t;		///< value of fictious time to integrate up to
	const size_t nt;	///< number of steps to take
	const double dt; 	///< step size resulting from t and nt
};

#endif /* INTEGRATOR_H_ */
