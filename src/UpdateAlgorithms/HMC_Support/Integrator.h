/*
 * Integrator.h
 *
 *  Created on: 01.03.2016
 *      Author: dschmidt
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <random>

#include "Field.h"
#include "BasicAction.h"

namespace FermiOwn {

class Integrator {
public:
	Integrator( Field<Real>& position, Field<Real>& momentum, const BasicAction& act, const double new_t, const size_t new_nt );
	virtual ~Integrator();

	void integrate();

private:
	Field<Real>& x;
	Field<Real>& p;
	const BasicAction& act;	///< the action function to evaluate during integration
	const double t;		///< value of fictious time to integrate up to
	const size_t nt;	///< number of steps to take
	const double dt; 	///< step size resulting from t and nt
};

} // namespace FermiOwn

#endif /* INTEGRATOR_H_ */
