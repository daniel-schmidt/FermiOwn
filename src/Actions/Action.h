/*
 * Action.h
 *
 *  Created on: 05.02.2016
 *      Author: dschmidt
 */

#ifndef ACTION_H_
#define ACTION_H_

#include "Field.h"
#include "Lattice.h"
#include "BasicAction.h"

namespace FermiOwn {

class Action: public BasicAction {
public:
	Action( const Lattice& new_lat, const double new_kappa, const double new_lambda);
	virtual ~Action();

	double getAction( const Field<Real>& phi ) const;

	Field<Real> getForce( const Field<Real>& phi ) const;

private:
	const Lattice& lat;
	const double kappa;
	const double lambda;
};

} // namespace FermiOwn

#endif /* ACTION_H_ */
