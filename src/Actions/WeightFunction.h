/*
 * WeightFunction.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_WEIGHTFUNCTION_H_
#define SRC_ACTIONS_WEIGHTFUNCTION_H_

#include "FieldBoolean.h"
#include "SlacOperatorMatrix.h"

namespace FermiOwn {

class WeightFunction {
public:
	WeightFunction( FieldBoolean& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling );
	virtual ~WeightFunction();


	/**
	 * @brief Calculates the full weight
	 * @return the value of the weight
	 */
	Complex calculateWeight();

	/**
	 * @brief Calculates the change in weight
	 * @return returns difference
	 */
	Complex updateWeight();

private:

	double getHypergeometricFactor( int flavour );	//TODO: merge this with the two-argument function
	double getHypergeometricFactor( int n1, int n2 );

	FieldBoolean & kxiab;							///< The field to work with.
	SlacOperatorMatrix slac;						///< The operator to work with.
	size_t V;										///< The lattice volume obtained from the field.
	double kappa;									///< The inverse coupling.
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_WEIGHTFUNCTION_H_ */
