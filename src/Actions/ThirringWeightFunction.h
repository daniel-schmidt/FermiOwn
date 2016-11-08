/*
 * WeightFunction.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_THIRRINGWEIGHTFUNCTION_H_
#define SRC_ACTIONS_THIRRINGWEIGHTFUNCTION_H_

#include <set>

#include "BasicWeightFunctionTemplate.h"
#include "ThirringKField.h"
#include "MatrixChanges.h"

namespace FermiOwn {

class ThirringWeightFunction : public BasicWeightFunctionTemplate<ThirringKField> {
public:
	ThirringWeightFunction( const ThirringKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling );
	virtual ~ThirringWeightFunction();

	/**
	 * @brief Calculates the full weight
	 * @return the value of the weight
	 */
	Complex calculateWeight();

	/**
	 * @brief Calculates the change in weight
	 * @return returns difference
	 */
	Complex updateWeight( const std::set< size_t > & changedAt );

private:
	double getHypergeometricFactor( int flavour );	//TODO: merge this with the two-argument function
	double getHypergeometricFactor( int n1, int n2 );

	MatrixChanges change;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_THIRRINGWEIGHTFUNCTION_H_ */
