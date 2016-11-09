/*
 * WeightFunctionGrossNeveu.h
 *
 *  Created on: 09.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_WEIGHTFUNCTIONGROSSNEVEU_H_
#define SRC_ACTIONS_WEIGHTFUNCTIONGROSSNEVEU_H_

#include <cmath>
#include "BasicWeightFunctionTemplate.h"
#include "GrossNeveuKField.h"
#include "MatrixChanges.h"

namespace FermiOwn {

class WeightFunctionGrossNeveu: public BasicWeightFunctionTemplate<GrossNeveuKField> {
public:
	WeightFunctionGrossNeveu( const GrossNeveuKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling );
	virtual ~WeightFunctionGrossNeveu();

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
	MatrixChanges<GrossNeveuKField> change;
	std::vector<Complex> locWeights;	//TODO: should be possible to use double weights together with hermitian slac operator...
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_WEIGHTFUNCTIONGROSSNEVEU_H_ */
