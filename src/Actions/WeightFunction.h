/*
 * WeightFunction.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_WEIGHTFUNCTION_H_
#define SRC_ACTIONS_WEIGHTFUNCTION_H_

#include <set>
#include "FieldBoolean.h"
#include "DSlashUpdater.h"

namespace FermiOwn {

class WeightFunction {
public:
	WeightFunction( const FieldBoolean& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling );
	virtual ~WeightFunction();


	/**
	 * @brief Calculates the full weight
	 * @return the value of the weight
	 */
	Complex calculateWeight();

	void saveState();

	/**
	 * @brief Calculates the change in weight
	 * @return returns difference
	 */
	Complex updateWeight( const std::set< size_t > & changedAt );

	inline void keep();
	inline void reset();

private:

	double getHypergeometricFactor( int flavour );	//TODO: merge this with the two-argument function
	double getHypergeometricFactor( int n1, int n2 );

	const FieldBoolean & kxiab;							///< The field to work with.
	FieldBoolean initialField;
	FieldBoolean savedState;
	DSlashUpdater dslash;							///< The operator to work with.
	size_t V;										///< The lattice volume obtained from the field.
	double kappa;									///< The inverse coupling.
	bool needsKeepDecision;
};

/*=========================================
 * Inline functions
 *========================================= */

inline void WeightFunction::saveState() {
	savedState = kxiab;
}

inline void WeightFunction::keep() {
	if( !needsKeepDecision ) {
		std::cerr << "WeightFunction keep called, although not necessary." << std::endl;
		exit(1);
	}
	dslash.keep();
	needsKeepDecision = false;
}

inline void WeightFunction::reset(){
	if( !needsKeepDecision ) {
		std::cerr << "WeightFunction reset called, although not necessary." << std::endl;
		exit(1);
	}
	dslash.reset();
	needsKeepDecision = false;
}


} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_WEIGHTFUNCTION_H_ */
