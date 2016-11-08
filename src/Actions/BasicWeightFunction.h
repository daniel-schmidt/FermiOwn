/*
 * BasicWeightFunction.h
 *
 *  Created on: 04.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_BASICWEIGHTFUNCTION_H_
#define SRC_ACTIONS_BASICWEIGHTFUNCTION_H_

#include "BasicKField.h"
#include "DSlashUpdater.h"

namespace FermiOwn {

class BasicWeightFunction {
public:
	BasicWeightFunction( const BasicKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling );
	virtual ~BasicWeightFunction();

	/**
	 * @brief Calculates the full weight
	 * @return the value of the weight
	 */
	virtual Complex calculateWeight()=0;

	/**
	 * @brief Calculates the change in weight
	 * @return returns difference
	 */
	virtual Complex updateWeight( const std::set< size_t > & changedAt )=0;

	void saveState();

	inline void keep();
	inline void reset();

protected:
	const BasicKField & kfield;							///< The field to work with.
	BasicKField* initialField;						//TODO: this should be const somehow
	BasicKField* savedState;
	DSlashUpdater dslash;							///< The operator to work with.
	size_t V;										///< The lattice volume obtained from the field.
	double kappa;									///< The inverse coupling.
	bool needsKeepDecision;
};

/*=========================================
 * Inline functions
 *========================================= */

BasicWeightFunction::BasicWeightFunction( const BasicKField& boolField, size_t timeSize, size_t spatialSize, size_t dim, size_t numFlavours, double coupling ) :
	kfield( boolField ),
	dslash( timeSize, spatialSize, dim, numFlavours ),
	V( timeSize*pow( spatialSize, dim-1 ) ),
	kappa( 1./coupling ),
	needsKeepDecision( false )
{
	initialField = kfield.clone();
	savedState = kfield.clone();
}

BasicWeightFunction::~BasicWeightFunction() {
	delete initialField;
	delete savedState;
}

inline void BasicWeightFunction::saveState() {
	BasicKField* tmp = savedState;
	savedState = kfield.clone();
	delete tmp;
}

inline void BasicWeightFunction::keep() {
	if( !needsKeepDecision ) {
		std::cerr << "WeightFunction keep called, although not necessary." << std::endl;
		exit(1);
	}
	dslash.keep();
	needsKeepDecision = false;
}

inline void ThirringWeightFunction::reset(){
	if( !needsKeepDecision ) {
		std::cerr << "WeightFunction reset called, although not necessary." << std::endl;
		exit(1);
	}
	dslash.reset();
	needsKeepDecision = false;
}

} /* namespace FermiOwn */



#endif /* SRC_ACTIONS_BASICWEIGHTFUNCTION_H_ */
