/*
 * FieldBoolean.h
 *
 *  Created on: 30.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_FIELDBOOLEAN_H_
#define SRC_DATA_FIELDBOOLEAN_H_

#include "Field.h"

namespace FermiOwn {

class FieldBoolean: public Field<bool> {
public:
	FieldBoolean( const size_t latticeVolume, const size_t numberOfSpin, const size_t numberOfFlavours, std::ranlux48 * rndGen, InitType init );
	virtual ~FieldBoolean();
	bool getValue( size_t x, size_t spin, size_t flavour1, size_t flavour2 ) const;
	void setValue( bool val, size_t x, size_t spin, size_t flavour1, size_t flavour2 );
	void invert(size_t x, size_t spin, size_t flavour1, size_t flavour2);
	bool constraintViolated( size_t x) const;

	/**
	 * @brief Sum over all components
	 *
	 * This is the number of entries set true, since the field values can be only 0 or 1.
	 * @return sum over all components: spacetime, spin, flavour (diagonal and offdiagonal)
	 */
	size_t sumAll() const;

	/**
	 * @brief Sums the field values for each flavour diagonal component over all spins and counts how often each possible value occurs
	 *
	 * @param x is the fixed spacetime point where to sum and count
	 * @return a vector with bins for 0, ... , Nf each containing the number of kxaa with this value
	 */
	Eigen::ArrayXi countSummedSpin( size_t x ) const;

private:
	size_t colIndex( size_t spin, size_t flavour1, size_t flavour2 ) const;

	size_t numSpins;
	size_t numFlavours;
	size_t numColsPerSpin;	///< number of columns for each spin degree: Nf*(Nf+1)/2, that is Nf for kxiaa and Nf*(Nf-1)/2 for kxiab with a!=b
};

} /* namespace FermiOwn */

#endif /* SRC_DATA_FIELDBOOLEAN_H_ */
