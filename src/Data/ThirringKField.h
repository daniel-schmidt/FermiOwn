/*
 * ThirringKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_THIRRINGKFIELD_H_
#define SRC_DATA_THIRRINGKFIELD_H_

#include "BasicCloneKField.h"

namespace FermiOwn {

class ThirringKField: public BasicCloneKField<ThirringKField> {
public:
	ThirringKField( const size_t latticeVolume, const std::vector<size_t>& internal );
	virtual ~ThirringKField();

	void enforceConstraint( size_t x, const std::vector<size_t>& internal );
	bool constraintViolated( size_t x) const;
	/**
	 * @brief Performs coefficient-wise compare with another field
	 *
	 * @param other is a second matrix to compare with
	 * @return boolean matrix, whose entries are "true", if this field differs from the other, and otherwise "false".
	 */
	ThirringKField different( ThirringKField other ) const;

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

	/**
	 * @brief Counts the number of offdiagonal flavours kxab set to 2
	 * @return number of offdiagonal flavour entries equal to 2
	 */
	size_t countOffdiagonal2() const;

protected:
	virtual inline size_t colIndex( const std::vector<size_t>& internal ) const;
};

} /* namespace FermiOwn */

#endif /* SRC_DATA_THIRRINGKFIELD_H_ */
