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

typedef std::function< size_t( size_t, size_t, size_t ) > matrixIndexFun;

class ThirringKField: public BasicCloneKField<ThirringKField> {
public:
	ThirringKField( const size_t latticeVolume, const idxVec& internal );
	virtual ~ThirringKField();

	void enforceConstraint( size_t x, const idxVec& internal );
	bool constraintViolated( size_t x) const;

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
