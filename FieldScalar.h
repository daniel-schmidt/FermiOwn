/*
 * FieldScalar.h
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#ifndef FIELDSCALAR_H_
#define FIELDSCALAR_H_

#include <iostream>
#include <Eigen/Dense>
#include "Lattice.h"
enum InitType {
	zeroInit,
	oneInit,
	randomInit
};

typedef std::complex< double > Complex;
typedef double Real;

template<class ScalarType> class FieldScalar {
public:
	FieldScalar( const Lattice & newLat, InitType init );
	virtual ~FieldScalar();

	/**
	 * add a field to this field
	 * @param rhs is the field to add to the current one.
	 * @return the sum of the current vector and rhs
	 */
	FieldScalar<ScalarType>& operator+=( const FieldScalar<ScalarType>& rhs );

	/**
	 * @brief subtract a field from this field
	 * @param rhs is the field to subtract from the current one.
	 * @return the current vector minus rhs
	 */
	FieldScalar<ScalarType>& operator-=( const FieldScalar<ScalarType>& rhs );

	/**
	 * @brief Coefficient-wise multiply
	 * @param rhs is the field to multiply with
	 * @return the current vector with each element divided by the corresponding one in rhs
	 */
	FieldScalar<ScalarType>& operator*=( const FieldScalar<ScalarType>& rhs );

	/**
	 * @brief Coefficient-wise divide
	 * @param rhs is the field to divide by
	 * @return the current vector with each element divided by the corresponding one in rhs
	 */
	FieldScalar<ScalarType>& operator/=( const FieldScalar<ScalarType>& rhs );

	/**
	 * @brief Scalar multiplication
	 * @param rhs scalar to multiply with
	 * @return the current field multiplied by rhs
	 */
	FieldScalar<ScalarType>& operator*=( const ScalarType& rhs );

	/**
	 * @brief Scalar division
	 * @param rhs scalar to divide by
	 * @return the current field divided by rhs
	 */
	FieldScalar<ScalarType>& operator/=( const ScalarType& rhs );

	/**
	 * @brief Calculate scalar product
	 * @param rhs other vector to perform scalar product with
	 * @return scalar product of this field with rhs
	 */
	ScalarType dot( const FieldScalar<ScalarType>& rhs );

	void Print();

private:
	Eigen::Matrix< ScalarType, Eigen::Dynamic, 1> data;
	const Lattice & lat;
};

/* ---------------------------------------------------------------------------------------------------
 * Implementation of member functions
 * ---------------------------------------------------------------------------------------------------*/

template<class ScalarType> FieldScalar<ScalarType>::FieldScalar(const Lattice& newLat, InitType init) :
				lat(newLat)
				{
	data = Eigen::Matrix< ScalarType, Eigen::Dynamic, 1>(lat.getVol());
	switch( init ) {
	case zeroInit:
		data.setZero();
		return;
	case oneInit:
		data.setOnes();
		return;
	case randomInit:
		data.setRandom();
		return;
	}

				}

template<class ScalarType> FieldScalar<ScalarType>::~FieldScalar() {}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator+=( const FieldScalar<ScalarType>& rhs ) {
	data += rhs.data;
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator-=( const FieldScalar<ScalarType>& rhs ) {
	data -= rhs.data;
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator*=( const FieldScalar<ScalarType>& rhs ) {
	data = data.cwiseProduct(rhs.data);
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator/=( const FieldScalar<ScalarType>& rhs ) {
	data = data.cwiseQuotient(rhs.data);
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator*=( const ScalarType& rhs ) {
	data *= rhs;
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator/=( const ScalarType& rhs ) {
	data /= rhs;
	return *this;
}

template<class ScalarType> ScalarType FieldScalar<ScalarType>::dot( const FieldScalar<ScalarType>& rhs ) {
//	ScalarType ret = data.dot(rhs.data);
	return data.dot(rhs.data);
}

template<class ScalarType> void FieldScalar<ScalarType>::Print() {
	std::cout << data << std::endl;
}

/* ---------------------------------------------------------------------------------------------------
 * Implementation of non-member functions
 * ---------------------------------------------------------------------------------------------------*/

template<class ScalarType> inline FieldScalar<ScalarType> operator+( FieldScalar<ScalarType> lhs, const FieldScalar<ScalarType>& rhs ) {
	lhs += rhs;
	return lhs;
}

template<class ScalarType> inline FieldScalar<ScalarType> operator-( FieldScalar<ScalarType> lhs, const FieldScalar<ScalarType>& rhs ) {
	lhs -= rhs;
	return lhs;
}

template<class ScalarType> inline FieldScalar<ScalarType> operator*( FieldScalar<ScalarType> lhs, const FieldScalar<ScalarType>& rhs ) {
	lhs *= rhs;
	return lhs;
}

template<class ScalarType> inline FieldScalar<ScalarType> operator/( FieldScalar<ScalarType> lhs, const FieldScalar<ScalarType>& rhs ) {
	lhs /= rhs;
	return lhs;
}

template<class ScalarType> FieldScalar<ScalarType>& operator*( FieldScalar<ScalarType> lhs, const ScalarType& rhs ) {
	lhs *= rhs;
	return lhs;
}

template<class ScalarType> FieldScalar<ScalarType>& operator/( FieldScalar<ScalarType> lhs, const ScalarType& rhs ) {
	lhs /= rhs;
	return lhs;
}

/* Point-wise product of fields */
//template<class ScalarType> FieldScalar<ScalarType>& operator*( FieldScalar<ScalarType> lhs, const FieldScalar<ScalarType>& rhs ) {
//	lhs *= rhs;
//	return lhs.cwise;
//}
//template<class ScalarType> ScalarType operator*( FieldScalar<ScalarType> lhs, FieldScalar<ScalarType> rhs) {
//	return lhs.data.dot(rhs.data);
//}

#endif /* FIELDSCALAR_H_ */
