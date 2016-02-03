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

	FieldScalar<ScalarType>& operator+=( const FieldScalar<ScalarType>& rhs );
	FieldScalar<ScalarType>& operator-=( const FieldScalar<ScalarType>& rhs );
	FieldScalar<ScalarType>& operator*=( const ScalarType& rhs );
	FieldScalar<ScalarType>& operator/=( const ScalarType& rhs );

	template<class T> friend T operator*( FieldScalar<T> lhs, FieldScalar<T> rhs);

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

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator*=( const ScalarType& rhs ) {
	data *= rhs;
	return *this;
}

template<class ScalarType> FieldScalar<ScalarType>& FieldScalar<ScalarType>::operator/=( const ScalarType& rhs ) {
	data /= rhs;
	return *this;
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

template<class ScalarType> FieldScalar<ScalarType>& operator/( FieldScalar<ScalarType> lhs, const ScalarType& rhs ) {
	lhs /= rhs;
	return lhs;
}

template<class ScalarType> FieldScalar<ScalarType>& operator*( FieldScalar<ScalarType> lhs, const ScalarType& rhs ) {
	lhs *= rhs;
	return lhs;
}

/* Point-wise product of fields */
template<class ScalarType> ScalarType operator*( FieldScalar<ScalarType> lhs, FieldScalar<ScalarType> rhs) {
	return lhs.data.dot(rhs.data);
}

#endif /* FIELDSCALAR_H_ */
