/*
 * FieldScalar.cpp
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#include "FieldScalar.h"

// template specialization for real numbers
template<> void FieldScalar<Real>::setGaussian() {
	auto gaussian = [&] (Real) {return normalDistribution01(randomGenerator); };
	data = Eigen::Matrix< Real, Eigen::Dynamic, 1>::NullaryExpr( data.size(), gaussian );
}
// template specialization for complex numbers: draw two reals
template<> void FieldScalar<Complex>::setGaussian() {
	auto gaussian = [&] (Real) {return normalDistribution01(randomGenerator); };
	data.real() = Eigen::Matrix< Real, Eigen::Dynamic, 1>::NullaryExpr( data.size(), gaussian );
	data.imag() = Eigen::Matrix< Real, Eigen::Dynamic, 1>::NullaryExpr( data.size(), gaussian );
}


