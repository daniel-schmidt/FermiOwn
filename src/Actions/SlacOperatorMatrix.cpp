/*
 * SlacOperatorMatrix.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include "SlacOperatorMatrix.h"

namespace FermiOwn {

SlacOperatorMatrix::SlacOperatorMatrix( size_t size ) :
		N(size),
		dslac(N,N)
{
	for( size_t i = 0; i < N; i++ ) {
		for( size_t j = 0; j < N; j++ ) {
			if( i == j ) {
				dslac(i,j) = 0.;
			} else {
				int t = i-j;
				dslac(i,j) = PI*pow(-1, t) / ( N * sin(PI*t/N) );
			}
		}
	}

}

SlacOperatorMatrix::~SlacOperatorMatrix() {
	// TODO Auto-generated destructor stub
}

const Eigen::MatrixXd SlacOperatorMatrix::getMatrix() const {
	return dslac;
}

} /* namespace FermiOwn */


