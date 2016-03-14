/*
 * SlacOperatorMatrix.h
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_SLACOPERATORMATRIX_H_
#define SRC_ACTIONS_SLACOPERATORMATRIX_H_

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include "Constants.h"
#include "CliffordAlgebra.h"

namespace FermiOwn {

class SlacOperatorMatrix {
public:
	SlacOperatorMatrix( size_t size );
	SlacOperatorMatrix( size_t Nt, size_t Ns, size_t dim );
	virtual ~SlacOperatorMatrix();

	const Eigen::MatrixXcd getMatrix() const;

private:

	Eigen::MatrixXcd make1D( size_t size );

	size_t N;
	Eigen::MatrixXcd dslac;
	CliffordAlgebra cliff;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
