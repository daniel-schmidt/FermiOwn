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
#include <assert.h>
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
	const Complex det() const;
	void deletePoint( size_t x );
	void eraseCols( VectorXb kx0, VectorXb kx1 );
	void eraseRows( VectorXb kx0, VectorXb kx1 );
	void setFull();
//	void addPoint( size_t x );
private:

	Eigen::MatrixXcd make1D( size_t size );

	size_t N;					///< total number of physical points
	size_t dimSpinor;
	Eigen::MatrixXcd dslac;
	Eigen::MatrixXcd fullSlac;
	CliffordAlgebra cliff;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
