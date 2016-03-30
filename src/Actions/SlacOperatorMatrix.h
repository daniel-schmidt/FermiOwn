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
#include "FieldBoolean.h"
#include "CliffordAlgebra.h"

namespace FermiOwn {

class SlacOperatorMatrix {
public:
	SlacOperatorMatrix( size_t size );
	SlacOperatorMatrix( size_t Nt, size_t Ns, size_t dim, size_t numFlavours );
	virtual ~SlacOperatorMatrix();

	const Eigen::MatrixXcd getMatrix() const;
	const Complex det() const;
//	void deletePoint( size_t x );
	void erase( const FieldBoolean& kxiab );
	void erase( size_t x, size_t spin, size_t flavour1, size_t flavour2 );

//	void eraseCol( size_t x, size_t spin );
//	void eraseRow( size_t x, size_t spin );
//	void eraseCols( VectorXb kx0, VectorXb kx1 );
//	void eraseRows( VectorXb kx0, VectorXb kx1 );
	void setFull();
//	void addPoint( size_t x );
private:

	void resetMaps();
	Eigen::MatrixXcd make1D( size_t size );

	size_t N;					///< total number of physical points
	size_t dimSpinor;			///< number of spin degrees of freedom
	size_t Nf;					///< number of flavours
	Eigen::MatrixXcd dslac;
	Eigen::MatrixXcd fullSlac;
	CliffordAlgebra cliff;
	Eigen::VectorXi rowMap;	///< mapping between indices of full matrix and non-deleted ones, if row is deleted, this is -1, otherwise the row number in not deleted submatrix
	Eigen::VectorXi colMap;	///< mapping between indices of full matrix and non-deleted ones, if col is deleted, this is -1, otherwise the col number in not deleted submatrix
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
