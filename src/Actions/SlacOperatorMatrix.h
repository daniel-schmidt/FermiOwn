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

	size_t matIndex( size_t x, size_t spin, size_t flavour );

	void erase( const FieldBoolean& kxiab );
	void erase( size_t x, size_t spin, size_t flavour1, size_t flavour2 );

	void update( size_t row1, size_t col1, size_t row2, size_t col2 );
	void deleteEntries( std::vector<size_t> rows, std::vector<size_t> cols );
	void addEntries( std::vector<size_t> rows, std::vector<size_t> cols );
	void update( FieldBoolean kxiab, FieldBoolean changed );

	void setFull();

	Eigen::MatrixXcd inverse;
private:

	Eigen::MatrixXcd make1D( size_t size );

	size_t N;					///< total number of physical points
	size_t dimSpinor;			///< number of spin degrees of freedom
	size_t Nf;					///< number of flavours
	Eigen::MatrixXcd dslac;
	Eigen::MatrixXcd fullSlac;
	Eigen::MatrixXcd fullInverse;

	Complex detVal;
	Complex fullDet;

	CliffordAlgebra cliff;

	std::vector<size_t> deletedRows;
	std::vector<size_t> deletedCols;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
