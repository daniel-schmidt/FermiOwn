/*
 * WoodburyMatrix.h
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_MATH_WOODBURYMATRIX_H_
#define SRC_MATH_WOODBURYMATRIX_H_

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "Constants.h"

namespace FermiOwn {

typedef std::vector< Eigen::Triplet<Complex> > MatCoeffList;
typedef Eigen::SparseMatrix<Complex, Eigen::ColMajor> SparseMat;
class WoodburyMatrix {
public:
	WoodburyMatrix( const size_t matSize, const MatCoeffList & list );
	virtual ~WoodburyMatrix();

	inline Complex getDet() const;

	inline void setUpdateMatrices( SparseMat colMatrix, SparseMat rowMatrix );

	void updateMatrix();
	Complex updateDet();
	void updateInverse();

	inline void Print() const;
	inline void PrintInverse() const;



private:
	SparseMat mat;
	const size_t size;

	SparseMat inv;
	Complex det;

	bool validState;

	SparseMat U;
	SparseMat V;

//	Eigen::FullPivLU< Eigen::MatrixXcd > capacitanceMatrixLU;	///< a full-pivoted LU decomposition of the matrix 1 + V*inv*U used to update determinant and inverse
	Eigen::SparseLU< SparseMat > capacitanceMatrixLU;	///< a full-pivoted LU decomposition of the matrix 1 + V*inv*U used to update determinant and inverse
	SparseMat VTimesInv;			///< temporary storage for the product V*inv, to prevent double evaluation if updating determinant and inverse
	SparseMat smallId;				///< identitiy matrix in the small subspace
};

/*============================================================
 * Inline Functions
 *============================================================*/

inline Complex WoodburyMatrix::getDet() const {
	return det;
}

inline void WoodburyMatrix::setUpdateMatrices( SparseMat colMatrix, SparseMat rowMatrix ) {
	if( colMatrix.rows() != int( size ) ) {
		std::cerr << "WoodburyMatrix got a colMatrix of wrong size, not matching the matrix dimension " << size << ", instead size " << colMatrix.rows() << std::endl;
		exit(1);
	}
	if( rowMatrix.cols() != int( size ) ) {
		std::cerr << "WoodburyMatrix got a rowMatrix of wrong size, not matching the matrix dimension " << size << ", instead size " << rowMatrix.cols() << std::endl;
		exit(1);
	}
	if( colMatrix.cols() != rowMatrix.rows() ) {
		std::cerr << "The two update matrices passed to WoodburyMatrix do not fit in size! They are col: " << colMatrix.cols() << " and row: " << rowMatrix.rows() << std::endl;
		exit(1);
	}
	U = colMatrix;
	V = rowMatrix;

	validState = false;
}

inline void WoodburyMatrix::Print() const {
	std::cout << mat << std::endl;
}

inline void WoodburyMatrix::PrintInverse() const {
	std::cout << inv << std::endl;
}


} /* namespace FermiOwn */

#endif /* SRC_MATH_WOODBURYMATRIX_H_ */
