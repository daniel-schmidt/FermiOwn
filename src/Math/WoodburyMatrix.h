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

/**
 * @brief A sparse square matrix, that can be updated using the Woodbury formula
 *
 * This class provides storage for a general matrix, with methods to update this matrix by low-rank update matrices.
 * This allows to keep the value of the inverse and the determinant updated with low computational cost and without
 * explicitly re-calculating it for the updated matrix.
 */
class WoodburyMatrix {
public:

	/**
	 * @brief Default constructor, constructing empty 0 x 0 matrix.
	 */
	WoodburyMatrix();

	/**
	 * @brief Constructs an uninitialized matrix of size (matSize x matSize).
	 * @param matSize is the size of the matrix: we construct a (matSize x matSize) matrix.
	 */
	WoodburyMatrix( const size_t matSize );

	/**
	 * @brief Constructs a new matrix with given size and coefficients.
	 *
	 * This calls setFromCoeffList() which gives more instructions.
	 *
	 * @param matSize is the size of the matrix: we construct a (matSize x matSize) matrix.
	 * @param list is a list of coefficients for the sparse matrix. The entries must have the form ( i, j, m_ij ).
	 * @param selfAdjointInit should be true, if the matrix is selfadjoint and the coefficients are given for the upper triangular part only.
	 */
	WoodburyMatrix( const size_t matSize, const MatCoeffList & list, const bool selfAdjointInit = false );


//	WoodburyMatrix( const WoodburyMatrix& other );

	/**
	 * @brief Initialize matrix with a list of coefficients.
	 *
	 * Initializes the internal sparse matrix from the coefficient list. If the matrix to be constructed is selfadjoint,
	 * only half of the non-zero coefficients are needed for initialization. They are expected to be in the upper
	 * triangular part of the matrix. The remaining coefficients are obtained by conjugate-transpose and copy to the
	 * final matrix storage.
	 *
	 * This also calculates the inverse and the determinant of the matrix, which may be time consuming.
	 *
	 * @param list is a list of coefficients for the sparse matrix. The entries must have the form ( i, j, m_ij ).
	 * @param selfAdjointInit should be true, if the matrix is selfadjoint and the coefficients are given for the upper triangular part only.
	 */
	void setFromCoeffList( const MatCoeffList & list, const bool selfAdjointInit = false );

	virtual ~WoodburyMatrix();

//	inline WoodburyMatrix& operator=( WoodburyMatrix other );

//	friend void swap( WoodburyMatrix& first, WoodburyMatrix& second );

	/**
	 * @brief Returns the determinant of the matrix.
	 *
	 * Performes calculations, if new update matrices were specified, but no update was executed.
	 *
	 * @return the determinant of the matrix.
	 */
	inline Complex getDet();

	/**
	 * @brief Returns the matrix itself.
	 *
	 * If needed, an update of the matrix with the current update matrices is performed.
	 *
	 * @return const reference to the matrix in an Eigen sparse matrix format.
	 */
	inline const SparseMat& getMatrix();

	/**
	 * @brief Returns the inverse of the matrix.
	 *
	 * If needed, an update of the inverse with the current update matrices is performed.
	 *
	 * @return const reference to the inverse of the matrix in an Eigen sparse matrix format.
	 */
	inline const SparseMat& getInverse();

	/**
	 * @brief Prepares update by passing update matrices to the class.
	 *
	 * No calculations are performed. The update is given by this + colMat*rowMat
	 *
	 * @param colMatrix is a matrix of dimensions ( size x m ), thus specifying m columns to update.
	 * @param rowMatrix is a matrix of dimensions ( m x size ), thus specifying m rows to update.
	 */
	void setUpdateMatrices( const SparseMat* colMatrix, const SparseMat* rowMatrix );

	void setUpdateMatricesDeleteOnly( const SparseMat* colMatrix, const SparseMat* rowMatrix, const std::vector<size_t>& colsToDelete, const std::vector<size_t>& rowsToDelete );

	/**
	 * @brief Update the matrix, the inverse and the determinant.
	 */
	void update();

	/**
	 * @brief Update only the matrix
	 * Determinant and inverse are not touched and may not be up to date afterwards.
	 */
	void updateMatrix();

	/**
	 * @brief Update only the determinant
	 *
	 * This may be useful to check, if the update leads to an invertible matrix,
	 * such that we do not need to proceed, if this is not the case
	 *
	 * @return The change in the determinant, which is a factor to multiply the old determinant with, to get the current one.
	 */
	Complex updateDet();

	/**
	 * @brief Update the inverse, which always needs the determinant to be updated.
	 *
	 * Thus, if the determinant is not updated jet, we call updateDet() to catch up on this.
	 */
	void updateInverse();

	/**
	 * @brief Print the matrix to console.
	 *
	 * If necessary, this performs updates to display the correct state.
	 */
	inline void Print();
	/**
	 * @brief Print the inverse to console.
	 *
	 * If necessary, this performs updates to display the correct state.
	 */
	inline void PrintInverse();

private:
	SparseMat mat;					///< the actual matrix storage
	size_t size;					///< the matrix is of dimension (size x size)

	SparseMat inv;					///< storage for the inverse matrix
	Complex det;					///< the determinant of the matrix

	bool matNeedsUpdate;			///< switch to check, if we have to update the matrix
	bool detNeedsUpdate;			///< switch to check, if we have to update the determinant
	bool invNeedsUpdate;			///< switch to check, if we have to update the inverse

	const SparseMat* U;					///< the left update vector for A + U*V
	const SparseMat* V;					///< the right update vector for A + U*V

	Eigen::SparseLU< SparseMat > capacitanceMatrixLU;	///< a full-pivoted LU decomposition of the matrix 1 + V*inv*U (the capacitance matrix) used to update determinant and inverse
	SparseMat VTimesInv;			///< temporary storage for the product V*inv, to prevent double evaluation if updating determinant and inverse
	SparseMat smallId;				///< identitiy matrix in the small subspace, needs to be stored for separate evaluation of determinant and inverse

	bool deleteOnly;				///< switch on to use faster updates, if only rows/cols are set to zero.
	std::vector<size_t> rows;
	std::vector<size_t> cols;
	Eigen::MatrixXcd submat;
};

/*============================================================
 * Inline Functions
 *============================================================*/

//inline WoodburyMatrix& WoodburyMatrix::operator=( WoodburyMatrix other ) {
//	swap( *this, other );
//	return *this;
//}

inline Complex WoodburyMatrix::getDet() {
	if( detNeedsUpdate ) updateDet();
	return det;
}

inline const SparseMat& WoodburyMatrix::getMatrix() {
	if( matNeedsUpdate ) updateMatrix();
	return mat;
}

inline const SparseMat& WoodburyMatrix::getInverse() {
	if( invNeedsUpdate ) updateInverse();
	return inv;
}

inline void WoodburyMatrix::Print() {
	if( matNeedsUpdate ) updateMatrix();
	std::cout << mat << std::endl;
}

inline void WoodburyMatrix::PrintInverse() {
	if( invNeedsUpdate ) updateInverse();
	std::cout << inv << std::endl;
}


} /* namespace FermiOwn */

#endif /* SRC_MATH_WOODBURYMATRIX_H_ */
