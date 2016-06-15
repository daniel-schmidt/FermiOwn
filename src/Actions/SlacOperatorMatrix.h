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
#include <unsupported/Eigen/KroneckerProduct>
#include "FieldBoolean.h"
#include "CliffordAlgebra.h"
#include "WoodburyMatrix.h"

namespace FermiOwn {

/**
 * @brief A class for a matrix-representation of the SLAC derivative
 *
 * This class provides methods to add and remove rows and columns from the full matrix
 * by a number of different methods.
 * It keeps track of changes and updates the determinant and the inverse of the matrix
 * using low-rank update formulae as the Woodbury lemma.
 */

class SlacOperatorMatrix {
public:

	enum updateType {
		eraseUpdate,
		separateUpdate,
		combinedUpdate
	};

	/**
	 * @brief Constructor for a 1D SLAC operator without flavour or spin structure
	 *
	 * @param size is the matrix size, equal to the number of lattice points
	 */
	SlacOperatorMatrix( size_t size );

	/**
	 * @brief Constructor for a general SLAC operator matrix with flavour and spin in any dimension
	 *
	 * Matrix size is Nt*Ns^(dim-1) * numFlavours * 2, where the factor of 2 is for the
	 * irreducible number of spin degrees of freedom.
	 *
	 * @param Nt is the temporal extend of the lattice, should be even for antiperiodic boundary conditions
	 * @param Ns is the spatial extend of the lattice, should be odd for periodic boundary conditions
	 * @param dim is the number of spacetime dimensions
	 * @param numFlavours is the number of flavours
	 */
	SlacOperatorMatrix( size_t Nt, size_t Ns, size_t dim, size_t numFlavours );

	virtual ~SlacOperatorMatrix();

	/**
	 * @brief Allowes direct read access to the SLAC operator matrix
	 * @return The matrix
	 */
	const Eigen::MatrixXcd getMatrix() const;

	/**
	 * @brief Allowes read access to the inverse of the current operator
	 * @return the inverse
	 */
	const Eigen::MatrixXcd getInverse() const;

	/**
	 * @brief Returns the determinant of the current matrix.
	 *
	 * @return The determinant.
	 */
	const Complex det() const;

	/**
	 * @brief returns linear matrix index for a given physical variables
	 *
	 * @param x is the spacetime point
	 * @param spin is the spinor component ( can be 0 or 1 )
	 * @param flavour is the number of flavour to access
	 * @return a linear index between 0 and the matrix size corresponding to the physical input values.
	 * 			Valid for addressing both rows and columns.
	 */
	size_t matIndex( size_t x, size_t spin, size_t flavour );

	/**
	 * @brief Method to add or remove column-row pairs from the matrix
	 *
	 * @param kxiab is boolean field where the current values indicate, which row/col pairs are to be deleted.
	 * @param changed is a boolean field indicating changes of the current kxiab to a previous state.
	 * @param upType is a selector for an update scheme.
	 */
	void update( FieldBoolean kxiab, FieldBoolean changed, updateType upType );

	/**
	 * @brief Reset the matrix to the initial state without any deleted rows or columns
	 */
	void setFull();

	/**
	 * @brief Deletes rows and columns corresponding to entries set in the input field.
	 *
	 * Translates the entries in the field to physical parameters and calls the single-
	 * point erase method. Afterwards, a full new calculation of determinant and inverse are performed.
	 *
	 * @param kxiab is a boolean field, where "true" indicates, that the corresponding row/column pair should be deleted from the matrix
	 */
	void erase( const FieldBoolean& kxiab );

	void resetState();
	void saveState();

private:

	/**
	 * @brief Returns a matrix for the one dimensional derivative via SLAC operator
	 *
	 * @param size is the number of lattice points in the direction, corresponding to the matrix size
	 * @return complex square matrix of size size representing the one dimensional derivative
	 */
	Eigen::MatrixXcd make1D( size_t size );

	/**
	 * @brief general implementation of the Woodbury formula for a low-rank update of determinant and inverse
	 *
	 * Calculates determinant and inverse of dslac + U.V
	 *
	 * @param U
	 * @param V
	 */
	void WoodburyUpdate( Eigen::MatrixXcd U, Eigen::MatrixXcd V );

	/**
	 * @brief Replace a single row and column by 0, with a 1 in the crossing entry
	 * The row and column to delete is determined by the physical parameters passed in.
	 *
	 * @param x is the spacetime point
	 * @param spin is the spinor component ( can be 0 or 1 )
	 * @param flavour1 is the number of flavour for the row to access
	 * @param flavour2 is the number of flavour for the column to access
	 */
	void erase( size_t x, size_t spin, size_t flavour1, size_t flavour2 );

	void combined( std::vector<size_t> addRows, std::vector<size_t> addCols, std::vector<size_t> delRows, std::vector<size_t> delCols);
	void deleteEntries( std::vector<size_t> rows, std::vector<size_t> cols );
	void addEntries( std::vector<size_t> rows, std::vector<size_t> cols );
	void update( size_t row1, size_t col1, size_t row2, size_t col2 );

	size_t N;						///< total number of physical points
	size_t dimSpinor;				///< number of spin degrees of freedom
	size_t Nf;						///< number of flavours
	Eigen::MatrixXcd dslac;			///< the current state of the operator matrix
	Eigen::MatrixXcd oldSlac;
	Eigen::MatrixXcd inverse;		///< the current state of the inverse of the operator, should multiply dslac to give identity
	Eigen::MatrixXcd oldInverse;
	Eigen::MatrixXcd fullSlac;		///< the initially constructed full matrix, used to restore elements or the full operator
	Eigen::MatrixXcd fullInverse;	//TODO: not needed?

	Complex detVal;					///< current value of the determinant
	Complex fullDet;				//TODO: not needed?
	Complex oldDet;

	CliffordAlgebra cliff;				///< The Clifford algebra used in the operator, necessary to get the gamma-matrices during construction of the operator.
	//TODO: needed only at construction, perhaps not necessary as field

	std::vector<size_t> deletedRows;	///< Vector keeping all indices of rows, that are deleted from the matrix
	std::vector<size_t> oldDelRows;
	std::vector<size_t> deletedCols;	///< Vector keeping all indices of columns, that are deleted from the matrix
	std::vector<size_t> oldDelCols;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
