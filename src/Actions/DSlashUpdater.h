/*
 * DSlashUpdater.h
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_DSLASHUPDATER_H_
#define SRC_ACTIONS_DSLASHUPDATER_H_

#include <set>
#include "CliffordAlgebra.h"
#include "FieldBoolean.h"
#include "WoodburyMatrix.h"

namespace FermiOwn {

/**
 * @brief This class represents the operator d-slash (gamma^mu partial_mu)
 * 	  	  and offers methods to delete rows and columns from it
 */
class DSlashUpdater {
public:
	DSlashUpdater();

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
	DSlashUpdater( size_t Nt, size_t Ns, size_t dim, size_t numFlavours );

	virtual ~DSlashUpdater();

	inline const SparseMat & getMatrix();
	inline const SparseMat & getInverse();
	inline Complex getDet();

	void calculateUpdateMatrices( const FieldBoolean& kxiab, const FieldBoolean& change );

	inline Complex updateDet();

	inline void keep();

	inline void reset();


	/**
	 * @brief returns linear matrix index for a given physical variables
	 *
	 * @param x is the spacetime point
	 * @param spin is the spinor component ( can be 0 or 1 )
	 * @param flavour is the number of flavour to access
	 * @return a linear index between 0 and the matrix size corresponding to the physical input values.
	 * 			Valid for addressing both rows and columns.
	 */
	inline size_t matIndex( size_t x, size_t spin, size_t flavour ) const;

	Eigen::MatrixXcd makeSlac1D( size_t size) const;

private:

	size_t N;						///< total number of physical points
	size_t dimSpinor;				///< number of spin degrees of freedom
	size_t Nf;						///< number of flavours
	size_t matSize;					///< the full operator is a (matSize x matSize)-matrix, usually matSize = N*dimSpinor*Nf

//	enum updateState {
//		initial,					///< operator is in the initial state, no update was prepared so far.
//		needsUpdate,				///< update matrices were set, but no updates executed
//		detUpdated,					///< only the determinant was updated with the currently set update matrices
////		fullUpdated					///< the full operator was updated, but the change was not accepted so far.
//	} status;						///< keep track of the internal state to prevent missuse

	WoodburyMatrix oldMat;
	WoodburyMatrix currMat;

	SparseMat fullOperator;

	std::vector<size_t> currentRows;
	std::vector<size_t> currentCols;
	std::vector<size_t> targetRows;
	std::vector<size_t> targetCols;

	bool changed;
};

inline const SparseMat & DSlashUpdater::getMatrix() {
	return currMat.getMatrix();
}

inline const SparseMat & DSlashUpdater::getInverse() {
	return currMat.getInverse();
}

inline Complex DSlashUpdater::getDet() {
	return currMat.getDet();
}

inline Complex DSlashUpdater::updateDet() {
//	if( changed == false ) std::cerr << "Warning: DSlacUpdater tries update of Determinant, but there are no changes. Perhaps you need to calculateUpdateMatrices first." << std::endl;
	return currMat.updateDet();
}

inline size_t DSlashUpdater::matIndex( size_t x, size_t spin, size_t flavour ) const {
	return N*( dimSpinor*flavour + spin ) + x;
}

inline void DSlashUpdater::keep() {
	if( std::fabs( getDet() ) < ZERO_TOL ) {
		std::cerr << "Error: you are trying to keep a DSlashUpdater state with determinant " << getDet() << std::endl;
		std::cerr << "This is not a valid state!" << std::endl;
		exit(1);
	}

	currMat.updateMatrix();
	currMat.updateInverse();
	currentRows = targetRows;	//TODO: can we save a copy here by exchanging pointers or something?
	currentCols = targetCols;
	changed = false;
}

inline void DSlashUpdater::reset() {
	currMat = oldMat;
	changed = false;
}

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_DSLASHUPDATER_H_ */
