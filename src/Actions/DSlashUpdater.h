/*
 * DSlashUpdater.h
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_DSLASHUPDATER_H_
#define SRC_ACTIONS_DSLASHUPDATER_H_

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

	void calculateUpdateMatrices( const FieldBoolean& kxiab, const FieldBoolean& change );
	Complex updateDet();
	void keep();
	void reset();

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

	WoodburyMatrix oldMat;
	WoodburyMatrix currMat;

	SparseMat fullOperator;

	std::vector<size_t> currentRows;
	std::vector<size_t> currentCols;

	bool updated;
};

inline size_t DSlashUpdater::matIndex( size_t x, size_t spin, size_t flavour ) const {
	return N*( dimSpinor*flavour + spin ) + x;
}

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_DSLASHUPDATER_H_ */
