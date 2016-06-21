/*
 * WoodburyMatrix.cpp
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#include "WoodburyMatrix.h"

namespace FermiOwn {

WoodburyMatrix::WoodburyMatrix() :
	size( 0 ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true )
{};

WoodburyMatrix::WoodburyMatrix( const size_t matSize ) :
	mat( matSize, matSize ),
	size( matSize ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true )
{
}

WoodburyMatrix::WoodburyMatrix( const size_t matSize, const MatCoeffList & coeffList, const bool selfAdjointInit ) :
	mat( matSize, matSize ),
	size( matSize ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true )
{
	setFromCoeffList( coeffList, selfAdjointInit );
}

WoodburyMatrix::~WoodburyMatrix() {
}

void WoodburyMatrix::setFromCoeffList( const MatCoeffList & coeffList, const bool selfAdjointInit ) {
	mat.setFromTriplets( coeffList.begin(), coeffList.end() );

	if( selfAdjointInit ) {
		SparseMat tmp = mat.selfadjointView<Eigen::Upper>();
		mat = tmp;
	}

	matNeedsUpdate = false;

	Eigen::SparseLU< SparseMat > solver;	//TODO: check, if we can use an SPD solver here for SLAC
	solver.compute( mat );
	if( solver.info() != Eigen::Success ) {
		std::cerr << "WoodburyMatrix cannot calculate matrix decomposition to obtain determinant." << std::endl;
		exit(1);
	}

	det =  solver.determinant();
	detNeedsUpdate = false;

	SparseMat id( size, size );
	id.setIdentity();
	inv = solver.solve( id );

	if( solver.info() != Eigen::Success ) {
		std::cerr << "WoodburyMatrix failed to invert the matrix." << std::endl;
		exit(1);
	}
	invNeedsUpdate = false;

}

void WoodburyMatrix::setUpdateMatrices( const SparseMat & colMatrix, const SparseMat & rowMatrix ) {
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
	if( colMatrix.nonZeros() == 0 || rowMatrix.nonZeros() == 0 ) {
		std::cerr << "Error: WoodburyMatrix got update matrix with no non-zero element!" << std::endl;
		exit(1);
	}
	U = colMatrix;
	V = rowMatrix;

	detNeedsUpdate = true;
	invNeedsUpdate = true;
	matNeedsUpdate = true;
}

void WoodburyMatrix::update() {
	updateMatrix();
	updateDet();
	updateInverse();
}

void WoodburyMatrix::updateMatrix() {
	if( matNeedsUpdate ) {
		mat += U * V;	//TODO: setting entries directly is most likely more efficient
//		mat.prune( Complex( 0., 0.) ); could be used to delete entries from matrix that got zero but are still listed as elements, but does a copy of the matrix...
		matNeedsUpdate = false;
	}
}

Complex WoodburyMatrix::updateDet() {
	Complex detChange = 1.;
	if( detNeedsUpdate ) {
		VTimesInv = V*inv;
		smallId.resize( V.rows(), U.cols() );
		smallId.setIdentity();
		capacitanceMatrixLU.compute( VTimesInv * U + smallId );
		if( capacitanceMatrixLU.info() == Eigen::Success ) {
			detChange = capacitanceMatrixLU.determinant();
		} else {
			detChange = 0.;
		}
		det *= detChange;
		detNeedsUpdate = false;
	}
	return detChange;
}

void WoodburyMatrix::updateInverse() {
	if( invNeedsUpdate ) {
		if( detNeedsUpdate ) updateDet();			// we need the updated LU of the capacitance matrix, which is calculated in the updateDet().
		SparseMat capInv = capacitanceMatrixLU.solve( smallId );
		if( capacitanceMatrixLU.info() == Eigen::Success ) {
			inv -= inv * U * capInv * VTimesInv;
//		inv.prune( Complex( 0., 0.) ); could be used to delete entries from matrix that got zero but are still listed as elements, but does a copy of the matrix...
		} else {
			std::cerr << "Warning: WoodburyMatrix was not able to invert capacitance matrix needed to update inverse matrix!" << std::endl;
		}
		invNeedsUpdate = false;
	}
}

} /* namespace FermiOwn */