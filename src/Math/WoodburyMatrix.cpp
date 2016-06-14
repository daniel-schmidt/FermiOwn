/*
 * WoodburyMatrix.cpp
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#include "WoodburyMatrix.h"

namespace FermiOwn {

WoodburyMatrix::WoodburyMatrix( const size_t matSize, const MatCoeffList & coeffList ) :
		mat( matSize, matSize ),
		size( matSize ),
		matNeedsUpdate( true ),
		detNeedsUpdate( true ),
		invNeedsUpdate( true )
{
	mat.setFromTriplets( coeffList.begin(), coeffList.end() );
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

WoodburyMatrix::~WoodburyMatrix() {
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
		} else {
			std::cerr << "Warning: WoodburyMatrix was not able to invert capacitance matrix needed to update inverse matrix!" << std::endl;
		}
		invNeedsUpdate = false;
	}
}

} /* namespace FermiOwn */
