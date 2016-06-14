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
		validState( false )
{
	mat.setFromTriplets( coeffList.begin(), coeffList.end() );

	Eigen::SparseLU< SparseMat > solver;	//TODO: check, if we can use an SPD solver here for SLAC
	solver.compute( mat );
	if( solver.info() != Eigen::Success ) {
		std::cerr << "WoodburyMatrix cannot calculate matrix decomposition to obtain determinant." << std::endl;
		exit(1);
	}

	det =  solver.determinant();

	SparseMat id( size, size );
	id.setIdentity();
	inv = solver.solve( id );

	if( solver.info() != Eigen::Success ) {
		std::cerr << "WoodburyMatrix failed to invert the matrix." << std::endl;
		exit(1);
	}
	validState = true;
}

WoodburyMatrix::~WoodburyMatrix() {
}


void WoodburyMatrix::updateMatrix() {
	mat += U * V;
}

Complex WoodburyMatrix::updateDet() {
	VTimesInv = V*inv;
	smallId.resize( V.rows(), U.cols() );
	smallId.setIdentity();
	capacitanceMatrixLU.compute( VTimesInv * U + smallId );
	Complex detChange;
	if( capacitanceMatrixLU.info() == Eigen::Success ) {
		detChange = capacitanceMatrixLU.determinant();
	} else {
		detChange = 0.;
	}
	det *= detChange;
	return detChange;
}

void WoodburyMatrix::updateInverse() {
	SparseMat capInv = capacitanceMatrixLU.solve( smallId );
	if( capacitanceMatrixLU.info() == Eigen::Success ) {
		inv -= inv * U * capInv * VTimesInv;
	} else {
		std::cerr << "Warning: WoodburyMatrix was not able to invert capacitance matrix needed to update inverse matrix!" << std::endl;
	}
}

} /* namespace FermiOwn */
