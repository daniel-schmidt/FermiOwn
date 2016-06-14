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
		size( matSize )
{
	mat.setFromTriplets( coeffList.begin(), coeffList.end() );

	Eigen::SparseLU< SparseMat > solver;	//TODO: check, if we can use an SPD solver here for SLAC
	solver.compute( mat );
	det =  solver.determinant();

	SparseMat id( size, size );
	id.setIdentity();
	inv = solver.solve( id );
}

WoodburyMatrix::~WoodburyMatrix() {
	// TODO Auto-generated destructor stub
}

} /* namespace FermiOwn */
