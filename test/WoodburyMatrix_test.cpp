/*
 * WoodburyMatrix_test.cpp
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/SparseCore>
#include "WoodburyMatrix.h"

int main( int argc, char** argv ) {

	using namespace FermiOwn;

	MatCoeffList coeffList;

	for( size_t i = 0; i < 4; i++ ) {
		coeffList.push_back( Eigen::Triplet<Complex>( i, i, Complex( 1.5, 2.*i ) ) );
	}

	WoodburyMatrix mat( 4, coeffList );
	mat.Print();
	std::cout << "Determinant: " << mat.getDet() << " should agree with (-93.9375, -31.5)" << std::endl;
	mat.PrintInverse();
}
