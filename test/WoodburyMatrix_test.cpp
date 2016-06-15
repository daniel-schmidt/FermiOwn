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
	std::cout << "Inverse:" << std::endl;
	mat.PrintInverse();

	SparseMat colUpdate( 4, 4 );
	MatCoeffList colUpCoeffs(4);
	colUpCoeffs[0] = Eigen::Triplet<Complex> ( 1, 0 , 1.);
	colUpCoeffs[1] = Eigen::Triplet<Complex> ( 1, 2 , 1.);
	colUpCoeffs[2] = Eigen::Triplet<Complex> ( 2, 1 , 1.);
	colUpCoeffs[3] = Eigen::Triplet<Complex> ( 2, 3 , 1.);
	colUpdate.setFromTriplets( colUpCoeffs.begin(), colUpCoeffs.end() );

	SparseMat rowUpdate( 4, 4 );
	MatCoeffList rowUpCoeffs(4);
	rowUpCoeffs[0] = Eigen::Triplet<Complex> ( 0, 1 , Complex( -1.5, -2 ) );
	rowUpCoeffs[1] = Eigen::Triplet<Complex> ( 1, 2 , Complex( -1.5, -4 ) );
	rowUpCoeffs[2] = Eigen::Triplet<Complex> ( 2, 2 , 1.);
	rowUpCoeffs[3] = Eigen::Triplet<Complex> ( 3, 1 , 1.);
	rowUpdate.setFromTriplets( rowUpCoeffs.begin(), rowUpCoeffs.end() );

	std::cout << "Performing update with " << std::endl << colUpdate << std::endl << "and" << std::endl << rowUpdate << std::endl;

	mat.setUpdateMatrices( colUpdate, rowUpdate );
	mat.updateDet();

	std::cout << "Det after update: " << mat.getDet() << " and should be (-2.25, -9)" << std::endl;

	mat.updateMatrix();
	std::cout << "Matrix after update: " << std::endl;
	mat.Print();

	mat.updateInverse();
	std::cout << "Inverse after update: " << std::endl;
	mat.PrintInverse();

	std::cout << " and should be " << std::endl << "(0.666667,0) 0 0 0\n0 (0,0) (1,0) 0\n0 (1,0) (0,0) 0\n0 0 0 (0.0392157,-0.156863)" << std::endl;

	std::cout << "Testing self-adjoint initialization." << std::endl;
	coeffList.clear();
	for( size_t i = 0; i < 4; i++ ) {
		coeffList.push_back( Eigen::Triplet<Complex>( 0, i, Complex( 1.5, 2.*i ) ) );
	}
	coeffList.push_back( Eigen::Triplet<Complex>( 1, 2, Complex( 0, 5 ) ) );
	WoodburyMatrix selfAdjointMat( 4, coeffList, true );
	selfAdjointMat.Print();
	std::cout << "Det of selfadjoint matrix is " << selfAdjointMat.getDet() << " and should be 956.25"<< std::endl;
}
