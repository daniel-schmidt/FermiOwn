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
		coeffList.push_back( CoeffTriplet( i, i, Complex( 1.5, 2.*i ) ) );
	}

	WoodburyMatrix mat( 4, coeffList );
	mat.Print();
	std::cout << "Determinant: " << mat.getDet() << " should agree with (-93.9375, -31.5)" << std::endl;
	std::cout << "Inverse:" << std::endl;
	mat.PrintInverse();

	SparseMat colUpdate( 4, 4 );
	MatCoeffList colUpCoeffs(4);
	colUpCoeffs[0] = CoeffTriplet ( 1, 0 , 1.);
	colUpCoeffs[1] = CoeffTriplet ( 1, 2 , 1.);
	colUpCoeffs[2] = CoeffTriplet ( 2, 1 , 1.);
	colUpCoeffs[3] = CoeffTriplet ( 2, 3 , 1.);
	colUpdate.setFromTriplets( colUpCoeffs.begin(), colUpCoeffs.end() );

	SparseMat rowUpdate( 4, 4 );
	MatCoeffList rowUpCoeffs(4);
	rowUpCoeffs[0] = CoeffTriplet ( 0, 1 , Complex( -1.5, -2 ) );
	rowUpCoeffs[1] = CoeffTriplet ( 1, 2 , Complex( -1.5, -4 ) );
	rowUpCoeffs[2] = CoeffTriplet ( 2, 2 , 1.);
	rowUpCoeffs[3] = CoeffTriplet ( 3, 1 , 1.);
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
		coeffList.push_back( CoeffTriplet( 0, i, Complex( 1.5, 2.*i ) ) );
	}
	coeffList.push_back( CoeffTriplet( 1, 2, Complex( 0, 5 ) ) );
	WoodburyMatrix selfAdjointMat( 4, coeffList, true );
	selfAdjointMat.Print();
	std::cout << "Det of selfadjoint matrix is " << selfAdjointMat.getDet() << " and should be 956.25"<< std::endl;

	std::cout << std::endl << "Testing delete-only update." << std::endl;

	coeffList.clear();
	coeffList.push_back( CoeffTriplet( 0, 3, 3 ) );
	coeffList.push_back( CoeffTriplet( 1, 1, 5 ) );
	coeffList.push_back( CoeffTriplet( 2, 0, 2 ) );
	coeffList.push_back( CoeffTriplet( 3, 2, 7 ) );
	coeffList.push_back( CoeffTriplet( 3, 3, 8 ) );
	WoodburyMatrix delOnlyMat( 4, coeffList );
	delOnlyMat.Print();
	std::cout << "Initial determinant is " << delOnlyMat.getDet() << " and should be " << 210 << std::endl;

	// update matrices calculated by mathematica
	colUpCoeffs[0] = CoeffTriplet ( 0, 0 , 1);
	colUpCoeffs[1] = CoeffTriplet ( 0, 3 , 1.);
	colUpCoeffs[2] = CoeffTriplet ( 2, 1 , 1.);
	colUpCoeffs[3] = CoeffTriplet ( 2, 2 , 1.);
	colUpCoeffs.push_back( CoeffTriplet( 3, 3, -7. ) );
	colUpdate.setFromTriplets( colUpCoeffs.begin(), colUpCoeffs.end() );

	rowUpCoeffs[0] = CoeffTriplet ( 0, 3 , -3. );
	rowUpCoeffs[1] = CoeffTriplet ( 1, 0 , -2. );
	rowUpCoeffs[2] = CoeffTriplet ( 2, 0 , 1.);
	rowUpCoeffs[3] = CoeffTriplet ( 3, 2 , 1.);
	rowUpdate.setFromTriplets( rowUpCoeffs.begin(), rowUpCoeffs.end() );

	std::cout << "Performing update with " << std::endl << colUpdate << std::endl << "and" << std::endl << rowUpdate << std::endl;

	std::vector< size_t > rowsToDel = { 0, 2 };
	std::vector< size_t > colsToDel = { 2, 0 };

	delOnlyMat.setUpdateMatricesDeleteOnly( colUpdate, rowUpdate, rowsToDel, colsToDel );
	delOnlyMat.updateDet();

	std::cout << "Det after update: " << delOnlyMat.getDet() << " and should be -40" << std::endl;

	delOnlyMat.updateMatrix();
	std::cout << "Matrix after update: " << std::endl;
	delOnlyMat.Print();

	delOnlyMat.updateInverse();
	std::cout << "Inverse after update: " << std::endl;
	delOnlyMat.PrintInverse();
	std::cout << " and should have diagonal entries 0.2 at(1,1) and 0.125 at (3,3) and ones a (0,2) and (2,0)" << std::endl;
	std::cout << "Matrix times inverse: " << std::endl << delOnlyMat.getMatrix() * delOnlyMat.getInverse() << std::endl;
}
