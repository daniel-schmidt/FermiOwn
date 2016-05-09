/*
 * SlacOperatorMatrix_test.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>
#include "FieldBoolean.h"
#include "SlacOperatorMatrix.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;

	if( argc != 2 ) {
		std::cerr << "SlacOperatorMatrix test needs 1 argument: matrix size" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1]);
	if( N > 5 || N < 2 ) {
		std::cerr << "No reference for matrix size " << N << " implemented." << std::endl;
		exit(1);
	}

	SlacOperatorMatrix dslac( N );

	Eigen::Matrix2cd sol2;
	sol2 << 0., 1.5708, -1.5708, 0.;

	Eigen::Matrix3cd sol3;
	sol3 << 0., 1.2092, -1.2092, -1.2092, 0., 1.2092, 1.2092, -1.2092, 0.;

	Eigen::Matrix4cd sol4;
	sol4 << 0., 1.11072, -0.785398, 1.11072,
		  -1.11072, 0., 1.11072, -0.785398,
		   0.785398, -1.11072, 0., 1.11072,
		  -1.11072, 0.785398, -1.11072, 0.;

	Eigen::MatrixXcd sol5(5,5);
	sol5 << 0., 1.06896, -0.660653, 0.660653, -1.06896,
		-1.06896, 0., 1.06896, -0.660653, 0.660653,
		 0.660653, -1.06896, 0., 1.06896, -0.660653,
		 -0.660653, 0.660653, -1.06896, 0., 1.06896,
		 1.06896, -0.660653, 0.660653, -1.06896, 0.;

	Eigen::MatrixXcd sol(N,N);
	switch ( N ) {
	case 2: sol = sol2; break;
	case 3: sol = sol3; break;
	case 4: sol = sol4; break;
	case 5: sol = sol5; break;
	}

	std::cout << dslac.getMatrix() << " should be " << sol << std::endl;
	if(dslac.getMatrix().isApprox(sol, 1e-5) ) std::cout << "Size " << N << " works!" << std::endl;

	std::cout << std::endl << "Testing 3d initialization:" << std::endl;
	size_t Nf = 1;
	size_t dim = 3;
	SlacOperatorMatrix dslac3d( N, N+1, dim, Nf );
	std::cout << "Full determinant:" << dslac3d.det() << std::endl;
//	dslac.setFull();
	for( size_t x = 0; x < N*(N+1)*(N+1); x++ ) {
		dslac3d.erase( x, 0, 0, 0);
		dslac3d.erase( x, 1, 0, 0);
	}
	std::cout << "Deleted everything:" << dslac3d.getMatrix().isIdentity() << std::endl;
	std::cout << "Deleted everything, det=" << dslac3d.det() << std::endl;

	SlacOperatorMatrix dslacNf1( N, 1, dim, 1 );
	std::cout << "Single Flavour: " << std::endl << dslacNf1.getMatrix() << std::endl;
	SlacOperatorMatrix dslacNf2( N, 1, dim, 2 );
	std::cout << "Two Flavour: " << std::endl << dslacNf2.getMatrix() << std::endl;
	std::cout << "Two Flavour determinant: " << std::endl << dslacNf2.det() << std::endl;
	dslacNf2.erase( 0, 0, 0, 1 );
	dslacNf2.erase( 1, 1, 0, 1 );
	dslacNf2.erase( 0, 0, 1, 0 );
	dslacNf2.erase( 1, 1, 1, 0 );
	std::cout << "Two Flavour deleted: " << std::endl << dslacNf2.getMatrix() << std::endl;
	std::cout << "Two Flavour deleted determinant: " << std::endl << dslacNf2.det() << std::endl;
	dslacNf2.setFull();

	std::ranlux48 rndGen;
	FieldBoolean fbool( N, 2, 2, &rndGen, zeroInit );
	fbool.Print();
	fbool.setValue( 1, 0, 0, 0, 1 );
	fbool.setValue( 1, 1, 1, 0, 1 );
	fbool.setValue( 1, 0, 0, 1, 0 );
	fbool.setValue( 1, 1, 1, 1, 0 );
	fbool.Print();
	dslacNf2.erase( fbool );
	std::cout << "Two Flavour deleted: " << std::endl << dslacNf2.getMatrix() << std::endl;
	std::cout << "Two Flavour deleted determinant: " << std::endl << dslacNf2.det() << std::endl;
//	std::cout << std::endl << "Truely erasing cols from slac matrix:" << std::endl;
//	VectorXb kx0 = VectorXb::Constant(N*(N-1)*(N-1), true );
//	VectorXb kx1 = VectorXb::Constant(N*(N-1)*(N-1), true );
//
//	kx0( 2 ) = false;
//	kx0( 4 ) = false;
//
//	dslac3d.setFull();
//	std::cout << "before:" << std::endl << dslac3d.getMatrix() << std::endl;
//	dslac3d.eraseCols( kx0, kx1 );
//	std::cout << "after erasing cols:" << std::endl << dslac3d.getMatrix() << std::endl;
//	kx0( 2 ) = true;
//	kx0( 4 ) = true;
//	kx0( 15 ) = false;
//	kx0( 16 ) = false;
//	dslac3d.eraseRows( kx0, kx1 );
//	std::cout << "after erasing rows:" << std::endl << dslac3d.getMatrix() << std::endl;
//	dslac3d.addPoint( 5 );
//	std::cout << dslac3d.getMatrix() << std::endl;
//	std::cout << "without 6: " << dslac3d.det() << std::endl;
//	dslac3d.addPoint( 6 );
//	std::cout << "full: " << dslac3d.det() << std::endl;
}
