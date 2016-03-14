/*
 * SlacOperatorMatrix_test.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>
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
	SlacOperatorMatrix dslac3d( N, N-1, 3 );
	std::cout << "Det=" << dslac3d.det() << std::endl;

	dslac3d.deletePoint( 5 );
	std::cout << dslac3d.getMatrix() << std::endl;
	std::cout << dslac3d.det() << std::endl;
	dslac3d.deletePoint( 6 );
	std::cout << dslac3d.det() << std::endl;
}
