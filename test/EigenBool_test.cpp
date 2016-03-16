/*
 * EigenBool_test.cpp
 *
 *  Created on: 15.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>

int linIndex(int N, int i, int j) {
	if( i == j) {
		return i;
	} else {
		if( i > j ) {
			int tmp = j;
			j = i;
			i = tmp;
		}
		return N+N*(N-1)/2 - (N-i)*(N-i-1)/2+j-i-1;
	}
}

int main( int argc, char** argv ) {
	if( argc != 2 ) {
		std::cerr << "Need 1 argument: vector size" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1]);

	using namespace Eigen;

	typedef Matrix< bool, 1, Dynamic > VectorXb;
	typedef Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > MatrixXb;

	VectorXb vecTrue = VectorXb::Constant(N,true);
	std::cout << vecTrue << std::endl;

	VectorXb vecFalse = VectorXb::Constant(N,false);
	std::cout << vecFalse << std::endl;

	//	std::cout << vecTrue*vecFalse << std::endl;
	//	std::cout << vecTrue*vecTrue << std::endl;
	//	std::cout << vecFalse*vecFalse << std::endl;
	std::cout << vecFalse.dot(vecTrue) << std::endl;
	std::cout << vecTrue.dot(vecTrue) << std::endl;

	std::cout << vecTrue+vecFalse << std::endl;
	std::cout << vecTrue+vecTrue << std::endl;
	std::cout << vecFalse + vecFalse << std::endl;
	std::cout << vecFalse - vecFalse << std::endl;
	std::cout << vecTrue - vecTrue << std::endl;
	std::cout << vecTrue.count() << std::endl;
	vecTrue(0) = false;
	std::cout << vecTrue.count() << std::endl;
	std::cout << vecTrue.size() << std::endl;

	std::cout << std::endl << "Testing matrix of bools..." << std::endl;
	MatrixXb mat = MatrixXb::Constant(N,N*(N+1)/2,true);
	mat(1,1) = false;
	mat(2,0) = false;
	mat(2,1) = false;
	std::cout << mat << std::endl;
	std::cout << mat.count() << std::endl;
	std::cout << mat.leftCols(2).rowwise().count() << std::endl;

	std::cout << std::endl << "Test index setting..." << std::endl;
	// set all entries with a=3 to false:
	for( size_t a = 0; a < N; a++ ) {
		mat(3,linIndex(N, a, 3)) = false;
	}

	// set all entries with a=2 to false:
	for( size_t a = 0; a < N; a++ ) {
		mat(4,linIndex(N, a, 2)) = false;
	}

	std::cout << mat.row(3) << std::endl;
	std::cout << mat.row(4) << std::endl;
}
