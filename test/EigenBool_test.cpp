/*
 * EigenBool_test.cpp
 *
 *  Created on: 15.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>

int main( int argc, char** argv ) {
	if( argc != 2 ) {
		std::cerr << "Need 1 argument: vector size" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1]);

	using namespace Eigen;

	typedef Matrix< bool, 1, Dynamic > VectorXb;

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
}
