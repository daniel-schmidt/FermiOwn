/*
 * CliffordAlgebra_test.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include "CliffordAlgebra.h"

int main() {
	using namespace FermiOwn;
	CliffordAlgebra cliff;
	auto gamma = cliff.getGammas();
	std::cout << gamma[0] * gamma[0] << std::endl;
	std::cout << gamma[1] * gamma[1] << std::endl;
	std::cout << gamma[2] * gamma[2] << std::endl;
	std::cout << gamma[0] * gamma[1] << std::endl;
	std::cout << gamma[1] * gamma[2] << std::endl;
	std::cout << gamma[0] * gamma[2] << std::endl;

	size_t passed = 0;
	if( (gamma[0]*gamma[0]).isApprox(Eigen::Matrix2cd::Identity())) passed++;
	if( (gamma[1]*gamma[1]).isApprox(Eigen::Matrix2cd::Identity())) passed++;
	if( (gamma[2]*gamma[2]).isApprox(Eigen::Matrix2cd::Identity())) passed++;
	std::cout << passed << " matrices squared to unity." << std::endl;

	passed = 0;
	if( (gamma[0]*gamma[1]).isApprox(I*gamma[2])) passed++;
	if( (gamma[2]*gamma[0]).isApprox(I*gamma[1])) passed++;
	if( (gamma[1]*gamma[2]).isApprox(I*gamma[0])) passed++;
	std::cout << passed << " correct products." << std::endl;

	passed = 0;
	if( (gamma[0] * gamma[1] + gamma[1] * gamma[0]).isApprox( Eigen::Matrix2cd::Zero(), 1e-5 ) ) passed++;
	if( (gamma[0] * gamma[2] + gamma[2] * gamma[0]).isApprox( Eigen::Matrix2cd::Zero(), 1e-5 ) ) passed++;
	if( (gamma[1] * gamma[2] + gamma[2] * gamma[1]).isApprox( Eigen::Matrix2cd::Zero(), 1e-5 ) ) passed++;
	std::cout << passed << " commutator relations worked." << std::endl;
}
