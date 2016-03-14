/*
 * CliffordAlgebra.h
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#ifndef CLIFFORDALGEBRA_H_
#define CLIFFORDALGEBRA_H_

#include <vector>
#include <Eigen/Dense>
#include "Constants.h"

namespace FermiOwn {

// TODO: make this templated by the dimenson of the algebra
class CliffordAlgebra {
public:
	CliffordAlgebra();
	virtual ~CliffordAlgebra();

	inline const std::vector< Eigen::Matrix2cd > getGammas() const;

private:
	std::vector< Eigen::Matrix2cd > gamma;
};

CliffordAlgebra::CliffordAlgebra() {
	Eigen::Matrix2cd mat;
	mat << 0., 1., 1., 0.;
	gamma.push_back(mat);
	mat << 0., -I, I, 0.;
	gamma.push_back(mat);
	mat << 1., 0., 0., -1.;
	gamma.push_back(mat);
}

CliffordAlgebra::~CliffordAlgebra() {};

const std::vector< Eigen::Matrix2cd > CliffordAlgebra::getGammas() const {
	return gamma;
}


} /* namespace FermiOwn */

#endif /* CLIFFORDALGEBRA_H_ */
