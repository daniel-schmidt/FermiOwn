/*
 * SlacOperatorMatrix.h
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_SLACOPERATORMATRIX_H_
#define SRC_ACTIONS_SLACOPERATORMATRIX_H_

#include <Eigen/Dense>
#include "Constants.h"

namespace FermiOwn {

class SlacOperatorMatrix {
public:
	SlacOperatorMatrix( size_t size );
	virtual ~SlacOperatorMatrix();

	const Eigen::MatrixXd getMatrix() const;

private:
	size_t N;
	Eigen::MatrixXd dslac;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_SLACOPERATORMATRIX_H_ */
