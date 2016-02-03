/*
 * FieldScalar.h
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#ifndef FIELDSCALAR_H_
#define FIELDSCALAR_H_

#include <Eigen/Dense>
#include "Lattice.h"

enum InitType {
	zeroInit,
	oneInit,
	randomInit
};

class FieldScalar {
public:
	FieldScalar( const Lattice & newLat, InitType init );
	virtual ~FieldScalar();

	void Print();

private:
	Eigen::Matrix<double, Eigen::Dynamic, 1> data;
	const Lattice & lat;
};

#endif /* FIELDSCALAR_H_ */
