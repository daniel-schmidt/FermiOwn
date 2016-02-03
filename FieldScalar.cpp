/*
 * FieldScalar.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include <iostream>
#include "FieldScalar.h"

FieldScalar::FieldScalar(const Lattice& newLat, InitType init) :
lat(newLat)
{
	data = Eigen::Matrix<double, Eigen::Dynamic, 1>(lat.getVol());
	switch( init ) {
	case zeroInit:
		data.setZero();
		return;
	case oneInit:
		data.setOnes();
		return;
	case randomInit:
		data.setRandom();
		return;
	}

}

FieldScalar::~FieldScalar() {}

void FieldScalar::Print() {
	std::cout << data << std::endl;
}
