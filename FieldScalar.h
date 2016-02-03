/*
 * FieldScalar.h
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#ifndef FIELDSCALAR_H_
#define FIELDSCALAR_H_

#include <iostream>
#include <Eigen/Dense>
#include "Lattice.h"
enum InitType {
	zeroInit,
	oneInit,
	randomInit
};

typedef std::complex< double > Complex;
typedef double Real;

template<class ScalarType> class FieldScalar {
public:
	FieldScalar( const Lattice & newLat, InitType init );
	virtual ~FieldScalar();

	void Print();

private:
	Eigen::Matrix< ScalarType, Eigen::Dynamic, 1> data;
	const Lattice & lat;
};

template<class ScalarType> FieldScalar<ScalarType>::FieldScalar(const Lattice& newLat, InitType init) :
lat(newLat)
{
	data = Eigen::Matrix< ScalarType, Eigen::Dynamic, 1>(lat.getVol());
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

template<class ScalarType> FieldScalar<ScalarType>::~FieldScalar() {}

template<class ScalarType> void FieldScalar<ScalarType>::Print() {
	std::cout << data << std::endl;
}

#endif /* FIELDSCALAR_H_ */
