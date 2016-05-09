/*
 * FermiBoolMetropolis.h
 *
 *  Created on: 31.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_UPDATEALGORITHMS_FERMIBOOLMETROPOLIS_H_
#define SRC_UPDATEALGORITHMS_FERMIBOOLMETROPOLIS_H_

#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include "Lattice.h"
#include "FieldBoolean.h"
#include "SlacOperatorMatrix.h"

namespace FermiOwn {

class FermiBoolMetropolis {
public:
	FermiBoolMetropolis( FieldBoolean& boolField, SlacOperatorMatrix& slacMat, const Lattice& lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator );
	virtual ~FermiBoolMetropolis();

	bool updateField();
	bool updateField( size_t x, size_t spin, size_t a, size_t b );

	Complex calculateWeight( int dk, int dntilde );
	bool accept( Complex weight );

	double getHypergeometricFactor( int flavour );
	double getHypergeometricFactor( int n1, int n2 );

	double kappa;
	FieldBoolean & kxiab;
	FieldBoolean oldField;
	SlacOperatorMatrix& slac;
	const Lattice& lat;
	std::ranlux48* rndGen;
	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist;
	std::uniform_int_distribution<int> int2mu_dist;
	std::uniform_int_distribution<int> intNf_dist;
	std::uniform_int_distribution<int> intSpin_dist;

	Eigen::ArrayXi nxOld;
	Eigen::ArrayXi nxNew;
	Eigen::ArrayXi nyOld;
	Eigen::ArrayXi nyNew;

	size_t acceptanceCounter;
	Complex det;
	Complex detOld;
	std::ofstream fWeight;
};

} /* namespace FermiOwn */

#endif /* SRC_UPDATEALGORITHMS_FERMIBOOLMETROPOLIS_H_ */
