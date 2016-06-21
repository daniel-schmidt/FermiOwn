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
#include "ConfigGenerator.h"
#include "WeightFunction.h"

namespace FermiOwn {

class FermiBoolMetropolis {
public:


	FermiBoolMetropolis( FieldBoolean& boolField, const Lattice& lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator );
	virtual ~FermiBoolMetropolis();



	bool updateField();
//	bool updateField( size_t x, size_t spin, size_t a, size_t b );
//	void updateNaive( size_t x );

//	Complex calculateWeight();
//	Complex calculateWeightChange();
//	Complex calculateWeight( int dk, int dntilde );
	bool accept( Complex weight );

	double kappa;
	size_t Nf;
	FieldBoolean & kxiab;
	const Lattice& lat;
	std::ranlux48* rndGen;
	ConfigGenerator confGen;
	WeightFunction weightFun;



	std::uniform_real_distribution<double> uni_real_dist;
	std::uniform_int_distribution<int> intV_dist;
	std::uniform_int_distribution<int> int2mu_dist;
	std::uniform_int_distribution<int> intNf_dist;
	std::uniform_int_distribution<int> intSpin_dist;
	std::uniform_int_distribution<int> intConfIndex;


	FieldBoolean oldField;

	Eigen::ArrayXi nxOld;
	Eigen::ArrayXi nxNew;
	Eigen::ArrayXi nyOld;
	Eigen::ArrayXi nyNew;

	size_t acceptanceCounter;
	Complex det;
	Complex detOld;
	Complex weightOld;
	std::ofstream fWeight;
};

} /* namespace FermiOwn */

#endif /* SRC_UPDATEALGORITHMS_FERMIBOOLMETROPOLIS_H_ */
