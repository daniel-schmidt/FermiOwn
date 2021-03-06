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
#include <set>
#include "Lattice.h"
#include "FieldBoolean.h"
#include "MetropolisStep.h"
#include "ConfigPerPointGenerator.h"
#include "WeightFunction.h"


namespace FermiOwn {

class FermiBoolMetropolis: public MetropolisStep {
public:
	FermiBoolMetropolis( FieldBoolean& boolField, const Lattice& lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator );
	virtual ~FermiBoolMetropolis();

	void initializeField();

protected:

	virtual void propose();
	virtual double change();
	virtual void accept();
	virtual void reject();

	double kappa;
	size_t Nf;
	FieldBoolean & kxiab;
	FieldBoolean oldField;
	const Lattice& lat;

	std::uniform_int_distribution<int> intV_dist;
	std::set<size_t> changedPoints;

	WeightFunction weightFun;
	ConfigPerPointGenerator confGen;
};

//class FermiBoolMetropolis {
//public:
//
//
//	FermiBoolMetropolis( FieldBoolean& boolField, const Lattice& lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator );
//	virtual ~FermiBoolMetropolis();
//
//	void initializeField();
//
//	bool updateField();
//	//	bool updateField( size_t x, size_t spin, size_t a, size_t b );
//	//	void updateNaive( size_t x );
//
//	//	Complex calculateWeight();
//	//	Complex calculateWeightChange();
//	//	Complex calculateWeight( int dk, int dntilde );
//	bool accept( Complex weight );
//
//	double kappa;
//	size_t Nf;
//	FieldBoolean & kxiab;
//	const Lattice& lat;
//	std::ranlux48* rndGen;
//	ConfigGenerator confGen;
//	WeightFunction weightFun;
//
//	std::uniform_real_distribution<double> uni_real_dist;
//	std::uniform_int_distribution<int> intV_dist;
//	std::uniform_int_distribution<int> int2mu_dist;
//	std::uniform_int_distribution<int> intNf_dist;
//	std::uniform_int_distribution<int> intSpin_dist;
//	std::uniform_int_distribution<int> intConfIndex;
//
//	FieldBoolean oldField;
//
//	size_t acceptanceCounter;
//	std::ofstream fWeight;
//
//	double phase;
//	Complex expPhase;
//};

} /* namespace FermiOwn */

#endif /* SRC_UPDATEALGORITHMS_FERMIBOOLMETROPOLIS_H_ */
