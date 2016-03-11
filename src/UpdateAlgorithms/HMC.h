/*
 * HMC.h
 *
 *  Created on: 02.03.2016
 *      Author: dschmidt
 */

#ifndef HMC_H_
#define HMC_H_

#include <random>
#include <cmath>

#include "FieldScalar.h"
#include "BasicAction.h"
#include "Integrator.h"


class HMC {
public:
	HMC( double new_t, size_t new_nt, const BasicAction& naction, FieldScalar<Real>& nphi, std::ranlux48* rndGen );
	virtual ~HMC();

	bool update();
private:

	double Hamiltonian();

	std::ranlux48* randomGenerator;
	const BasicAction& act;
	FieldScalar<Real>& phi;
	FieldScalar<Real> momentum;
	Integrator integrator;
	std::uniform_real_distribution<Real> uniformDistribution01;
};

#endif /* HMC_H_ */