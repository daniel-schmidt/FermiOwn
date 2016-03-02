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
	HMC( std::ranlux48& rndGen );
	virtual ~HMC();

	void update();

private:

	double Hamiltonian();

	BasicAction& act;
	Integrator& integrator;
	FieldScalar<Real>& phi;
	FieldScalar<Real> momentum;
	std::ranlux48& randomGenerator;
	std::uniform_real_distribution<Real> uniformDistribution01;
};

#endif /* HMC_H_ */
