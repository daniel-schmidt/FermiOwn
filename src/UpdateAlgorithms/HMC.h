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

#include "Field.h"
#include "BasicAction.h"
#include "Integrator.h"

namespace FermiOwn {

class HMC {
public:
	HMC( double new_t, size_t new_nt, const BasicAction& naction, Field<Real>& nphi, std::ranlux48* rndGen );
	virtual ~HMC();

	bool update();
private:

	double Hamiltonian();

	std::ranlux48* randomGenerator;
	const BasicAction& act;
	Field<Real>& phi;
	Field<Real> momentum;
	Integrator integrator;
	std::uniform_real_distribution<Real> uniformDistribution01;
};

} // namespace FermiOwn

#endif /* HMC_H_ */
