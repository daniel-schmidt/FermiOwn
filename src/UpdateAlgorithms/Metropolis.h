/*
 * Metropolis.h
 *
 *  Created on: 10.03.2016
 *      Author: dschmidt
 */

#ifndef METROPOLIS_H_
#define METROPOLIS_H_

#include <random>
#include <cmath>

#include "Field.h"
#include "BasicAction.h"

namespace FermiOwn {

class Metropolis {
public:
	/**
	 *
	 * @param ndelta is the variation witdh for newly proposed configurations:
	 * 		  phi[n+1]=phi[n] + delta[n], where delta[n] is drawn from the interval [-delta,delta]
	 * @param naction
	 * @param nphi
	 * @param rndGen
	 */
	Metropolis( const double ndelta, const BasicAction& naction, Field<Real>& nphi, std::ranlux48* rndGen );
	virtual ~Metropolis();

	/** @brief Lattice sweep
	 *
	 * @return acceptance rate of this sweep
	 */
	double update();
private:
	const double delta;
	std::ranlux48* randomGenerator;
	const BasicAction& act;
	Field<Real>& phi;
	std::uniform_real_distribution<Real> uniformDistribution01;
};

} // namespace FermiOwn

#endif /* METROPOLIS_H_ */
