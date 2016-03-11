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

#include "FieldScalar.h"
#include "BasicAction.h"
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
	Metropolis( const double ndelta, const BasicAction& naction, FieldScalar<Real>& nphi, std::ranlux48* rndGen );
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
	FieldScalar<Real>& phi;
	std::uniform_real_distribution<Real> uniformDistribution01;
};

#endif /* METROPOLIS_H_ */
