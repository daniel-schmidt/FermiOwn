/*
 * Metropolis.cpp
 *
 *  Created on: 10.03.2016
 *      Author: dschmidt
 */

#include "Metropolis.h"

namespace FermiOwn {

Metropolis::Metropolis( const double ndelta, const BasicAction& naction, Field<Real>& nphi, std::ranlux48* rndGen) :
	delta(ndelta),
	randomGenerator(rndGen),
	act(naction),
	phi(nphi)
{}

Metropolis::~Metropolis() {
	// TODO Auto-generated destructor stub
}

double Metropolis::update() {
	size_t accepted = 0;
	for( size_t x = 0; x < phi.getSize(); x++ ) {

		// store old field and value of action
		double dAction = act.getAction(phi);
		Field<Real> old = phi;

		// generate random number in the interval [-delta, delta] and update field
		double dphi = delta*( 2*uniformDistribution01(*randomGenerator)-1);
		phi(x) += dphi;

		// calculate action difference
		dAction -= act.getAction(phi);

		// accept-reject step
		double r = uniformDistribution01(*randomGenerator);
//		std::cout << std::exp(dAction);
		if( r < std::exp(dAction) ) {
			accepted++;
//			std::cout << "\t1" << std:: endl;
		} else {
			phi = old;
//			std::cout << "\t0" << std:: endl;
		}
	}
	return double(accepted)/phi.getSize();
}

} // namespace FermiOwn
