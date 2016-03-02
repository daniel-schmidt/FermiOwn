/*
 * Integrator_test.cpp
 *
 *  Created on: 01.03.2016
 *      Author: dschmidt
 */

#include <random>
#include "Lattice.h"
#include "BasicAction.h"
#include "Action.h"
#include "Integrator.h"

class HarmonicAction : public BasicAction {
public:
	HarmonicAction(){};
	~HarmonicAction(){};

	double getAction( const FieldScalar<Real>& phi ) const {
		return (-1.)*phi.dot(phi);
	}

	FieldScalar<Real> getForce( const FieldScalar<Real>& phi ) const {
		return phi;
	}
};


int main() {
	std::cout << "Testing the integrator for the HMC algorithm..." << std::endl;

	std::ranlux48 rndGen;

	size_t Nt = 2, Ns = 2, dim = 3;
	double kappa = 0.5, lambda = 0.1;
	Lattice lat(Nt, Ns, dim);

	FieldScalar<Real> fs0(lat, rndGen, zeroInit);
	Action act(kappa, lambda);
	fs0.Print();

	// time reversal test
	// TODO: automate this test
	double t = 0.6;
	double nt = 5;
	Integrator integrator(fs0, act, rndGen, t, nt);
	integrator.integrate();
	std::cout << "After first trajectory:" << std::endl;
	fs0.Print();

	integrator.integrate();
	std::cout << "After second trajectory:" << std::endl;
	fs0.Print();

	std::cout << "Performing time reversal check, result should be the same as before." << std::endl;
	integrator.invertMomentum();
	integrator.integrate();
	fs0.Print();

	std::cout << "Another step back should lead to zero: " << std::endl;
	integrator.integrate();
	fs0.Print();

	std::cout << "Trying harmonic action..." << std::endl;
	HarmonicAction hact;
	nt=100;
	size_t imax = 50;
	t=25.1327/imax; // 8 pi divided in imax parts, each nt integration steps inbetween
	Integrator hint(fs0, hact, rndGen, t, nt);
	std::cout << "data={";
	for( size_t i = 0; i < imax; i++ )
	{
		std::cout << "{" << i*t << ", " << fs0(0) << "}," << std::endl;
		hint.integrate();
	}
	std::cout << "{" << imax*t << ", " << fs0(0) << "}}" << std::endl;

}
