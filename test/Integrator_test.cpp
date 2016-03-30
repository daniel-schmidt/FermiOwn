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

namespace FermiOwn {
class HarmonicAction : public BasicAction {

public:
	HarmonicAction(){};
	~HarmonicAction(){};

	double getAction( const Field<Real>& phi ) const {
		return phi.cwiseMultAndSum(phi);
	}

	Field<Real> getForce( const Field<Real>& phi ) const {
		return phi;
	}
};
} // namespace FermiOwn

int main() {
	using namespace FermiOwn;
	std::cout << "Testing the integrator for the HMC algorithm..." << std::endl;

	std::ranlux48 rndGen;

	size_t Nt = 2, Ns = 2, dim = 3;
	double kappa = 0.5, lambda = 0.1;
	Lattice lat(Nt, Ns, dim);

	Field<Real> fs0(lat.getVol(), 1, &rndGen, zeroInit);
	Field<Real> pRand(lat.getVol(), 1, &rndGen, gaussianInit);
	Action act(lat, kappa, lambda);
	fs0.Print();

	// time reversal test
	// TODO: automate this test
	double t = 0.6;
	double nt = 5;
	Integrator integrator(fs0, pRand, act, t, nt);
	integrator.integrate();
	std::cout << "After first trajectory:" << std::endl;
	fs0.Print();

	integrator.integrate();
	std::cout << "After second trajectory:" << std::endl;
	fs0.Print();

	std::cout << "Performing time reversal check, result should be the same as before." << std::endl;
	pRand *= -1.;
	integrator.integrate();
	fs0.Print();

	std::cout << "Another step back should lead to zero: " << std::endl;
	integrator.integrate();
	fs0.Print();

	std::cout << "Trying harmonic action..." << std::endl;
	HarmonicAction hact;
	Field<Real> pOnes(lat.getVol(), 1, &rndGen, oneInit);
	nt=5;
	size_t imax = 200;
	t=25.1327/imax; // 8 pi divided in imax parts, each nt integration steps inbetween
	Integrator hint(fs0, pOnes, hact, t, nt);
	std::cout << "data={";
	for( size_t i = 0; i < imax; i++ )
	{
//		std::cout << "{" << i*t << ", " << fs0(0) << "}," << std::endl;
		std::cout << (pOnes.cwiseMultAndSum(pOnes)+hact.getAction(fs0)) << "," << std::endl;
		hint.integrate();
	}
//	std::cout << "{" << imax*t << ", " << fs0(0) << "}}" << std::endl;
	std::cout << pOnes.cwiseMultAndSum(pOnes) + hact.getAction(fs0) << "};" << std::endl;

}
