/*
 * FieldScalar_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include "FieldScalar.h"
#include "Lattice.h"
#include "Action.h"

int main() {
	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	FieldScalar<Real> fs0(lat, zeroInit);
	fs0.Print();
	FieldScalar<Complex> fs1(lat, oneInit);
	fs1.Print();
	FieldScalar<Complex> fsrnd(lat, randomInit);
	fsrnd.Print();

	std::cout << std::endl << "Trying +=" << std::endl;
	fs1 += fsrnd;
	fs1.Print();

	std::cout << std::endl << "Trying +" << std::endl;
	(fs1 + fsrnd).Print();

	FieldScalar<Complex> fsmult(lat, oneInit);
	std::cout << std::endl << "Trying *=Scalar" << std::endl;
	fsmult *= Complex(0.,1.);
	fsmult.Print();
	fsmult /= Complex(0.,1.);
	fsmult.Print();

	FieldScalar<Complex> fsmult2(lat, oneInit);
	fsmult2 *= 2.;
	fsmult *= 4.;
	std::cout<< "Result of coeff. wise product is " << std::endl;
	(fsmult*fsmult2).Print();
	std::cout << " and should be a vector of 8s."<< std::endl;
	std::cout<< "Result of coeff. wise divide is " << std::endl;
	(fsmult/fsmult2).Print();
	std::cout << " and should be a vector of 2."<< std::endl;
	fsmult2 /= fsmult;
	std::cout<< "Result of /= should be a vector of 0.5." << std::endl;
	fsmult2.Print();

	std::cout<<"Scalar product = " << fsmult2.dot(fsmult) << " should be (16,0)" << std::endl;
	fsmult2 *= Complex(1.,0.5);
	std::cout<<"Scalar product = " << fsmult.dot(fsmult2) << " should be (16,8)" << std::endl;

	FieldScalar<Real> phi(lat, oneInit);
	double lambda = 0.5;
	double kappa = 0.1;
	double Spot = kappa*phi.dot(phi) + lambda*(phi*phi).dot(phi*phi);
	std::cout << "Potential is: " << Spot  << " and should be 4.8" << std::endl;

	Action act( phi, kappa, lambda);
	std::cout << "Action is: " << act.getAction()  << " and should be 3.2" << std::endl;
	std::cout << "Force is: ";
	act.getForce().Print();
	std::cout<< " and should be vector of 14.8s." << std::endl;
	lambda = 0.9;
	kappa = 0.3;
	Action act2( phi, kappa, lambda );
	std::cout << "Action is: " << act2.getAction()  << " and should be -6.4" << std::endl;



}
