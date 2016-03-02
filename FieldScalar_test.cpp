/*
 * FieldScalar_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include <random>
#include "FieldScalar.h"
#include "Lattice.h"
#include "Action.h"

int main() {

	std::ranlux48 rndGen;

	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	FieldScalar<Real> fs0(lat, rndGen, zeroInit);
	bool correct = true;
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( 0 !=fs0(x) ) correct = false;
	}
	if( !correct ) {
		std::cout << "Initialization with 0s failed! Field has values:" << std::endl;
		fs0.Print();
		exit(1);
	}

	FieldScalar<Complex> fs1(lat, rndGen, oneInit);
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( Complex(1.,0) !=fs1(x) ) correct = false;
	}
	if( !correct ) {
		std::cout << "Initialization with complex 1s failed! Field has values:" << std::endl;
		fs1.Print();
		exit(1);
	} else {
		std::cout << "Initialization with 0 and 1 works!" << std::endl;
	}

	// Random initializations of large vectors
	Nt = 16, Ns = 15, dim = 3;
	Lattice latLarge(Nt, Ns, dim);
	FieldScalar<Complex> fsrnd(latLarge, rndGen, randomInit);
	fsrnd.writeToFile("randomUniform.dat");
	std::cout << "Writing data to randomUniform.dat, check with Octave, if it looks uniformly distributed." << std::endl;

	FieldScalar<Real> fsGauss(latLarge, rndGen, gaussianInit);
	fsGauss.writeToFile("randomGauss1.dat");
	std::cout << "Writing data to randomGauss1.dat, check with Octave, if it looks gaussian distributed." << std::endl;
	fsGauss.setGaussian();
	fsGauss.writeToFile("randomGauss2.dat");
	std::cout << "Writing data to randomGauss2.dat, check with Octave, if it looks gaussian distributed and is different from randomGauss1.dat." << std::endl;
	FieldScalar<Complex> fsGaussC(latLarge, rndGen, gaussianInit);
	fsGaussC.writeToFile("randomGaussC.dat");
	std::cout << "Writing complex data to randomGaussC.dat" << std::endl;
	std::cout << "Check with Octave if both real and imaginary parts are gaussian distributed and independent from the other files." << std::endl;
//	std::cout << std::endl << "Trying +=" << std::endl;
//	fs1 += fsrnd;
//	fs1.Print();
//
//	std::cout << std::endl << "Trying +" << std::endl;
//	(fs1 + fsrnd).Print();
//
//	FieldScalar<Complex> fsmult(lat, oneInit);
//	std::cout << std::endl << "Trying *=Scalar" << std::endl;
//	fsmult *= Complex(0.,1.);
//	fsmult.Print();
//	fsmult /= Complex(0.,1.);
//	fsmult.Print();
//
//	FieldScalar<Complex> fsmult2(lat, oneInit);
//	fsmult2 *= 2.;
//	fsmult *= 4.;
//	std::cout<< "Result of coeff. wise product is " << std::endl;
//	(fsmult*fsmult2).Print();
//	std::cout << " and should be a vector of 8s."<< std::endl;
//	std::cout<< "Result of coeff. wise divide is " << std::endl;
//	(fsmult/fsmult2).Print();
//	std::cout << " and should be a vector of 2."<< std::endl;
//	fsmult2 /= fsmult;
//	std::cout<< "Result of /= should be a vector of 0.5." << std::endl;
//	fsmult2.Print();
//
//	std::cout<<"Scalar product = " << fsmult2.dot(fsmult) << " should be (16,0)" << std::endl;
//	fsmult2 *= Complex(1.,0.5);
//	std::cout<<"Scalar product = " << fsmult.dot(fsmult2) << " should be (16,8)" << std::endl;

}
