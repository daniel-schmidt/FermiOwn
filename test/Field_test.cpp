/*
 * Field_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include <random>

#include "Field.h"
#include "FieldBoolean.h"
#include "Lattice.h"

int main() {
	using namespace FermiOwn;

	std::ranlux48 rndGen;

	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	Field<Real> fs0(lat.getVol(), 1, &rndGen, zeroInit);
	bool correct = true;
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( 0 !=fs0(x) ) correct = false;
	}
	if( !correct ) {
		std::cout << "Initialization with 0s failed! Field has values:" << std::endl;
		fs0.Print();
		exit(1);
	}

	Field<Complex> fs1(lat.getVol(), 1, &rndGen, oneInit);
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

	Field<Complex> fs1copy( fs1 );
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( Complex(1.,0) !=fs1copy(x) ) correct = false;
	}
	if( !correct ) {
		std::cout << "Initialization by copy failed! Field has values:" << std::endl;
		fs1copy.Print();
		exit(1);
	} else {
		std::cout << "Initialization by copy works!" << std::endl;
	}

	Field<Complex> fsSum(lat.getVol(), 1, &rndGen, randomInit);
	fsSum = fs1+fs1copy;
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( Complex(2.,0) !=fsSum(x) ) correct = false;
	}
	if( !correct ) {
		std::cout << "Copy assignment or sum failed! Field has values:" << std::endl;
		fsSum.Print();
		exit(1);
	} else {
		std::cout << "Copy assignment works!" << std::endl;
	}

	// Random initializations of large vectors
	Nt = 16, Ns = 15, dim = 3;
	Lattice latLarge(Nt, Ns, dim);
	Field<Complex> fsrnd(latLarge.getVol(), 1, &rndGen, randomInit);
	fsrnd.writeToFile("randomUniform.dat");
	std::cout << "Writing data to randomUniform.dat, check with Octave, if it looks uniformly distributed." << std::endl;

	Field<Real> fsGauss(latLarge.getVol(), 1, &rndGen, gaussianInit);
	fsGauss.writeToFile("randomGauss1.dat");
	std::cout << "Writing data to randomGauss1.dat, check with Octave, if it looks gaussian distributed." << std::endl;
	fsGauss.setGaussian();
	fsGauss.writeToFile("randomGauss2.dat");
	std::cout << "Writing data to randomGauss2.dat, check with Octave, if it looks gaussian distributed and is different from randomGauss1.dat." << std::endl;
	Field<Complex> fsGaussC(latLarge.getVol(), 1, &rndGen, gaussianInit);
	fsGaussC.writeToFile("randomGaussC.dat");
	std::cout << "Writing complex data to randomGaussC.dat" << std::endl;
	std::cout << "Check with Octave if both real and imaginary parts are gaussian distributed and independent from the other files." << std::endl;

	std::cout << std::endl << "Testing boolean field..." << std::endl;
	FieldBoolean fbool1( 2*3*3, 2, 3, &rndGen, oneInit );
	fbool1.Print();
	FieldBoolean fbool0( 2*3*3, 2, 3, &rndGen, zeroInit );
	fbool0.Print();
	std::cout << "Inverting, should be at row 3, col 6 and 10..." << std::endl;
	fbool1.setValue( false, 2, 1, 0, 1 );
	fbool1.setValue( false, 2, 0, 2, 1 );
	fbool1.Print();
	std::cout << "Inverting, should be at row 7, col 2..." << std::endl;
	fbool0.invert(6,0,1,1);
	fbool0.Print();

	std::cout << "Testing count, field 0 at x=6: " << std::endl << fbool0.countSummedSpin(6) << std::endl << " field 1 at x=3: " << std::endl << fbool1.countSummedSpin(3) << std::endl;

	std::cout << std::endl << "Testing copy assignment:" << std::endl;
	fbool1 = fbool0;
	fbool1.Print();

	return 0;
}
