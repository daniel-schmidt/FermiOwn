/*
 * WeightFunction_test.cpp
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#include "../src/Actions/ThirringWeightFunction.h"

#include <set>
#include "ThirringKField.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;

	size_t Nt = 4;
	size_t Ns = 1;
	size_t dim = 3;
	size_t dimSpinor = 2;
	size_t Nf = 2;

	RowVectorXb conf96( 8 );
	conf96 << 0, 1, 0, 0, 0, 0, 1, 0;

	RowVectorXb row0( 8 );
	row0 << 0, 1, 0, 0, 0, 0, 1, 0;

	RowVectorXb row1( 8 );
	row1 << 0, 1, 1, 0, 0, 1, 1, 0;

	RowVectorXb row2( 8 );
	row2 << 0, 0, 1, 0, 0, 1, 0, 0;

	std::vector<Complex> testConf( 20 );
	double dlambda = 0.1;

	std::cout << "Testing weight for Configuration 96" << std::endl << "lambda\tRe(weight)\tIm(weight)" << std::endl;
	for( int i = 1; i <= 20; i++ ) {
		double lambda = dlambda*i;
		ThirringKField bool0( Nt*Ns*Ns, {dimSpinor, Nf, Nf} );
		ThirringWeightFunction weight( bool0, Nt, Ns, dim, Nf, lambda );

		bool0.setRow( conf96, 0 );
		bool0.setRow( conf96, 1 );
		Complex weightVal = weight.calculateWeight();
		std::cerr << lambda << "\t" << std::real(weightVal) << "\t" << std::imag(weightVal) << std::endl;
		testConf[i-1] = weightVal;
	}

	// reference values from Mathematica for Config 96 on lambda = 0.1 .. 2 with delta 0.1
	std::vector<Complex> refConf96 = { 1784.94, 446.235, 198.327, 111.559, 71.3976,	49.5817, 36.4273, 27.8897, 22.0363,	17.8494,
			14.7516, 12.3954, 10.5618, 9.10684, 7.93307, 6.97242, 6.17626, 5.50907,	4.94443, 4.46235};

	for( size_t i = 0; i < testConf.size(); i++ ) {
		if( std::fabs( refConf96[i] - testConf[i] ) > 1e-2 ) {
			std::cerr << "Deviation from reference value found for Config96." << std::endl;
			std::cerr << "Value is " << testConf[i] << " but should be " << refConf96[i] << std::endl;
		}
	}

	// Testing config 1251
//	testConf.clear();
	std::cout << "Testing weight for Configuration 1251" << std::endl << "lambda\tRe(weight)\tIm(weight)" << std::endl;
	for( int i = 1; i <= 20; i++ ) {
		double lambda = dlambda*i;
		ThirringKField bool0( Nt*Ns*Ns, {dimSpinor, Nf, Nf} );
		ThirringWeightFunction weight( bool0, Nt, Ns, dim, Nf, lambda );

		bool0.setRow( row0, 0 );
		bool0.setRow( row1, 1 );
		bool0.setRow( row2, 2 );

		Complex weightVal = weight.calculateWeight();
		std::cerr << lambda << "\t" << std::real(weightVal) << "\t" << std::imag(weightVal) << std::endl;
		testConf[i-1] = weightVal;
	}

	// reference values from Mathematica for Config 1251 on lambda = 0.1 .. 2 with delta 0.1
	std::vector<Complex> refConf1251 = {46330.7, 2895.67, 571.984, 180.979, 74.1291, 35.749, 19.2964,
	11.3112, 7.06153, 4.63307, 3.16445, 2.23431, 1.62217, 1.20603,	0.915175, 0.706951, 0.554719, 0.441346, 0.355512, 0.289567};
	for( size_t i = 0; i < testConf.size(); i++ ) {
		if( std::fabs( refConf1251[i] - testConf[i] ) > 1e-2 ) {
			std::cerr << "Deviation from reference value found for Config1251." << std::endl;
			std::cerr << "Value is " << testConf[i] << " but should be " << refConf1251[i] << std::endl;
		}
	}

	std::cout << std::endl << "Testing update sequence" << std::endl;
	double lambda = dlambda;
	ThirringKField bool0( Nt*Ns*Ns, {dimSpinor, Nf, Nf} );
	ThirringWeightFunction weight( bool0, Nt, Ns, dim, Nf, lambda );

	bool0.setRow( conf96, 0 );
	bool0.setRow( conf96, 1 );

	std::cout << "First weight: " << weight.calculateWeight() << " should be " << refConf96[0] << std::endl;
	weight.keep();
	weight.saveState();

	bool0.setRow( row0, 0 );
	bool0.setRow( row1, 1 );
	bool0.setRow( row2, 2 );

	std::set<size_t> changed = { 0, 1, 2 };

	std::cout << "weight change: " << weight.updateWeight( changed ) << " should be " << refConf1251[0]/refConf96[0] << std::endl;
	weight.keep();
//	std::cout << "new weight: " << weight.calculateWeight() << " should be " << refConf1251[0] << std::endl;

	weight.saveState();
	bool0.setRow( row1, 0 );
	bool0.setRow( row2, 1 );
	bool0.setRow( row0, 2 );

	std::cout << "weight change: " << weight.updateWeight( changed ) << " should be " << 0.25 << std::endl;

	ThirringKField bool1( Nt*Ns*Ns, {dimSpinor, Nf, Nf} );
	ThirringWeightFunction weight1( bool1, Nt, Ns, dim, Nf, lambda );

	bool1.setRow( row1, 0 );
	bool1.setRow( row2, 1 );
	bool1.setRow( row0, 2 );

	std::cout << "First weight: " << weight1.calculateWeight() << " should be " << refConf1251[0]/4. << std::endl;

//	weight.reset();
//	std::cout << "new weight: " << weight.calculateWeight() << " should be " << refConf1251[0] << std::endl;
}
