/*
 * DSlashUpdater_test.cpp
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#include "FieldBoolean.h"
#include "DSlashUpdater.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;
	FieldBoolean bool1( 4, 2, 2, NULL, zeroInit );
	DSlashUpdater dslash( 4, 1, 3, 2 );

	RowVectorXb conf96( 8 );
	conf96 << 0, 1, 0, 0, 0, 0, 1, 0;

	bool1.setRow( conf96, 0 );
	bool1.setRow( conf96, 1 );

	FieldBoolean bool0( 4, 2, 2, NULL, zeroInit );
	FieldBoolean diff = bool1.different( bool0 );

	dslash.calculateUpdateMatrices( bool1, diff );

	Complex detChange = dslash.updateDet();
	Complex currDet = dslash.getDet();
	std::cout << "Det changed by newdet=" << detChange << "*olddet and the factor should be 0.129782. Det is now: " << currDet << std::endl;

	dslash.keep();
	std::cout << "Matrix after update: " << std::endl << dslash.getMatrix() << std::endl << std::endl;

	SparseMat product = (dslash.getMatrix() * dslash.getInverse());
	SparseMat id( product.rows(), product.cols() );
	id.setIdentity();
	std::cout << "Matrix times its inverse is approximately identity: " << product.isApprox( id ) << std::endl;

	FieldBoolean bool2( 4, 2, 2, NULL, zeroInit );
	RowVectorXb conf80( 8 );
	conf80 << 0, 0, 1, 0, 0, 1, 0, 0;
	bool2.setRow( conf80, 0 );
	bool2.setRow( conf80, 1 );

	diff = bool2.different( bool1 );

	dslash.calculateUpdateMatrices( bool2, diff );
	detChange = dslash.updateDet();
	currDet = dslash.getDet();
	std::cout << "Det changed by newdet=" << detChange << "*olddet and the factor should be 1. Det is now: " << currDet << std::endl;

	dslash.keep();
	std::cout << "Matrix after update: " << std::endl << dslash.getMatrix() << std::endl << std::endl;

	product = (dslash.getMatrix() * dslash.getInverse());
	id = SparseMat( product.rows(), product.cols() );
	id.setIdentity();
	std::cout << "Matrix times its inverse is approximately identity: " << product.isApprox( id ) << std::endl;
	std::cout << "Non-Zeros: " << dslash.getInverse().nonZeros() << std::endl;
}
