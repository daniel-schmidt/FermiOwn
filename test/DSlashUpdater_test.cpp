/*
 * DSlashUpdater_test.cpp
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#include "FieldBoolean.h"
#include "DSlashUpdater.h"
#include "../src/Actions/ConfigPerPointGenerator.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;
	FieldBoolean bool1( 4, 2, 2, NULL, zeroInit );
	DSlashUpdater dslash( 4, 1, 3, 2 );
	DSlashUpdater initialDSlash( 4, 1, 3, 2 );

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
	std::cout << "Inverse after update: " << std::endl << dslash.getInverse() << std::endl << std::endl;

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
	std::cout << "Non-Zeros of inverse: " << dslash.getInverse().nonZeros() << " of " << dslash.getInverse().size() << " entries."<< std::endl;

	std::cout << std::endl << std::endl << "Further tests with generated configs" << std::endl << "========================================="<< std::endl;

	ConfigPerPointGenerator confGen( 2, 2 );
	confGen.generateAllowedConfs();
	MatrixXb allowedConfs = confGen.getAllConfs();

	FieldBoolean fbool( 4, 2, 2, NULL, zeroInit );

	dslash = initialDSlash;
	DSlashUpdater alwaysReset( 4, 1, 3, 2 );
	FieldBoolean initialField = fbool;

	for( int row = 0; row < allowedConfs.rows(); row++ ) {
		std::cout << "Testing configuration " << std::endl;
		FieldBoolean oldField = fbool;
		fbool.setRow( allowedConfs.row( row ), 0 );
		fbool.setRow( allowedConfs.row( row ), 1 );
		fbool.Print();

		std::cout << "Diff to last accepted: " << std::endl;
		fbool.different( oldField ).Print();

		dslash.calculateUpdateMatrices( fbool, fbool.different( oldField) );
		alwaysReset.calculateUpdateMatrices( fbool, fbool.different( initialField ) );

		std::cout << dslash.getDet() << " " << alwaysReset.getDet() << std::endl;
		if( std::fabs( dslash.getDet() ) > ZERO_TOL ) {
			dslash.keep();
			std::cout << "keep" << std::endl;
		} else {
			dslash.reset();
			std::cout << "reset" << std::endl;
			fbool = oldField;
		}
		alwaysReset.reset();
		product = (dslash.getMatrix() * dslash.getInverse());
		id = SparseMat( product.rows(), product.cols() );
		id.setIdentity();
		std::cout << "Matrix times its inverse is approximately identity: " << product.isApprox( id ) << std::endl << std::endl;
	}

	//	for( int i = 0; i < 3; i++ ) {
	//		std::cout << std::endl << "Testing algorithm " <<  i << std::endl;
	//		dslacNf2.setFull();
	//		fbool = FieldBoolean( 4, 2, 2, NULL, zeroInit );
	//		for( int row = 0; row < allowedConfs.rows(); row++ ) {
	//			FieldBoolean oldField = fbool;
	//			//			SlacOperatorMatrix oldSlac = dslacNf2;
	//			dslacNf2.saveState();
	//			std::cout << "Det before: " << dslacNf2.det() << std::endl;
	//			fbool.Print();
	//			fbool.setRow( allowedConfs.row( row ), 0 );
	//			fbool.setRow( allowedConfs.row( row ), 1 );
	//			dslacNf2.update( fbool, fbool.different( oldField ), static_cast<SlacOperatorMatrix::updateType>( i ) );
	//			std::cout << "Det after: " << dslacNf2.det() << std::endl;
	//			fbool.Print();
	//			if( abs(dslacNf2.det()) < 1e-10 ) {
	//				fbool = oldField;
	//				dslacNf2.resetState();
	//				//				dslacNf2 = oldSlac;
	//			}
	//			std::cout << std::endl;
	//		}
	//	}

	std::cout << std::endl << std::endl << "Testing fixed update sequence" << std::endl << "========================================="<< std::endl;;

	fbool = FieldBoolean( 4, 2, 2, NULL, zeroInit );
	FieldBoolean fboolInitial = fbool;
	fbool.setValue( 1, 0, 0, 0, 1 );
	fbool.setValue( 1, 1, 1, 0, 1 );
	fbool.setValue( 1, 0, 0, 1, 0 );
	fbool.setValue( 1, 1, 1, 1, 0 );

	dslash = initialDSlash;
	std::cout << "Det before: " << dslash.getDet() << std::endl;
	diff = fbool.different( fboolInitial );
	fbool.Print();
	//		diff.Print();
	dslash.calculateUpdateMatrices( fbool, diff );
	std::cout << "Det after first delete: " << dslash.getDet() << std::endl;
	dslash.keep();

	fboolInitial = fbool;
	fbool.invert( 0, 1, 0, 0 );
	fbool.invert( 0, 1, 1, 1 );
	fbool.invert( 1, 0, 0, 0 );
	fbool.invert( 1, 0, 1, 1 );
	fbool.invert( 0, 0, 0, 1 );
	fbool.invert( 0, 0, 1, 0 );
	fbool.invert( 1, 1, 0, 1 );
	fbool.invert( 1, 1, 1, 0 );
	fbool.Print();
	diff = fbool.different( fboolInitial );

	dslash.calculateUpdateMatrices( fbool, diff );
	dslash.updateDet();
	std::cout << "switched diag->off, det: " << dslash.getDet() << std::endl;
	dslash.keep();

	fboolInitial = fbool;
	fbool.invert( 0, 0, 0, 1 );
	fbool.invert( 0, 0, 1, 0 );
	fbool.invert( 1, 1, 0, 1 );
	fbool.invert( 1, 1, 1, 0 );
	fbool.Print();

	diff = fbool.different( fboolInitial );
	dslash.calculateUpdateMatrices( fbool, diff );
	dslash.updateDet();
	std::cout << "deleted more, det: " << dslash.getDet() << std::endl;
	dslash.keep();

	fboolInitial = fbool;
	fbool.invert( 0, 1, 0, 0 );
	fbool.invert( 0, 1, 1, 1 );
	fbool.invert( 1, 0, 0, 0 );
	fbool.invert( 1, 0, 1, 1 );
	fbool.Print();
	diff = fbool.different( fboolInitial );
	dslash.calculateUpdateMatrices( fbool, diff );
	std::cout << "Re-adding, det again: " << dslash.getDet() << std::endl;
	dslash.keep();

	fboolInitial = fbool;
	fbool.invert( 0, 0, 0, 1 );
	fbool.invert( 0, 0, 1, 0 );
	fbool.invert( 1, 1, 0, 1 );
	fbool.invert( 1, 1, 1, 0 );
	fbool.Print();
	diff = fbool.different( fboolInitial );
	dslash.calculateUpdateMatrices( fbool, diff );
	std::cout << "Full det again: " << dslash.getDet() << std::endl;
	dslash.keep();
}
