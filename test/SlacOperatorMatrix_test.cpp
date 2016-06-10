/*
 * SlacOperatorMatrix_test.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include <iostream>
#include <Eigen/Dense>
#include "FieldBoolean.h"
#include "ConfigGenerator.h"
#include "SlacOperatorMatrix.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;

	if( argc != 2 ) {
		std::cerr << "SlacOperatorMatrix test needs 1 argument: matrix size" << std::endl;
		exit(1);
	}

	size_t N = atoi( argv[1]);
	if( N > 5 || N < 2 ) {
		std::cerr << "No reference for matrix size " << N << " implemented." << std::endl;
		exit(1);
	}

	SlacOperatorMatrix dslac( N );

	Eigen::Matrix2cd sol2;
	sol2 << 0., 1.5708, -1.5708, 0.;

	Eigen::Matrix3cd sol3;
	sol3 << 0., 1.2092, -1.2092, -1.2092, 0., 1.2092, 1.2092, -1.2092, 0.;

	Eigen::Matrix4cd sol4;
	sol4 << 0., 1.11072, -0.785398, 1.11072,
			-1.11072, 0., 1.11072, -0.785398,
			0.785398, -1.11072, 0., 1.11072,
			-1.11072, 0.785398, -1.11072, 0.;

	Eigen::MatrixXcd sol5(5,5);
	sol5 << 0., 1.06896, -0.660653, 0.660653, -1.06896,
			-1.06896, 0., 1.06896, -0.660653, 0.660653,
			0.660653, -1.06896, 0., 1.06896, -0.660653,
			-0.660653, 0.660653, -1.06896, 0., 1.06896,
			1.06896, -0.660653, 0.660653, -1.06896, 0.;

	Eigen::MatrixXcd sol(N,N);
	switch ( N ) {
	case 2: sol = sol2; break;
	case 3: sol = sol3; break;
	case 4: sol = sol4; break;
	case 5: sol = sol5; break;
	}

	std::cout << dslac.getMatrix() << " should be " << sol << std::endl;
	if(dslac.getMatrix().isApprox(sol, 1e-5) ) std::cout << "Size " << N << " works!" << std::endl;

	std::cout << std::endl << "Testing 3d initialization:" << std::endl;
	size_t Nf = 1;
	size_t dim = 3;
	SlacOperatorMatrix dslac3d( N, N+1, dim, Nf );
	std::cout << "Full determinant:" << dslac3d.det() << std::endl;
	//	dslac.setFull();
	//	for( size_t x = 0; x < N*(N+1)*(N+1); x++ ) {
	//		dslac3d.erase( x, 0, 0, 0);
	//		dslac3d.erase( x, 1, 0, 0);
	//	}
	//	std::cout << "Deleted everything:" << dslac3d.getMatrix().isIdentity() << std::endl;
	//	std::cout << "Deleted everything, det=" << dslac3d.det() << std::endl;
	//
	//	SlacOperatorMatrix dslacNf1( N, 1, dim, 1 );
	//	std::cout << "Single Flavour: " << std::endl << dslacNf1.getMatrix() << std::endl;
	SlacOperatorMatrix dslacNf2( 4, 1, dim, 2 );
	std::cout << "Two Flavour full matrix: " << std::endl << dslacNf2.getMatrix() << std::endl;
	std::cout << "Two Flavour determinant: " << std::endl << dslacNf2.det() << std::endl;
	//	dslacNf2.erase( 0, 0, 0, 1 );
	//	dslacNf2.erase( 1, 1, 0, 1 );
	//	dslacNf2.erase( 0, 0, 1, 0 );
	//	dslacNf2.erase( 1, 1, 1, 0 );
	//	std::cout << "Two Flavour deleted: " << std::endl << dslacNf2.getMatrix() << std::endl;
	//	std::cout << "Two Flavour deleted determinant: " << std::endl << dslacNf2.det() << std::endl;
	//	dslacNf2.setFull();

	ConfigGenerator confGen( 2, 2 );
	confGen.generateAllowedConfs();
	MatrixXb allowedConfs = confGen.getAllConfs();

	FieldBoolean fbool( 4, 2, 2, NULL, zeroInit );

	for( int row = 0; row < allowedConfs.rows(); row++ ) {
		std::cout << "Testing configuration " << std::endl;
		fbool.Print();
		dslacNf2.setFull();
		for( int i = 0; i < 3; i++ ) {
			FieldBoolean oldField = fbool;
			fbool.setRow( allowedConfs.row( row ), 0 );
			fbool.setRow( allowedConfs.row( row ), 1 );
			dslacNf2.update( fbool, fbool.different( oldField), static_cast<SlacOperatorMatrix::updateType>( i ) );
			std::cout << dslacNf2.det() << " ";
		}

		std::cout << std::endl;
	}

	for( int i = 0; i < 3; i++ ) {
		std::cout << std::endl << "Testing algorithm " <<  i << std::endl;
		dslacNf2.setFull();
		for( int row = 0; row < allowedConfs.rows(); row++ ) {
			FieldBoolean oldField = fbool;
			SlacOperatorMatrix oldSlac = dslacNf2;
			fbool.setRow( allowedConfs.row( row ), 0 );
			fbool.setRow( allowedConfs.row( row ), 1 );
			dslacNf2.update( fbool, fbool.different( oldField ), static_cast<SlacOperatorMatrix::updateType>( i ) );
			//			fbool.Print();
			std::cout << dslacNf2.det() << std::endl;
			if( abs(dslacNf2.det()) < 1e-10 ) {
				fbool = oldField;
				dslacNf2 = oldSlac;
			}
		}
	}

	std::cout << "Testing fixed update sequence" << std::endl;

	for( int i = 0; i < 3; i++ ) {
		std::cout << std::endl << "Testing algorithm " <<  i << std::endl;
		SlacOperatorMatrix::updateType upType = static_cast<SlacOperatorMatrix::updateType>( i );
		fbool = FieldBoolean( 4, 2, 2, NULL, zeroInit );
		FieldBoolean fboolInitial = fbool;
		fbool.setValue( 1, 0, 0, 0, 1 );
		fbool.setValue( 1, 1, 1, 0, 1 );
		fbool.setValue( 1, 0, 0, 1, 0 );
		fbool.setValue( 1, 1, 1, 1, 0 );
		fbool.Print();

		dslacNf2.setFull();
		std::cout << "Det before: " << dslacNf2.det() << std::endl;
		FieldBoolean diff = fbool.different( fboolInitial );
//		diff.Print();
		dslacNf2.update( fbool, diff, upType );
		std::cout << "Det after first delete: " << dslacNf2.det() << std::endl;

		fboolInitial = fbool;
		fbool.invert( 0, 1, 0, 0 );
		fbool.invert( 0, 1, 1, 1 );
		fbool.invert( 1, 0, 0, 0 );
		fbool.invert( 1, 0, 1, 1 );
		fbool.invert( 0, 0, 0, 1 );
		fbool.invert( 0, 0, 1, 0 );
		fbool.invert( 1, 1, 0, 1 );
		fbool.invert( 1, 1, 1, 0 );
//		fbool.Print();
		diff = fbool.different( fboolInitial );
		diff.Print();

		dslacNf2.update( fbool, diff, upType );
		std::cout << "switched diag->off, det: " << dslacNf2.det() << std::endl;

		fboolInitial = fbool;
		fbool.invert( 0, 0, 0, 1 );
		fbool.invert( 0, 0, 1, 0 );
		fbool.invert( 1, 1, 0, 1 );
		fbool.invert( 1, 1, 1, 0 );

		diff = fbool.different( fboolInitial );
		dslacNf2.update( fbool, diff, upType );
		diff.Print();
		std::cout << "deleted more, det: " << dslacNf2.det() << std::endl;

		fboolInitial = fbool;
		fbool.invert( 0, 1, 0, 0 );
		fbool.invert( 0, 1, 1, 1 );
		fbool.invert( 1, 0, 0, 0 );
		fbool.invert( 1, 0, 1, 1 );
		fbool.invert( 0, 0, 0, 1 );
		fbool.invert( 0, 0, 1, 0 );
		fbool.invert( 1, 1, 0, 1 );
		fbool.invert( 1, 1, 1, 0 );

		diff = fbool.different( fboolInitial );
		dslacNf2.update( fbool, diff, upType );
		diff.Print();
		std::cout << "Full det again: " << dslacNf2.det() << std::endl;
	}

	//		std::vector< Complex > dets;
	//	for( int i = 0; i < 3; i++ ) {
	//		std::cout << std::endl << "Testing algorithm " <<  i << std::endl;
	//		dslacNf2.update( fbool, diff, static_cast<SlacOperatorMatrix::updateType>( 0  ) );
	//		std::cout << "Two Flavour deleted: " << std::endl << dslacNf2.getMatrix() << std::endl;
	//		std::cout << "Two Flavour deleted determinant: " << std::endl << dslacNf2.det() << std::endl;
	//		dets.push_back( dslacNf2.det() );
	//		dslacNf2.setFull();
	//	}
	//
	//	if( abs( dets[0] - dets[1] ) > 1e-10 ) std::cout << "Difference between erase and separate." << std::endl;
	//	if( abs( dets[0] - dets[2] ) > 1e-10 ) std::cout << "Difference between erase and combined." << std::endl;
	//
	//	std::cout << "Testing simultaneous add and delete" << std::endl;
	//	fboolInitial = fbool;
	//
	//	for( size_t x = 0; x < 2*3*3; x++ ) {
	//		for( size_t spin = 0; spin < 2; spin++ ) {
	//			for( size_t flavour1 = 0; flavour1 <=1; flavour1++ ) {
	//				for( size_t flavour2 = 0; flavour2 <= 1; flavour2 ++ ) {
	//					fbool.invert( x, spin, flavour1, flavour2 );
	//					fbool.enforceConstraint(x, spin, flavour1, flavour2 );
	//					fbool.Print();
	//					diff = fbool.different( fboolInitial );
	//
	//
	//					for( int i = 0; i < 3; i++ ) {
	//						std::cout << std::endl << "Testing algorithm " <<  i << std::endl;
	//						dslacNf2.update( fbool, diff, static_cast<SlacOperatorMatrix::updateType>( 0  ) );
	////						std::cout << "Two Flavour deleted: " << std::endl << dslacNf2.getMatrix() << std::endl;
	//						std::cout << "Two Flavour deleted determinant: " << std::endl << dslacNf2.det() << std::endl;
	//						dets.push_back( dslacNf2.det() );
	//						dslacNf2.setFull();
	//					}
	//
	//					if( abs( dets[0] - dets[1] ) > 1e-10 ) std::cout << "Difference between erase and separate." << std::endl;
	//					if( abs( dets[0] - dets[2] ) > 1e-10 ) std::cout << "Difference between erase and combined." << std::endl;
	//				}
	//			}
	//		}
	//	}
}
