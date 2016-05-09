/*
 * SlacOperatorMatrix.cpp
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#include "SlacOperatorMatrix.h"

namespace FermiOwn {

SlacOperatorMatrix::SlacOperatorMatrix( size_t size ) :
														N(size),
														dslac(make1D(size))
{}

SlacOperatorMatrix::SlacOperatorMatrix( size_t Nt, size_t Ns, size_t dim, size_t numFlavours ) :
														N(Nt*Ns*Ns),
														dimSpinor(2),
														Nf(numFlavours),
														rowMap(N*numFlavours*dimSpinor),
														colMap(N*numFlavours*dimSpinor)
{
	using namespace Eigen;
	if( dim != 3 ) {
		std::cerr << "SlacOperatorMatrix implemented only in 3D." << std::endl;
		exit(1);
	}
	std::vector< Matrix2cd > gamma = cliff.getGammas();
	std::vector< MatrixXcd > Gammas;
	for( auto gammaMu : gamma) {
		MatrixXcd GammaMu = KroneckerProduct<Matrix2cd, MatrixXcd>( gammaMu, MatrixXcd::Identity(N,N) );
		Gammas.push_back(GammaMu);
	}

	MatrixXcd dx = KroneckerProduct<MatrixXd, MatrixXcd>( MatrixXd::Identity(dimSpinor*Ns*Ns,dimSpinor*Ns*Ns), make1D(Nt) );
	MatrixXcd dy( KroneckerProduct<MatrixXcd, MatrixXd>( make1D(Ns), MatrixXd::Identity(Nt,Nt) ) );
	dy = KroneckerProduct<MatrixXd, MatrixXcd>( MatrixXd::Identity(dimSpinor*Ns,dimSpinor*Ns), dy ).eval();
	MatrixXcd dz( KroneckerProduct<MatrixXcd, MatrixXd>( make1D(Ns), MatrixXd::Identity(Ns*Nt, Ns*Nt) ) );
	dz = KroneckerProduct<MatrixXd, MatrixXcd>( MatrixXd::Identity(dimSpinor,dimSpinor), dz ).eval();
	dslac = Gammas[0]*dx+Gammas[1]*dy+Gammas[2]*dz;
	dslac = KroneckerProduct<MatrixXd, MatrixXcd>( MatrixXd::Identity(Nf,Nf), dslac ).eval();
	fullSlac = dslac;
	resetMaps();
}

SlacOperatorMatrix::~SlacOperatorMatrix() {}


const Eigen::MatrixXcd SlacOperatorMatrix::getMatrix() const {
	return dslac;
}

const Complex SlacOperatorMatrix::det() const {
	return (dslac*I).determinant();
}

//void SlacOperatorMatrix::deletePoint(size_t x) {
//	dslac.row(x) = Eigen::RowVectorXcd::Zero( dslac.cols() );
//	dslac.col(x) = Eigen::VectorXcd::Zero( dslac.rows() );
//	dslac(x,x) = 1.;
//	dslac.row(x+N) = Eigen::RowVectorXcd::Zero( dslac.cols() );
//	dslac.col(x+N) = Eigen::VectorXcd::Zero( dslac.rows() );
//	dslac(x+N,x+N) = 1.;
//}

void SlacOperatorMatrix::erase( const FieldBoolean& kxiab ) {
	for( size_t flavour1 = 0; flavour1 < Nf; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < Nf; flavour2++ ) {
			for( size_t spin = 0; spin < dimSpinor; spin++ ) {
				for( size_t x = 0; x < N; x++ ) {
					if( kxiab.getValue( x, spin, flavour1, flavour2 ) ) {
						erase( x, spin, flavour1, flavour2 );
					}
				}
			}
		}
	}
}

void SlacOperatorMatrix::erase( size_t x, size_t spin, size_t flavour1, size_t flavour2 ) {
	// check for correct input
	if( x >= N ) {
		std::cerr << "Trying to erase col/row at x=" << x << ", which is larger or equal than the lattice volume!" << std::endl;
		exit(1);
	} else if ( spin >= dimSpinor ) {
		std::cerr << "Trying to erase col/row at spin=" << spin << ", which is larger or equal than the spinor size!" << std::endl;
		exit(1);
	} else if ( flavour1 >= Nf || flavour2 >= Nf ) {
		std::cerr << "Trying to erase col/row at flavour1=" << flavour1 << " or flavour2=" << flavour2 <<  ", which is larger or equal than Nf=" << Nf << std::endl;
		exit(1);
	}
	size_t rowIndex = N*( dimSpinor*flavour1 + spin ) + x;
	size_t colIndex = N*( dimSpinor*flavour2 + spin ) + x;

//	if( rowMap( rowIndex ) == -1 ) {
//		std::cerr << "Row already deleted!" << std::endl;
////		exit(1);
//	} else if ( colMap( colIndex ) == -1 ) {
//		std::cerr << "Col already deleted!" << std::endl;
////		exit(1);
//	}

	dslac.row(rowIndex) = Eigen::RowVectorXcd::Zero( dslac.cols() );
	dslac.col(colIndex) = Eigen::VectorXcd::Zero( dslac.rows() );

	// set crossing entry correct to get correct full determinant
//	if( ( rowMap( rowIndex ) + colMap( colIndex ) ) % 2 == 0 )
		dslac(rowIndex,colIndex) = 1.;
//	else
//		dslac(rowIndex,colIndex) = -1.;

	// update row and col maps
//	rowMap( rowIndex ) = -1;
//	int tailSize = rowMap.size() - rowIndex;
//	rowMap.tail( tailSize ) -= Eigen::VectorXi::Ones(tailSize);
//	colMap( colIndex ) = -1;
//	tailSize = colMap.size() - colIndex;
//	colMap.tail( tailSize ) -= Eigen::VectorXi::Ones(tailSize);
}

void SlacOperatorMatrix::setFull() {
	dslac = fullSlac;
	resetMaps();
}

void SlacOperatorMatrix::resetMaps() {
	for( int i = 0; i < rowMap.size(); i++ ) {
		rowMap(i) = i;
		colMap(i) = i;
	}
}

//void SlacOperatorMatrix::addPoint(size_t x) {
//	dslac.row(x) = fullSlac.row(x);
//	dslac.row(x+N) = fullSlac.row(x+N);
//	dslac.col(x) = fullSlac.col(x);
//	dslac.col(x+N) = fullSlac.col(x+N);
//}

Eigen::MatrixXcd SlacOperatorMatrix::make1D( size_t size) {
	Eigen::MatrixXcd ret( size, size );
	for( size_t i = 0; i < size; i++ ) {
		for( size_t j = 0; j < size; j++ ) {
			if( i == j ) {
				ret(i,j) = 0.;
			} else {
				int t = i-j;
				ret(i,j) = PI*pow(-1, t) / ( size * sin(PI*t/size) );
			}
		}
	}
	return ret;
}

//void SlacOperatorMatrix::eraseCol( size_t x, size_t spin ) {
//	// check for correct input
//	if( x >= N ) {
//		std::cerr << "Trying to erase a column at x=" << x << ", which is larger than the lattice volume!" << std::endl;
//		exit(1);
//	} else if ( spin >= dimSpinor ) {
//		std::cerr << "Trying to erase a column at spin=" << spin << ", which is larger than the spinor size!" << std::endl;
//		exit(1);
//	}
//}
//
//void SlacOperatorMatrix::eraseCols( VectorXb kx0, VectorXb kx1 ) {
//	assert( kx0.size() == N && kx1.size() == N );
//
//	int cols = kx0.size() - kx0.count() + kx1.size() - kx1.count(); // new number of columns
//	Eigen::MatrixXcd newSlac( dslac.rows(), cols );
//
//	int colIndex = 0;
//
//	for( size_t x = 0; x < N; x++ ) {
//		if( !kx0(x) ) {
//			newSlac.col( colIndex ) = dslac.col( x );
//			colIndex++;
//		}
//	}
//	for( size_t x = 0; x < N; x++ ) {
//		if( !kx1(x) ) {
//			newSlac.col( colIndex ) = dslac.col( x + N );
//			colIndex++;
//		}
//	}
//	dslac = newSlac;
//}
//
//void SlacOperatorMatrix::eraseRows( VectorXb kx0, VectorXb kx1 ) {
//	assert( kx0.size() == N && kx1.size() == N );
//
//	int rows = kx0.size() - kx0.count() + kx1.size() - kx1.count(); // new number of rows
//	Eigen::MatrixXcd newSlac( rows, dslac.cols() );
//
//	int rowIndex = 0;
//
//	// first spinor component
//	for( size_t x = 0; x < N; x++ ) {
//		if( !kx0(x) ) {
//			newSlac.row( rowIndex ) = dslac.row( x );
//			rowIndex++;
//		}
//	}
//	// second spinor component
//	for( size_t x = 0; x < N; x++ ) {
//		if( !kx1(x) ) {
//			newSlac.row( rowIndex ) = dslac.row( x + N );
//			rowIndex++;
//		}
//	}
//	dslac = newSlac;
//}

} /* namespace FermiOwn */

