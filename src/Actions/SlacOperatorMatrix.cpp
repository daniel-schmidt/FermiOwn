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
														Nf(numFlavours)
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
	dslac *= I;
	fullSlac = dslac;

	detVal = dslac.determinant();
	fullDet = detVal;
	inverse = dslac.inverse();
	fullInverse = inverse;
}

SlacOperatorMatrix::~SlacOperatorMatrix() {}


const Eigen::MatrixXcd SlacOperatorMatrix::getMatrix() const {
	return dslac;
}

const Complex SlacOperatorMatrix::det() const {
	return detVal;
}

size_t SlacOperatorMatrix::matIndex( size_t x, size_t spin, size_t flavour ) {
	return N*( dimSpinor*flavour + spin ) + x;
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
//	size_t fx = 0;
//	size_t fspin = 0;
//	size_t ff1 = 0;
//	size_t ff2 = 0;
	std::vector< size_t > cols;
	std::vector< size_t > rows;
	bool firstSet = false;
	std::cout << "Index list: " << std::endl;
	for( size_t flavour1 = 0; flavour1 < Nf; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < Nf; flavour2++ ) {
			for( size_t spin = 0; spin < dimSpinor; spin++ ) {
				for( size_t x = 0; x < N; x++ ) {
					if( kxiab.getValue( x, spin, flavour1, flavour2 ) ) {
						std::cout << "(" << matIndex(x, spin, flavour1) << ", " << matIndex( x, spin, flavour2 ) << ")" << std::endl;
						cols.push_back( matIndex(x, spin, flavour1) );
						rows.push_back( matIndex(x, spin, flavour2) );
//						erase( x, spin, flavour1, flavour2 );
//						if( firstSet ) {
////							std::cout << "firstSet" << std::endl;
//							update( matIndex( x, spin, flavour1), matIndex( x, spin, flavour2 ), matIndex( fx, fspin, ff1 ), matIndex( fx, fspin, ff2 ) );
//							firstSet = false;
//						} else {
////							std::cout << "not firstSet" << std::endl;
//							fx = x;
//							fspin = spin;
//							ff1 = flavour1;
//							ff2 = flavour2;
//							firstSet = true;
//						}
					}
				}
			}
		}
	}
	update( rows, cols );
	exit(0);
	if( firstSet ) std::cout << "First still set, thats bad..." << std::endl;
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

	dslac.row( matIndex( x, spin, flavour1 ) ) = Eigen::RowVectorXcd::Zero( dslac.cols() );
	dslac.col( matIndex( x, spin, flavour2 ) ) = Eigen::VectorXcd::Zero( dslac.rows() );

	// set crossing entry correct to get correct full determinant
	dslac(  matIndex( x, spin, flavour1 ), matIndex( x, spin, flavour2 ) ) = 1.;
}
void SlacOperatorMatrix::update( std::vector<size_t> rows, std::vector<size_t> cols ) {
	Eigen::MatrixXcd submat( rows.size(), cols.size() );
	Eigen::MatrixXcd colmat( rows.size(), dslac.cols() );
	Eigen::MatrixXcd rowmat( dslac.rows(), cols.size() );

	for( size_t colIndex = 0; colIndex < cols.size(); colIndex++ ) {
		rowmat.col( colIndex ) = inverse.col( cols[colIndex] );
		for( size_t rowIndex = 0; rowIndex < rows.size(); rowIndex++ ) {
			// fill colmat at first iteraton through the outer loop
			if( colIndex == 0 ) {
				colmat.row( rowIndex ) = inverse.row( rows[rowIndex] );
			}
			submat( rowIndex, colIndex ) = inverse( cols[rowIndex], rows[colIndex] );
		}
	}
	std::cout << submat.inverse() << std::endl;
	detVal *= submat.determinant();
	std::cout << "Det: " << detVal << std::endl;

	std::cout << "rowmat: " << std::endl << rowmat << std::endl;
	std::cout << "colmat: " << std::endl << colmat << std::endl;

	inverse -= rowmat * (submat.inverse() ) * colmat;
	for( size_t index = 0; index < cols.size(); index++ ) {
			inverse( rows[index], cols[index] ) += 1.;
	}
	std::cout << "Inverse: " << std::endl << inverse << std::endl;

}
void SlacOperatorMatrix::update( size_t a, size_t b, size_t m, size_t n ) {

	Complex detCoeff = inverse( n, m ) * inverse( b, a ) - inverse( n, a ) * inverse( b, m );
	if( n==b && m==a ) detCoeff += inverse( b, a ); // Kronecker delta

	std::cout << "detCoeff: " << detCoeff;
	detVal *= detCoeff;
	std::cout << " new detVal: " << detVal << std::endl;
	Eigen::VectorXcd colM = inverse.col( m );
	if( a==m ) colM(b) += 1.;
	Eigen::RowVectorXcd rowN = inverse.row( n );
	if( n==b ) rowN( a ) += 1.;
	Eigen::VectorXcd colA = inverse.col( a );
	Eigen::RowVectorXcd rowB = inverse.row( b );


	Eigen::MatrixXcd updateMat = Eigen::KroneckerProduct<Eigen::VectorXcd, Eigen::RowVectorXcd>( colM, rowN );
	updateMat *= inverse( b, a );

	updateMat -= inverse( n, a ) * Eigen::KroneckerProduct<Eigen::VectorXcd, Eigen::RowVectorXcd>( colM, rowB );
	updateMat -= inverse( b, m ) * Eigen::KroneckerProduct<Eigen::VectorXcd, Eigen::RowVectorXcd>( colA, rowN );

	Complex coeff = inverse( n, m );
	if( n==b && m==a ) coeff += 1.;

	updateMat += coeff * Eigen::KroneckerProduct<Eigen::VectorXcd, Eigen::RowVectorXcd>( colA, rowB );

	updateMat /= detCoeff;

	inverse -= updateMat;
	inverse( n, m ) += 1.;
	inverse( b, a ) += 1.;

}

void SlacOperatorMatrix::setFull() {
	dslac = fullSlac;
	detVal = fullDet;
	inverse = fullInverse;
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

