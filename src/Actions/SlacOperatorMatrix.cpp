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

SlacOperatorMatrix::SlacOperatorMatrix(size_t Nt, size_t Ns, size_t dim) :
		N(Nt*Ns*Ns)
{
	using namespace Eigen;
	if( dim != 3 ) {
		std::cerr << "SlacOperatorMatrix implemented only in 3D." << std::endl;
		exit(1);
	}
	std::vector< Matrix2cd > gamma = cliff.getGammas();
	dimSpinor = 2;
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
	fullSlac = dslac;
}

SlacOperatorMatrix::~SlacOperatorMatrix() {
	// TODO Auto-generated destructor stub
}


const Eigen::MatrixXcd SlacOperatorMatrix::getMatrix() const {
	return dslac;
}

const Complex SlacOperatorMatrix::det() const {
	return dslac.determinant();
}

void SlacOperatorMatrix::deletePoint(size_t x) {
	dslac.row(x) = Eigen::RowVectorXcd::Zero( dslac.cols() );
	dslac.col(x) = Eigen::VectorXcd::Zero( dslac.rows() );
	dslac(x,x) = 1.;
	dslac.row(x+N) = Eigen::RowVectorXcd::Zero( dslac.cols() );
	dslac.col(x+N) = Eigen::VectorXcd::Zero( dslac.rows() );
	dslac(x+N,x+N) = 1.;
}

void SlacOperatorMatrix::addPoint(size_t x) {
}

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

} /* namespace FermiOwn */


