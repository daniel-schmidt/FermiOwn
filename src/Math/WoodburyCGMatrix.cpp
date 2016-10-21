/*
 * WoodburyMatrix.cpp
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#include "WoodburyCGMatrix.h"

namespace FermiOwn {

WoodburyCGMatrix::WoodburyCGMatrix() :
	size( 0 ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true ),
	deleteOnly( false )
{};

WoodburyCGMatrix::WoodburyCGMatrix( const size_t matSize ) :
	mat( matSize, matSize ),
	size( matSize ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true ),
	deleteOnly( false )
{
}

WoodburyCGMatrix::WoodburyCGMatrix( const size_t matSize, const MatCoeffList & coeffList, const bool selfAdjointInit ) :
	mat( matSize, matSize ),
	size( matSize ),
	matNeedsUpdate( true ),
	detNeedsUpdate( true ),
	invNeedsUpdate( true ),
	deleteOnly( false )
{
	setFromCoeffList( coeffList, selfAdjointInit );
}

//WoodburyMatrix::WoodburyMatrix( const WoodburyMatrix& other ) :
//		mat( other.mat ),
//		size( other.size),
//		inv( other.inv ),
//		det( other.det ),
//		matNeedsUpdate( false ),
//		detNeedsUpdate( false ),
//		invNeedsUpdate( false ),
//		deleteOnly( false )
//{}


WoodburyCGMatrix::~WoodburyCGMatrix() {
}

//void swap( WoodburyMatrix& first, WoodburyMatrix& second ) {
//	using std::swap;
//	swap( first.mat, second.mat );
//	swap( first.size, second.size );
//	swap( first.inv, second.inv );
//	swap( first.det, second.det );
//	first.matNeedsUpdate = false;
//	first.detNeedsUpdate = false;
//	first.invNeedsUpdate = false;
//	second.matNeedsUpdate = false;
//	second.detNeedsUpdate = false;
//	second.invNeedsUpdate = false;
//}

void WoodburyCGMatrix::setFromCoeffList( const MatCoeffList & coeffList, const bool selfAdjointInit ) {
	mat.setFromTriplets( coeffList.begin(), coeffList.end() );

	if( selfAdjointInit ) {
		SparseMat tmp = mat.selfadjointView<Eigen::Upper>();
		mat = tmp;
	}

	matNeedsUpdate = false;

	Eigen::SparseLU< SparseMat > solver;	//TODO: check, if we can use an SPD solver here for SLAC
	solver.compute( mat );
	if( solver.info() != Eigen::Success ) {
		std::cerr << "WoodburyMatrix cannot calculate matrix decomposition to obtain determinant." << std::endl;
		exit(1);
	}

	det =  solver.determinant();
	detNeedsUpdate = false;
}

void WoodburyCGMatrix::setUpdateMatrices( const SparseMat* colMatrix, const SparseMat* rowMatrix ) {
	if( colMatrix->rows() != int( size ) ) {
		std::cerr << "WoodburyMatrix got a colMatrix of wrong size, not matching the matrix dimension " << size << ", instead size " << colMatrix->rows() << std::endl;
		exit(1);
	}
	if( rowMatrix->cols() != int( size ) ) {
		std::cerr << "WoodburyMatrix got a rowMatrix of wrong size, not matching the matrix dimension " << size << ", instead size " << rowMatrix->cols() << std::endl;
		exit(1);
	}
	if( colMatrix->cols() != rowMatrix->rows() ) {
		std::cerr << "The two update matrices passed to WoodburyMatrix do not fit in size! They are col: " << colMatrix->cols() << " and row: " << rowMatrix->rows() << std::endl;
		exit(1);
	}
	if( colMatrix->nonZeros() == 0 || rowMatrix->nonZeros() == 0 ) {
		std::cerr << "Error: WoodburyMatrix got update matrix with no non-zero element!" << std::endl;
		exit(1);
	}
	U = colMatrix;
	V = rowMatrix;

	detNeedsUpdate = true;
	invNeedsUpdate = true;
	matNeedsUpdate = true;
	deleteOnly = false;
}

void WoodburyCGMatrix::setUpdateMatricesDeleteOnly( const SparseMat* colMatrix, const SparseMat* rowMatrix, const std::vector<size_t>& colsToDelete, const std::vector<size_t>& rowsToDelete ) {
	setUpdateMatrices( colMatrix, rowMatrix );
	rows = rowsToDelete;
	cols = colsToDelete;
	deleteOnly = true;
}

void WoodburyCGMatrix::update() {
	updateDet();
	updateMatrix();
	updateInverse();
}

void WoodburyCGMatrix::updateMatrix() {
	if( matNeedsUpdate ) {
		mat += (*U) * (*V);	//TODO: setting entries directly is most likely more efficient
		mat.prune( Complex( 0., 0.) ); //could be used to delete entries from matrix that got zero but are still listed as elements, but does a copy of the matrix...
		matNeedsUpdate = false;
	}
}

Complex WoodburyCGMatrix::updateDet() {
	Complex detChange = 1.;

	if( detNeedsUpdate ) {
		Eigen::BiCGSTAB<SparseMat> solver;
		solver.compute( mat );

		submat.resize( rows.size(), cols.size() );
		submat.setIdentity();
//		Eigen::MatrixXcd AU( U->rows(), U->cols() );
//		for( int i = 0; i < U->cols(); i++ ) {
//			Eigen::VectorXcd b = U->col(i);
//			AU.col(i) = solver.solve( b );
//		}
		Eigen::MatrixXcd AU = solver.solve( Eigen::MatrixXcd(*U) );

		submat += (*V) * AU;
		detChange = submat.determinant();
//		std::cout << "AU: " << std::endl << AU << std::endl << "submat:" << std::endl << submat << std::endl << "detChange: " << detChange << std::endl;
		det *= detChange;
		detNeedsUpdate = false;
	}
	return detChange;
}

void WoodburyCGMatrix::updateInverse() {
}

} /* namespace FermiOwn */
