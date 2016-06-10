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

const Eigen::MatrixXcd SlacOperatorMatrix::getInverse() const {
	return inverse;
}

const Complex SlacOperatorMatrix::det() const {
	return detVal;
}

size_t SlacOperatorMatrix::matIndex( size_t x, size_t spin, size_t flavour ) {
	return N*( dimSpinor*flavour + spin ) + x;
}

void SlacOperatorMatrix::update( FieldBoolean kxiab, FieldBoolean changed, updateType upType ) {

	// check what changed, and if it requires adding or deleting row/col pairs
	std::vector<size_t> rowsAdd;
	std::vector<size_t> rowsDel;
	std::vector<size_t> colsAdd;
	std::vector<size_t> colsDel;

	for( size_t flavour1 = 0; flavour1 < Nf; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < Nf; flavour2++ ) {
			for( size_t spin = 0; spin < dimSpinor; spin++ ) {
				for( size_t x = 0; x < N; x++ ) {
					if( changed.getValue( x, spin, flavour1, flavour2 ) ) {
						if( kxiab.getValue( x, spin, flavour1, flavour2 ) ) {
							colsDel.push_back( matIndex(x, spin, flavour1) );
							rowsDel.push_back( matIndex(x, spin, flavour2) );
						} else {
							colsAdd.push_back( matIndex(x, spin, flavour1) );
							rowsAdd.push_back( matIndex(x, spin, flavour2) );
						}
					}
				}
			}
		}
	}

	std::cout << "Delete: " << std::endl;
	for( size_t col : colsDel) std::cout << col << " ";
	std::cout << std::endl;
	for( size_t row : rowsDel ) std::cout << row << " ";
	std::cout << std::endl << "Add: " << std::endl;
	for( size_t col : colsAdd ) std::cout << col << " ";
	std::cout << std::endl;
	for( size_t row : rowsAdd ) std::cout << row << " ";
	std::cout << std::endl;

	switch( upType ) {
	case eraseUpdate:
	{
		setFull();
		erase( kxiab );
		break;
	}
	case separateUpdate:
		deleteEntries( rowsDel, colsDel );
		if( abs( detVal ) > 10e-10 )
			addEntries( rowsAdd, colsAdd );
		else
			std::cerr << "Warning: separate Update of deleting and adding failed to add..." << std::endl;
		break;
	case combinedUpdate:
		combined( rowsAdd, colsAdd, rowsDel, colsDel );
		break;
	}
}

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

	Eigen::FullPivLU< Eigen::MatrixXcd > LUdecomp( dslac );
	if( LUdecomp.isInvertible() ) {
		inverse = LUdecomp.inverse();
		detVal = LUdecomp.determinant();
	} else {
		detVal = 0.;
	}
}




/* =======================================================================================
 * Private methods
 * =======================================================================================
 */


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

void SlacOperatorMatrix::combined( std::vector<size_t> addRows, std::vector<size_t> addCols, std::vector<size_t> delRows, std::vector<size_t> delCols) {
	if( addRows.size() != addCols.size() || delRows.size() != delCols.size() ) {
		std::cerr << "Trying to add not the same number of rows/cols. This should be impossible!" << std::endl;
		exit(1);
	}

	// number of rows/cols to modify
	size_t updateRank = addCols.size() + delCols.size();

//	if( addRows.size() != 0 && delRows.size() != 0 ) {
		Eigen::MatrixXcd subMat( updateRank, updateRank );
		for( size_t colIndex = 0; colIndex < updateRank; colIndex++ ) {
			for( size_t rowIndex = 0; rowIndex < updateRank; rowIndex++ ) {
				if( colIndex < addRows.size() && rowIndex < addRows.size() ) {
					subMat( rowIndex, colIndex ) = -fullSlac( addRows[rowIndex], addCols[colIndex] )-dslac( addRows[rowIndex], addCols[colIndex] );
				} else if( colIndex >= addRows.size() && rowIndex >= addRows.size() ) {
					subMat( rowIndex, colIndex ) = dslac( delRows[ rowIndex-addRows.size() ], delCols[ colIndex-addRows.size() ] );
					if( colIndex == rowIndex ) subMat( rowIndex, colIndex ) += 1.;
				} else {
					subMat( rowIndex, colIndex ) = 0.;
				}
			}
		}

		for( size_t index = 0; index < addCols.size(); index++ ) {
			deletedCols.erase( std::remove( deletedCols.begin(), deletedCols.end(), addCols[index] ), deletedCols.end() );
			deletedRows.erase( std::remove( deletedRows.begin(), deletedRows.end(), addRows[index] ), deletedRows.end() );
		}

		for( size_t index = 0; index < delCols.size(); index++ ) {
			deletedCols.push_back( delCols[index] );		//TODO: check return value, if we had this element already in the set
			deletedRows.push_back( delRows[index] );
		}

		Eigen::MatrixXcd slacOld = dslac;

		Eigen::MatrixXcd colUpdate( dslac.cols(), 2*updateRank );
		Eigen::MatrixXcd rowUpdate( 2*updateRank, dslac.rows() );
		Eigen::MatrixXcd onesCol = Eigen::MatrixXcd::Zero( updateRank, dslac.rows() );
		Eigen::MatrixXcd onesRow = Eigen::MatrixXcd::Zero( dslac.cols(), updateRank );
		for( size_t index = 0; index < updateRank; index++ ) {

			Eigen::VectorXcd colVec;
			Eigen::RowVectorXcd rowVec;

			if( index < addRows.size() ) {
				// replace added cols/rows by their original value in the full slac operator
				rowVec = fullSlac.row( addRows[index] );
				colVec = fullSlac.col( addCols[index] );
				// do not update, if cols/rows are still deleted, so replace them by their current values.
//				for( auto deletedCol : deletedCols ) {
//					rowVec( deletedCol ) = dslac( addRows[index], deletedCol );
//				}
//				for( auto deletedRow : deletedRows ) {
//					colVec( deletedRow ) = dslac( deletedRow, addCols[index]);
//				}
				onesRow( addRows[index], index ) = 1.;
				onesCol( index, addCols[index] ) = 1.;

//				dslac.col( addCols[index] ) = colVec;
//				dslac.row( addRows[index] ) = rowVec;

			} else {
//				colVec = Eigen::VectorXcd::Zero( dslac.rows() );
//				rowVec = Eigen::RowVectorXcd::Zero( dslac.cols() );

				onesRow( delRows[ index-addRows.size() ], index ) = 1.;
				onesCol( index, delCols[ index-addRows.size() ] ) = 1.;

//				dslac.col( delCols[ index-addRows.size() ] ) = colVec;
//				dslac.row( delRows[ index-addRows.size() ] ) = rowVec;

//				dslac( delRows[index-addRows.size()], delCols[index-addRows.size()] ) = 1.;

				rowVec = -dslac.row( delRows[ index-addRows.size() ] );
				colVec = -dslac.col( delCols[ index-addRows.size() ] );
			}

			colUpdate.col( index+updateRank ) = colVec;
			rowUpdate.row( index ) = rowVec;
		}

//		colUpdate += 0.5*onesRow*subMat;
		std::cout << "first rowUp" << std::endl << rowUpdate << std::endl << std::endl;
		rowUpdate.topRows( updateRank ) = rowUpdate.topRows( updateRank ) + subMat*onesCol;

//		std::cout << "onesRow" << std::endl << onesRow << std::endl << std::endl;
//		std::cout << "onesCol" << std::endl << onesCol << std::endl << std::endl;

		for( size_t i = updateRank; i < 2*updateRank; i++ ) {
			colUpdate.col( i - updateRank ) = onesRow.col( i - updateRank );
			rowUpdate.row( i ) = onesCol.row( i - updateRank );
		}
		std::cout << "subMat" << std::endl << subMat << std::endl << std::endl;
		std::cout << "rowUp" << std::endl << rowUpdate << std::endl << std::endl;
		std::cout << "colUp" << std::endl << colUpdate << std::endl << std::endl;

		std::cout << "outer: " << std::endl << colUpdate * rowUpdate << std::endl << std::endl;
		WoodburyUpdate( colUpdate, rowUpdate );
//		WoodburyUpdate( colUpdate, onesCol );
//		WoodburyUpdate( onesRow, rowUpdate );	// seems like we have to do this update first to avoid det=0, likely because of the subMat update
//		WoodburyUpdate( onesRow*subMat, onesCol );

//		dslac += onesRow*rowUpdate + colUpdate*onesCol;
		dslac += colUpdate * rowUpdate;
//		std::cout << "slacOld:" << std::endl << slacOld << std::endl << std::endl;

//	} else if( delRows.size() != 0 ) {
//		deleteEntries( delRows, delCols );
//	} else if( addRows.size() != 0 ){
//		addEntries( addRows, addCols );
//	}
}

void SlacOperatorMatrix::deleteEntries( std::vector<size_t> rows, std::vector<size_t> cols ) {
	if( rows.size() != cols.size() ) {
		std::cerr << "Trying to delete not the same number of rows/cols. This should be impossible!" << std::endl;
		exit(1);
	}

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

	detVal *= submat.determinant();

	inverse -= rowmat * (submat.inverse() ) * colmat;
	for( size_t index = 0; index < cols.size(); index++ ) {
		inverse( rows[index], cols[index] ) += 1.;

		dslac.col( cols[index] ) = Eigen::VectorXcd::Zero( dslac.rows() );
		dslac.row( rows[index] ) = Eigen::RowVectorXcd::Zero( dslac.cols() );
		dslac( rows[index], cols[index] ) = 1.;

		deletedCols.push_back( cols[index] );		//TODO: check return value, if we had this element already in the set
		deletedRows.push_back( rows[index] );
	}
}

void SlacOperatorMatrix::addEntries( std::vector<size_t> rows, std::vector<size_t> cols ) {
	if( rows.size() != cols.size() ) {
		std::cerr << "Trying to add not the same number of rows/cols. This should be impossible!" << std::endl;
		exit(1);
	}

	if( rows.size() != 0 ) {
		Eigen::MatrixXcd subMat( rows.size(), cols.size() );
		for( size_t colIndex = 0; colIndex < cols.size(); colIndex++ ) {
			for( size_t rowIndex = 0; rowIndex < rows.size(); rowIndex++ ) {
				subMat( rowIndex, colIndex ) = -fullSlac( rows[rowIndex], cols[colIndex] )-dslac( rows[rowIndex], cols[colIndex] );
			}
		}

		for( size_t index = 0; index < cols.size(); index++ ) {
			deletedCols.erase( std::remove( deletedCols.begin(), deletedCols.end(), cols[index] ), deletedCols.end() );
			deletedRows.erase( std::remove( deletedRows.begin(), deletedRows.end(), rows[index] ), deletedRows.end() );
		}

		Eigen::MatrixXcd colUpdate( dslac.cols(), rows.size() );
		Eigen::MatrixXcd rowUpdate( cols.size(), dslac.rows() );
		Eigen::MatrixXcd onesCol = Eigen::MatrixXcd::Zero( cols.size(), dslac.rows() );
		Eigen::MatrixXcd onesRow = Eigen::MatrixXcd::Zero( dslac.cols(), rows.size() );
		for( size_t index = 0; index < cols.size(); index++ ) {
			// replace added cols/rows by their original value in the full slac operator
			Eigen::VectorXcd colVec = fullSlac.col( cols[index] );
			Eigen::RowVectorXcd rowVec = fullSlac.row( rows[index] );

			// do not update, if cols/rows are still deleted, so replace them by their current values.
			for( auto deletedCol : deletedCols ) {
				rowVec( deletedCol ) = dslac( rows[index], deletedCol );
			}
			for( auto deletedRow : deletedRows ) {
				colVec( deletedRow ) = dslac( deletedRow, cols[index]);
			}

			colUpdate.col( index ) = colVec;
			rowUpdate.row( index ) = rowVec;

			onesRow( rows[index], index ) = 1.;
			onesCol( index, cols[index] ) = 1.;

			dslac.col( cols[index] ) = colVec;
			dslac.row( rows[index] ) = rowVec;
		}

		colUpdate += onesRow*subMat;

		WoodburyUpdate( onesRow, rowUpdate );	// seems like we have to do this update first to avoid det=0, likely because of the subMat update
		WoodburyUpdate( colUpdate, onesCol );
	}
}

void SlacOperatorMatrix::WoodburyUpdate( Eigen::MatrixXcd U, Eigen::MatrixXcd V ) {
	Eigen::MatrixXcd Vdinverse = V*inverse;
	Eigen::MatrixXcd subMat = Vdinverse * U + Eigen::MatrixXcd::Identity( V.rows(), U.cols() );
	Eigen::FullPivLU< Eigen::MatrixXcd > subLU( subMat );
	if( subLU.isInvertible() ) {
		detVal *= subLU.determinant();
		inverse -= inverse * U * subLU.inverse() * Vdinverse; // this is much faster than using *= operator with id - U*...
	} else {
		detVal = 0.;
	}
//	std::cout << "Det>" << detVal << std::endl;
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
	deletedCols.clear();
	deletedRows.clear();
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

