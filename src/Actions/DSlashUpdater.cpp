/*
 * DSlashUpdater.cpp
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#include "DSlashUpdater.h"

namespace FermiOwn {

DSlashUpdater::DSlashUpdater() {
	// TODO Auto-generated constructor stub

}

DSlashUpdater::DSlashUpdater( size_t Nt, size_t Ns, size_t dim, size_t numFlavours ) :
					N(Nt*Ns*Ns),
					dimSpinor(2),
					Nf(numFlavours),
					matSize( N*dimSpinor*Nf ),
					oldMat( matSize ),
					currMat( matSize ),
					currentRows( matSize ),
					currentCols( matSize )
{
	using namespace Eigen;
	if( dim != 3 ) {
		std::cerr << "DSlashUpdater implemented only in 3D." << std::endl;
		exit(1);
	}

	CliffordAlgebra cliff;
	std::vector< Matrix2cd > gamma = cliff.getGammas();
	MatrixXcd dx1D= makeSlac1D(Nt);
	MatrixXcd dyz = makeSlac1D(Ns);
	MatCoeffList coeffs;
	for( size_t flavour = 0; flavour < Nf; flavour++ ) {
		for( size_t spin1 = 0; spin1 < dimSpinor; spin1++ ) {
			for( size_t spin2 = 0; spin2 < dimSpinor; spin2++ ) {
				size_t internRowIndex = N*( dimSpinor*flavour + spin1 );
				size_t internColIndex = N*( dimSpinor*flavour + spin2 );

				// Matrix elements in x-direction
				Complex matElem = gamma[0]( spin1, spin2 );
				size_t rowIndex;
				size_t colIndex;
				if( std::fabs( matElem ) > ZERO_TOL ) {
					for( size_t r = 0; r < Ns*Ns; r++ ) {
						for( size_t v = 0; v < Nt; v++ ) {
							for( size_t u = v+1; u < Nt; u++ ) {
								rowIndex = internRowIndex + Nt*r + v;
								colIndex = internColIndex + Nt*r + u;
								// ensure, we write on the upper triangular part only
								if( rowIndex > colIndex ) {
									coeffs.push_back( Triplet<Complex>( colIndex, rowIndex, std::conj( Complex( I*matElem*dx1D(v,u) ) ) ) );
								} else {
									coeffs.push_back( Triplet<Complex>( rowIndex, colIndex, I*matElem*dx1D(v,u) ) );
								}
							}
						}
					}
				}

				matElem = gamma[1]( spin1, spin2 );
				if( std::fabs( matElem ) > ZERO_TOL ) {
					for( size_t q = 0; q < Ns; q++ ) {
						for( size_t v = 0; v < Nt; v++ ) {
							for( size_t r = 0; r < Ns; r++ ) {
								for( size_t s = r+1; s < Ns; s++ ) {
									rowIndex = internRowIndex + Ns*Nt*q + Nt*r + v;
									colIndex = internColIndex + Ns*Nt*q + Nt*s + v;
									// ensure, we write on the upper triangular part only
									if( rowIndex > colIndex ) {
										coeffs.push_back( Triplet<Complex>( colIndex, rowIndex, std::conj( Complex( I*matElem*dyz(r,s) ) ) ) );
									} else {
										coeffs.push_back( Triplet<Complex>( rowIndex, colIndex, I*matElem*dyz(r,s) ) );
									}
								}
							}
						}
					}
				}

				matElem = gamma[2]( spin1, spin2 );
				if( std::fabs( matElem ) > ZERO_TOL ) {
					for( size_t v = 0; v < Ns*Nt; v++ ) {
						for( size_t r = 0; r < Ns; r++ ) {
							for( size_t s = r+1; s < Ns; s++ ) {
								rowIndex = internRowIndex + Ns*Nt*r + v;
								colIndex = internColIndex + Ns*Nt*s + v;
								// ensure, we write on the upper triangular part only
								if( rowIndex > colIndex ) {
									coeffs.push_back( Triplet<Complex>( colIndex, rowIndex, std::conj( Complex( I*matElem*dyz(r,s) ) ) ) );
								} else {
									coeffs.push_back( Triplet<Complex>( rowIndex, colIndex, I*matElem*dyz(r,s) ) );
								}
							}
						}
					}
				}

			}
		}
	}

	currMat.setFromCoeffList( coeffs, true );
	fullOperator = currMat.getMatrix();

	size_t i = 0;
	std::iota( currentRows.begin(), currentRows.end(), i );
	std::iota( currentCols.begin(), currentCols.end(), i );

	for( auto num : currentRows ) std::cout << num << " ";
	for( auto num : currentCols ) std::cout << num << " ";

	updated = true;
}


DSlashUpdater::~DSlashUpdater() {
	// TODO Auto-generated destructor stub
}

void DSlashUpdater::calculateUpdateMatrices( const FieldBoolean& kxiab, const FieldBoolean& change ) {

	// check what changed, and if it requires adding or deleting row/col pairs
	std::vector<size_t> addRows;
	std::vector<size_t> delRows;
	std::vector<size_t> addCols;
	std::vector<size_t> delCols;
	std::vector<size_t> targetRows = currentRows;
	std::vector<size_t> targetCols = currentCols;

	for( size_t flavour1 = 0; flavour1 < Nf; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < Nf; flavour2++ ) {
			for( size_t spin = 0; spin < dimSpinor; spin++ ) {
				for( size_t x = 0; x < N; x++ ) {
					if( change.getValue( x, spin, flavour1, flavour2 ) ) {
						if( kxiab.getValue( x, spin, flavour1, flavour2 ) ) {
							delRows.push_back( matIndex(x, spin, flavour2) );
							delCols.push_back( matIndex(x, spin, flavour1) );
							targetRows.erase( std::remove( targetRows.begin(), targetRows.end(), delRows.back() ) );
							targetCols.erase( std::remove( targetCols.begin(), targetCols.end(), delCols.back() ) );
						} else {
							addRows.push_back( matIndex(x, spin, flavour2) );
							addCols.push_back( matIndex(x, spin, flavour1) );
							targetRows.push_back( addRows.back() );
							targetCols.push_back( addCols.back() );
						}
					}
				}
			}
		}
	}

	//TODO: implement something to find changes that add and delete the same row/col
	// std::set_intersection seems to need sorted vectors, but we have to keep the order...

	size_t delRank = delRows.size();
	//TODO: use a more efficient update, if addRank == 0, to be implemented in WoodburyMatrix
	size_t addRank = addRows.size();
	size_t updateRank = delRank + addRank;

	MatCoeffList rowCoeffs;
	MatCoeffList colCoeffs;

	SparseMat rowUpdate( 2*updateRank, matSize );
	SparseMat colUpdate( matSize, 2*updateRank );

	const SparseMat & curr = currMat.getMatrix();

	for( size_t i = 0; i < 2*updateRank; i++ ) {
		if( i < updateRank ) {
			// we construct the matrix elements from current and full operator, iterating over cols to be replaced
			if( i < delRank ) {
				for( size_t j = 0; j < currentCols.size(); j++ ) {
					// we perform deletions: elements are the negative of the current values.
					Complex insertElem = curr.coeff( delCols[i], currentCols[j] );

					// coeff should return an exact complex 0 if the element does not exist, so the following checks should be working.
					if( insertElem != Complex( 0., 0. ) ) rowCoeffs.push_back( CoeffTriplet( i, currentCols[j], -insertElem ) );
				}

				// we fill the first part of col update with 1s as counterpart for row matrix
				colCoeffs.push_back( CoeffTriplet( delRows[i], i, Complex( 1., 0. ) ) );
			} else {
				for( size_t j = 0; j < targetCols.size(); j++ ) {
					// we restore rows, so we insert the full value and remove the current
					Complex insertElem = fullOperator.coeff( addRows[i-delRank], targetCols[j] ) - curr.coeff( addRows[i-delRank], targetCols[j] );
					if( insertElem != Complex( 0., 0. ) ) rowCoeffs.push_back( CoeffTriplet( i, targetCols[j], insertElem ) );
				}

				colCoeffs.push_back( CoeffTriplet( addRows[i-delRank], i, Complex( 1., 0. ) ) );
			}
		} else {
			// we fill the second part of the row matrix with 1s as counterpart for the row updates.
			//TODO: can be simplified by pushing add/del in same list and only keeping delRank to distinguish add/del, possibly saves next if
			if( i < updateRank + delRank ) {
				rowCoeffs.push_back( CoeffTriplet( i, delCols[i-updateRank], Complex(1.,0.) ) );

				// the second part of the col matrix gets entries from current and full operator
				for( size_t j = 0; j < currentRows.size(); j++ ) {
					// check, if we already put the update into the row matrix
					if( std::find( delRows.begin(), delRows.end(), currentRows[j] ) != delRows.end() ) {
						// crossing entries need to get a 1
						if( delRows[i-updateRank] == currentRows[j] ) colCoeffs.push_back( CoeffTriplet( currentRows[j], i, Complex( 1., 0.) ) );
					} else {
						Complex insertElem = curr.coeff( currentRows[j], delCols[i-updateRank] );
						if( insertElem != Complex( 0., 0. ) ) colCoeffs.push_back( CoeffTriplet( currentRows[j], i, -insertElem ) );
					}
				}
			} else {
				rowCoeffs.push_back( CoeffTriplet( i, addCols[i-updateRank-delRank], Complex(1.,0.) ) );

				for( size_t j = 0; j < targetRows.size(); j++ ) {
					// only insert entries, if we didn't do it in the row update
					if( std::find( addRows.begin(), addRows.end(), targetRows[j]) == addRows.end() ) {
						Complex insertElem = fullOperator.coeff( targetRows[j], addCols[i-updateRank-delRank] ) - curr.coeff( targetRows[j], addCols[i-updateRank-delRank] );
						if( insertElem != Complex( 0., 0. ) ) colCoeffs.push_back( CoeffTriplet( targetRows[j], i ) );
					}
				}
			}
		}
	}

	rowUpdate.setFromTriplets( rowCoeffs.begin(), rowCoeffs.end() );
	colUpdate.setFromTriplets( colCoeffs.begin(), colCoeffs.end() );

	std::cout << rowUpdate << std::endl << std::endl << colUpdate << std::endl;
}

Eigen::MatrixXcd DSlashUpdater::makeSlac1D( size_t size ) const {
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
