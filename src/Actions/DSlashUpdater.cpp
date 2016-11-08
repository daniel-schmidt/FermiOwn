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
	currentCols( matSize ),
	changed( false ),
	oldMatNeedsUpdate( true )
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
	oldMat = currMat;
	oldMatNeedsUpdate = false;
	size_t i = 0;
	std::iota( currentRows.begin(), currentRows.end(), i );
	std::iota( currentCols.begin(), currentCols.end(), i );
}


DSlashUpdater::~DSlashUpdater() {
	// TODO Auto-generated destructor stub
}

//void DSlashUpdater::calculateUpdateMatrices( const ThirringKField& kxiab, const ThirringKField& change ) {
void DSlashUpdater::calculateUpdateMatrices( const AddDelRowCol& changes ) {
	if( changed == true ) std::cerr << "Warning: DSlashUpdater overwrites previously not executed changes. Better use keep/reset before!" << std::endl;
	// check what changed, and if it requires adding or deleting row/col pairs

//	std::vector<size_t> addRows;
//	std::vector<size_t> delRows;
//	std::vector<size_t> addCols;
//	std::vector<size_t> delCols;
	targetRows = currentRows;
	targetCols = currentCols;


//	for( size_t flavour1 = 0; flavour1 < Nf; flavour1++ ) {
//		for( size_t flavour2 = 0; flavour2 < Nf; flavour2++ ) {
//			for( size_t spin = 0; spin < dimSpinor; spin++ ) {
//				for( size_t x = 0; x < N; x++ ) {
//					if( change.getValue( x, {spin, flavour1, flavour2} ) ) {
//						if( kxiab.getValue( x, {spin, flavour1, flavour2} ) ) {
//							delRows.push_back( matIndex(x, spin, flavour1) );
//							delCols.push_back( matIndex(x, spin, flavour2) );
//						} else {
//							addRows.push_back( matIndex(x, spin, flavour1) );
//							addCols.push_back( matIndex(x, spin, flavour2) );
//						}
//					}
//				}
//			}
//		}
//	}


	//TODO: implement something to find changes that add and delete the same row/col
	// std::set_intersection seems to need sorted vectors, but we have to keep the order...

	std::vector<size_t> addRowsSorted = changes.addRows;
	std::vector<size_t> delRowsSorted = changes.delRows;
	std::vector<size_t> addColsSorted = changes.addCols;
	std::vector<size_t> delColsSorted = changes.delCols;

	std::sort( addRowsSorted.begin(), addRowsSorted.end() );
	std::sort( delRowsSorted.begin(), delRowsSorted.end() );
	std::sort( addColsSorted.begin(), addColsSorted.end() );
	std::sort( delColsSorted.begin(), delColsSorted.end() );

	std::vector<size_t> intersection;
	std::vector<size_t> trueAddRows = changes.addRows;
	std::vector<size_t> trueDelRows = changes.delRows;
	std::vector<size_t> trueAddCols = changes.addCols;
	std::vector<size_t> trueDelCols = changes.delCols;

	std::set_intersection( addRowsSorted.begin(), addRowsSorted.end(), delRowsSorted.begin(), delRowsSorted.end(), std::back_inserter( intersection ) );

	typedef std::pair< size_t, size_t > idxPair;
	std::set< idxPair > addOnes;
	std::set< idxPair > rmOnes;

	if( intersection.size() != 0 ) {
//		std::cout << "Values in rows for add and delete: ";
		for( auto val : intersection ) {
//			std::cout << val << " ";
			trueDelRows.erase( std::remove( trueDelRows.begin(), trueDelRows.end(), val ) );
			trueAddRows.erase( std::remove( trueAddRows.begin(), trueAddRows.end(), val ) );
			size_t index = std::find( changes.delRows.begin(), changes.delRows.end(), val ) - changes.delRows.begin();
			addOnes.insert( idxPair( changes.delRows[index], changes.delCols[index] ) );

			index = std::find( changes.addRows.begin(), changes.addRows.end(), val ) - changes.addRows.begin();
			rmOnes.insert( idxPair( changes.addRows[index], changes.addCols[index] ) );
		}
//		std::cout << std::endl;
	}

	intersection.clear();
	std::set_intersection( addColsSorted.begin(), addColsSorted.end(), delColsSorted.begin(), delColsSorted.end(), std::back_inserter( intersection ) );

	if( intersection.size() != 0 ) {
//		std::cout << "Values in cols for add and delete: ";
		for( auto val : intersection ) {
//			std::cout << val << " ";
			trueDelCols.erase( std::remove( trueDelCols.begin(), trueDelCols.end(), val ) );
			trueAddCols.erase( std::remove( trueAddCols.begin(), trueAddCols.end(), val ) );

			size_t index = std::find( changes.delCols.begin(), changes.delCols.end(), val ) - changes.delCols.begin();
			addOnes.insert( idxPair( changes.delRows[index], changes.delCols[index] ) );

			index = std::find( changes.addCols.begin(), changes.addCols.end(), val ) - changes.addCols.begin();
			rmOnes.insert( idxPair( changes.addRows[index], changes.addCols[index] ) );

		}
//		std::cout << std::endl;
	}

	for( auto col : changes.addCols ) {
		targetCols.push_back( col );
	}
	for( auto row : changes.addRows ) {
		targetRows.push_back( row );
	}
	for( auto col : changes.delCols ) {
		targetCols.erase( std::remove( targetCols.begin(), targetCols.end(), col ) );
	}
	for( auto row : changes.delRows ) {
		targetRows.erase( std::remove( targetRows.begin(), targetRows.end(), row) );
	}

//	std::cout << "addRows: ";
//	for( auto ar : addRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "addCols: ";
//	for( auto ar : addCols ) std::cout << ar << " ";
//	std::cout << std::endl;
//
//	std::cout << "delRows: ";
//	for( auto ar : delRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "delCols: ";
//	for( auto ar : delCols ) std::cout << ar << " ";
//	std::cout << std::endl;
//
//	std::cout << "trueAddRows: ";
//	for( auto ar : trueAddRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "trueAddCols: ";
//	for( auto ar : trueAddCols ) std::cout << ar << " ";
//	std::cout << std::endl;
//
//	std::cout << "trueDelRows: ";
//	for( auto ar : trueDelRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "trueDelCols: ";
//	for( auto ar : trueDelCols ) std::cout << ar << " ";
//	std::cout << std::endl;
//
//	std::cout << "addOnes: ";
//	for( auto ar : addOnes ) std::cout << "(" << ar.first << "," << ar.second << ") ";
//	std::cout << std::endl;
//	std::cout << "rmOnes: ";
//	for( auto ar : rmOnes ) std::cout << "(" << ar.first << "," << ar.second << ") ";
//	std::cout << std::endl;
//
//	std::cout << "currentRows: ";
//	for( auto ar : currentRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "currentCols: ";
//	for( auto ar : currentCols ) std::cout << ar << " ";
//	std::cout << std::endl;
//
//	std::cout << "targetRows: ";
//	for( auto ar : targetRows ) std::cout << ar << " ";
//	std::cout << std::endl;
//	std::cout << "targetCols: ";
//	for( auto ar : targetCols ) std::cout << ar << " ";
//	std::cout << std::endl;

	// Starting actual update

	//	size_t delRank = delRows.size();
	//TODO: use a more efficient update, if addRank == 0, to be implemented in WoodburyMatrix
	size_t addRank = changes.addRows.size() + changes.addCols.size();
	//	size_t updateRank = delRank + addRank;

	size_t updateRank = trueDelRows.size() + trueAddRows.size() + trueDelCols.size() + trueAddCols.size() + rmOnes.size() + addOnes.size();

	if( updateRank > 0 ) {
		MatCoeffList rowCoeffs;
		MatCoeffList colCoeffs;

		rowUpdate.resize( updateRank, matSize );
		colUpdate.resize( matSize, updateRank );

		const SparseMat & curr = currMat.getMatrix();

		size_t t = 0;

		for( size_t i = 0; i < trueDelRows.size(); i++ ) {
			for( size_t j = 0; j < currentCols.size(); j++ ) {
				// we perform deletions: elements are the negative of the current values.
				Complex insertElem = curr.coeff( trueDelRows[i], currentCols[j] );

				// coeff should return an exact complex 0 if the element does not exist, so the following checks should be working.
				if( insertElem != Complex( 0., 0. ) ) rowCoeffs.push_back( CoeffTriplet( t, currentCols[j], -insertElem ) );
			}
			// we fill the first part of col update with 1s as counterpart for row matrix
			colCoeffs.push_back( CoeffTriplet( trueDelRows[i], t, Complex( 1., 0. ) ) );
			t++;
		}

		for( size_t i = 0; i < trueAddRows.size(); i++ ) {
			for( size_t j = 0; j < targetCols.size(); j++ ) {
				// we restore rows, so we insert the full value and remove the current
				Complex insertElem = fullOperator.coeff( trueAddRows[i], targetCols[j] ) - curr.coeff( trueAddRows[i], targetCols[j] );
				if( insertElem != Complex( 0., 0. ) ) {
					rowCoeffs.push_back( CoeffTriplet( t, targetCols[j], insertElem ) );
				}
			}

			colCoeffs.push_back( CoeffTriplet( trueAddRows[i], t, Complex( 1., 0. ) ) );
			t++;
		}

		for( size_t i = 0; i < trueDelCols.size(); i++ ) {
			rowCoeffs.push_back( CoeffTriplet( t, trueDelCols[i], Complex(1.,0.) ) );

			// the second part of the col matrix gets entries from current and full operator
			for( size_t j = 0; j < currentRows.size(); j++ ) {
				// check, if we already put the update into the row matrix
				if( std::find( trueDelRows.begin(), trueDelRows.end(), currentRows[j] ) != trueDelRows.end() ) {
					// find crossing entries from original list, needs to get a 1
					size_t pairIdx = std::find( changes.delCols.begin(), changes.delCols.end(), trueDelCols[i] ) - changes.delCols.begin();
					if( changes.delRows[pairIdx] == currentRows[j] )

						colCoeffs.push_back( CoeffTriplet( currentRows[j], t, Complex( 1., 0.) ) );
				} else {
					Complex insertElem = curr.coeff( currentRows[j], trueDelCols[i] );
					if( insertElem != Complex( 0., 0. ) ) colCoeffs.push_back( CoeffTriplet( currentRows[j], t, -insertElem ) );
				}
			}
			t++;
		}

		for( size_t i = 0; i < trueAddCols.size(); i++ ) {
			rowCoeffs.push_back( CoeffTriplet( t, trueAddCols[i], Complex(1.,0.) ) );

			for( size_t j = 0; j < targetRows.size(); j++ ) {
				// only insert entries, if we didn't do it in the row update
				if( std::find( trueAddRows.begin(), trueAddRows.end(), targetRows[j]) == trueAddRows.end() ) {
					Complex insertElem = fullOperator.coeff( targetRows[j], trueAddCols[i] ) - curr.coeff( targetRows[j], trueAddCols[i] );
					if( insertElem != Complex( 0., 0. ) ) colCoeffs.push_back( CoeffTriplet( targetRows[j], t, insertElem ) );
				}
			}
			t++;
		}

		for( idxPair rmIdx : rmOnes ) {
			rowCoeffs.push_back( CoeffTriplet( t, rmIdx.second, Complex( -1., 0. ) ) );
			colCoeffs.push_back( CoeffTriplet( rmIdx.first, t, Complex( 1., 0. ) ) );
			t++;
		}

		for( idxPair addIdx : addOnes ) {
			rowCoeffs.push_back( CoeffTriplet( t, addIdx.second, Complex( 1., 0. ) ) );
			colCoeffs.push_back( CoeffTriplet( addIdx.first, t, Complex( 1., 0. ) ) );
			t++;
		}

		rowUpdate.setFromTriplets( rowCoeffs.begin(), rowCoeffs.end() );
		colUpdate.setFromTriplets( colCoeffs.begin(), colCoeffs.end() );

//			std::cout << "rowUpdate" << std::endl << rowUpdate << std::endl << std::endl << "colUpdate:"<< std::endl << colUpdate << std::endl;

		if( oldMatNeedsUpdate ) {
			oldMat = currMat;
			oldMatNeedsUpdate = false;
		}

		if( addRank > 0 ) {
			currMat.setUpdateMatrices( &colUpdate, &rowUpdate );
		} else {
			currMat.setUpdateMatricesDeleteOnly( &colUpdate, &rowUpdate, changes.delCols, changes.delRows );
		}

		changed = true;
	} // updateRank > 0
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
