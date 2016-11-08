/*
 * MatrixChanges.cpp
 *
 *  Created on: 08.11.2016
 *      Author: dschmidt
 */

#include "MatrixChanges.h"

namespace FermiOwn {

MatrixChanges::MatrixChanges( const DSlashUpdater& dslashUp, const ThirringKField& booleanField ) :
		dslash( dslashUp ),
		kfield( booleanField ),
		internalRanges( kfield.getInternalRanges() ),
		V( kfield.getVolume() )
	{}

MatrixChanges::~MatrixChanges() {}

const AddDelRowCol&  MatrixChanges::calculateDifference( const ThirringKField& other ) {
	// construct a vector filled with range 0 ... V-1 and construct a set from it
	idxVec range( V );
	std::iota( range.begin(), range.end(), 0 );
	std::set<size_t> changedAt( range.begin(), range.end() );
	return calculateDifference( changedAt, other );
}

const AddDelRowCol&  MatrixChanges::calculateDifference(const std::set<size_t>& changedAt, const ThirringKField& other ) {
	indices.delRows.clear();
	indices.delCols.clear();
	indices.addRows.clear();
	indices.addCols.clear();

	for( size_t flavour1 = 0; flavour1 < internalRanges[1]; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < internalRanges[2]; flavour2++ ) {
			for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
				for( size_t x : changedAt ) {
					size_t currVal = kfield.getValue( x, {spin, flavour1, flavour2} );
					if( currVal != other.getValue( x, {spin, flavour1, flavour2} ) ) {
						if( currVal ) {
							indices.delRows.push_back( dslash.matIndex(x, spin, flavour1) );
							indices.delCols.push_back( dslash.matIndex(x, spin, flavour2) );
						} else {
							indices.addRows.push_back( dslash.matIndex(x, spin, flavour1) );
							indices.addCols.push_back( dslash.matIndex(x, spin, flavour2) );
						}
					}
				}
			}
		}
	}
	return indices;
}

void MatrixChanges::Print() const {
	std::cout << "Current indices to change: " << std::endl << "delCols:";
	for( size_t col : indices.delCols ) std::cout << col << " ";
	std::cout << std::endl;
	std::cout << "delRows:";
	for( size_t row : indices.delRows ) std::cout << row << " ";
	std::cout << std::endl;
	std::cout << "addCols:";
	for( size_t col : indices.addCols ) std::cout << col << " ";
	std::cout << std::endl;
	std::cout << "addRows:";
	for( size_t row : indices.addRows ) std::cout << row << " ";
	std::cout << std::endl;
}

} /* namespace FermiOwn */
