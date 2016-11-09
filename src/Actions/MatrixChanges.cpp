/*
 * MatrixChangesBase.cpp
 *
 *  Created on: 09.11.2016
 *      Author: dschmidt
 */

#include "MatrixChanges.h"

namespace FermiOwn {

template<> const AddDelRowCol&  MatrixChanges<ThirringKField>::calculateDifference(const std::set<size_t>& changedAt, const ThirringKField& other ) {
	clearIndices();
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

template<> const AddDelRowCol& MatrixChanges<GrossNeveuKField>::calculateDifference( const std::set<size_t>& changedAt, const GrossNeveuKField& other ) {
	clearIndices();
	for( size_t flavour = 0; flavour < internalRanges[1]; flavour++ ) {
		for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
			for( size_t x: changedAt ) {
				idxVec internal = {spin, flavour};
				size_t currVal = kfield.getValue( x, internal );
				if( currVal != other.getValue( x, internal ) ) {
					size_t idx = dslash.matIndex( x, spin, flavour );
//					std::cout << "x=" << x << " spin=" << spin << " flavour=" << flavour << " idx=" << idx << std::endl;
					if( currVal ) {
						indices.delRows.push_back( idx );
						indices.delCols.push_back( idx );
					} else {
						indices.addRows.push_back( idx );
						indices.addCols.push_back( idx );
					}
				}
			}
		}
	}
	return indices;
}

} /* end namespace FermiOwn */
