/*
 * MatrixChangesBase.cpp
 *
 *  Created on: 09.11.2016
 *      Author: dschmidt
 */

#include "MatrixChanges.h"

namespace FermiOwn {

	template<> const AddDelRowCol&  MatrixChanges<ThirringKField>::calculateDifference(const std::set<size_t>& changedAt, const ThirringKField& other ) {
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

} /* end namespace FermiOwn */
