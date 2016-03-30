/*
 * FieldBoolean.cpp
 *
 *  Created on: 30.03.2016
 *      Author: dschmidt
 */

#include "FieldBoolean.h"

namespace FermiOwn {

FieldBoolean::FieldBoolean( const size_t latticeVolume,
		const size_t numberOfSpins, const size_t numberOfFlavours,
		std::ranlux48 * rndGen, InitType init ) :
Field<bool>( latticeVolume, numberOfFlavours*(numberOfFlavours+1)/2 * numberOfSpins, rndGen, init ),
numSpins(numberOfSpins),
numFlavours(numberOfFlavours)
{
}

FieldBoolean::~FieldBoolean() {
}

bool FieldBoolean::getValue(size_t x, size_t spin, size_t flavour1, size_t flavour2) const {
	return data( x, colIndex( spin, flavour1, flavour2 ) );
}

void FieldBoolean::setValue(bool val, size_t x, size_t spin, size_t flavour1, size_t flavour2) {
	data( x, colIndex( spin, flavour1, flavour2 ) ) = val;
}

void FieldBoolean::invert(size_t x, size_t spin, size_t flavour1, size_t flavour2) {
	data( x, colIndex( spin, flavour1, flavour2 ) ) = !data( x, colIndex( spin, flavour1, flavour2 ) );
}

bool FieldBoolean::constraintViolated( size_t x ) const {
	for( size_t i = 0; i < numSpins; i++ ) {
		for( size_t m = 0; m < numFlavours; m++ ) {
			int constrSum = 0;
			for( size_t n = 0; n < numFlavours; n++ )
			{
				constrSum += data( x, colIndex( i, m, n ) );
			}
			if( constrSum > 1 ) return true;
		}
	}
	return false;
}

size_t FieldBoolean::sumAll() const {
	return data.count();
}

/*
 * Private Functions
 *****************************************************************************/

size_t FieldBoolean::colIndex(size_t spin, size_t flavour1, size_t flavour2) const {
	// catching errors
	if( spin >= numSpins ) {
		std::cerr << "Spin index of FieldBoolean out of bounds!" << std::endl;
		exit(1);
	} else if ( flavour1 >= numFlavours || flavour2 >= numFlavours) {
		std::cerr << "Flavour index of FieldBoolean out of bounds!" << std::endl;
		exit(1);
	}
	size_t index;
	// calculating flavour-part of the index
	if( flavour1 == flavour2) {
		index = flavour1;
	} else {
		if( flavour1 > flavour2 ) {
			int tmp = flavour2;
			flavour2 = flavour1;
			flavour1 = tmp;
		}
		index = numFlavours+numFlavours*(numFlavours-1)/2 - (numFlavours-flavour1)*(numFlavours-flavour1-1)/2+flavour2-flavour1-1;
	}
	// shifting to correct spin part
	return index + spin*numFlavours*(numFlavours+1)/2;
}

} /* namespace FermiOwn */
