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
numFlavours(numberOfFlavours),
numColsPerSpin(numFlavours*(numFlavours+1)/2)
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

Eigen::ArrayXi FieldBoolean::countSummedSpin(size_t x) const {
	auto row = data.row(x);
	Eigen::ArrayXi count = Eigen::ArrayXi::Zero( numSpins+1 );
	for( size_t flavour = 0; flavour < numFlavours; flavour++ ) {
		size_t spinSum = 0;
		for( size_t spin = 0; spin < numSpins; spin++ ) {
			spinSum += row( spin*numColsPerSpin + flavour );
		}
		if( spinSum > numSpins ) {
			std::cerr << spinSum << " is too much spins counted in FieldBoolean. Should be less than " << numSpins << std::endl;
			exit(1);
		}
		count[spinSum]++;
	}
	return count;
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
		index = numColsPerSpin - (numFlavours-flavour1)*(numFlavours-flavour1-1)/2+flavour2-flavour1-1;		//TODO: index can be simplified!
	}
	// shifting to correct spin part
	return index + spin*numColsPerSpin;
}

} /* namespace FermiOwn */
