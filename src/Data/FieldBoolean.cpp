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
								Field<bool>( latticeVolume, numberOfFlavours*numberOfFlavours* numberOfSpins, rndGen, init ),
								numSpins(numberOfSpins),
								numFlavours(numberOfFlavours),
								numColsPerSpin(numFlavours*numFlavours)
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

	//checking kxab = kxba

	for( size_t flavour1 = 0; flavour1 < numFlavours; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < flavour1; flavour2++ ) {
			int abSum = 0;
			int baSum = 0;
			for( size_t spin = 0; spin < numSpins; spin++ ) {
				abSum += data( x, colIndex( spin, flavour1, flavour2 ) );
				baSum += data( x, colIndex( spin, flavour2, flavour1 ) );
			}
			if( abSum != baSum ) {
				std::cerr << "Constraint kxab=kxba violated!" << std::endl;
				return true;
			}
		}
	}

	// checking sum over second flavour for all other flavours and spins
	for( size_t spin = 0; spin < numSpins; spin++ ) {
		for( size_t flavour1 = 0; flavour1 < numFlavours; flavour1++ ) {
			int constrSum = 0;
			for( size_t flavour2 = 0; flavour2 < numFlavours; flavour2++ )
			{
				constrSum += data( x, colIndex( spin, flavour1, flavour2 ) );
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
			spinSum += row( colIndex( spin, flavour, flavour) );
		}
		if( spinSum > numSpins ) {
			std::cerr << spinSum << " is too much spins counted in FieldBoolean. Should be less than " << numSpins << std::endl;
			exit(1);
		}
		count[spinSum]++;
	}
	return count;
}

size_t FieldBoolean::countOffdiagonal2() const {
	size_t count = 0;
	for( size_t x = 0; x < latVol; x++ ) {
		for( size_t flavour2 = 0; flavour2 < numFlavours; flavour2++ ) {
			for( size_t flavour1 = 0; flavour1 < flavour2; flavour1++ ) {
				size_t spinSum = 0;
				for( size_t spin = 0; spin < numSpins; spin++ ) {
					spinSum += data( x, colIndex( spin, flavour1, flavour2 ) );
				}
				if( spinSum == 2 ) count++;
			}
		}
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
	index = flavour2 + numFlavours * flavour1;
	// shifting to correct spin part
	index += spin*numColsPerSpin;
	if( data.cols() <= index ) {
		std::cerr << "Col index too large!" << std::endl;
		exit(1);
	}
	return index;
}

} /* namespace FermiOwn */
