/*
 * ThirringKField.cpp
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#include "ThirringKField.h"

namespace FermiOwn {

ThirringKField::ThirringKField( const size_t latticeVolume, const std::vector<size_t>& internal ) :
		BasicCloneKField( latticeVolume, internal )
{}

ThirringKField::~ThirringKField() {}

void ThirringKField::enforceConstraint( size_t x, const std::vector<size_t>& internal ) {

	if( getValue( x, internal ) ) {
		if( internal[1]==internal[2] ) {
			for( size_t flavour = 0; flavour < internalRanges[1]; flavour++ ) {
				if( flavour != internal[1] ) {
					setValue( false, x, internal );
					setValue( false, x, internal );
				}
			}
		} else {
			setValue( false, x, internal );
			setValue( false, x, internal );
		}
	}
}
bool ThirringKField::constraintViolated( size_t x ) const {

	//checking kxab = kxba
	for( size_t flavour1 = 0; flavour1 < internalRanges[1]; flavour1++ ) {
		for( size_t flavour2 = 0; flavour2 < flavour1; flavour2++ ) {
			int abSum = 0;
			int baSum = 0;
			for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
				abSum += getValue( x, {spin, flavour1, flavour2} );
				baSum += getValue( x, {spin, flavour2, flavour1} );
			}
			if( abSum != baSum ) {
				//				std::cerr << "Constraint kxab=kxba violated!" << std::endl;
				return true;
			}
		}
	}

	// checking sum over first flavour for all other flavours and spins
	for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
		for( size_t flavour2 = 0; flavour2 < internalRanges[2]; flavour2++ ) {
			int constrSum = 0;
			for( size_t flavour1 = 0; flavour1 < internalRanges[1]; flavour1++ ) {
				constrSum += getValue( x, {spin, flavour1, flavour2} );
			}
			if( constrSum > 1 ) return true;
		}
	}

	// checking sum over second flavour for all other flavours and spins
	for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
		for( size_t flavour1 = 0; flavour1 < internalRanges[1]; flavour1++ ) {
			int constrSum = 0;
			for( size_t flavour2 = 0; flavour2 < internalRanges[2]; flavour2++ )
			{
				constrSum += getValue( x, {spin, flavour1, flavour2} );
			}
			if( constrSum > 1 ) return true;
		}
	}

	int count = 0;
	for( size_t flavour1 = 0; flavour1 < internalRanges[1]; flavour1++ ) {
		int constrSum = 0;
		for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
			constrSum += getValue( x, {spin, flavour1, flavour1} );
		}
		if( 1 == constrSum ) {
			count++;
		}
	}
//	std::cout << "count = " << count << std::endl;
	if( count % 2 != 0 ) {
//		std::cout << "returning true" << std::endl;
		return true;
	}

	return false;
}

size_t ThirringKField::sumAll() const {
	return data.count();
}

Eigen::ArrayXi ThirringKField::countSummedSpin(size_t x) const {
	if( x >= V ) {
		std::cerr << "Error: Current point x=" << x << " is outside the lattice volume of " << V << std::endl;
		exit(1);
	}
	RowVectorXb row = data.row(x);
	Eigen::ArrayXi count = Eigen::ArrayXi::Zero( internalRanges[0] + 1 );
	for( size_t flavour = 0; flavour < internalRanges[1]; flavour++ ) {
		size_t spinSum = 0;
		for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
			spinSum += row( colIndex( {spin, flavour, flavour} ) );
		}
		if( spinSum > internalRanges[0] ) {
			std::cerr << spinSum << " is too much spins counted in ThirringKField. Should be less than " << internalRanges[0] << std::endl;
			exit(1);
		}
		count[spinSum]++;
	}
	return count;
}

size_t ThirringKField::countOffdiagonal2() const {
	size_t count = 0;
	for( size_t x = 0; x < V; x++ ) {
		for( size_t flavour2 = 0; flavour2 < internalRanges[2]; flavour2++ ) {
			for( size_t flavour1 = 0; flavour1 < flavour2; flavour1++ ) {
				size_t spinSum = 0;
				for( size_t spin = 0; spin < internalRanges[0]; spin++ ) {
					spinSum += getValue( x, {spin, flavour1, flavour2} );
				}
				if( spinSum == 2 ) count++;
			}
		}
	}
	return count;
}

ThirringKField ThirringKField::different( ThirringKField other ) const {
	ThirringKField difference = other;
	difference.setZero();
	difference.data = data.cwiseNotEqual(other.data);
	return difference;
}

size_t ThirringKField::colIndex( const std::vector<size_t>& internal ) const {
	// catching errors
	if( internal[0] >= internalRanges[0] ) {
		std::cerr << "Spin index of FieldBoolean out of bounds!" << std::endl;
		exit(1);
	} else if ( internal[1] >= internalRanges[1] || internal[2] >= internalRanges[2] ) {
		std::cerr << "Flavour index of FieldBoolean out of bounds!" << std::endl;
		exit(1);
	}
	// calculating flavour-part of the index
	size_t index = internal[2] + internalRanges[1] * internal[1] + internalRanges[1]*internalRanges[2]*internal[0];
	if( index >= DoFperX ) {
		std::cerr << "Col index too large!" << std::endl;
		exit(1);
	}
	return index;
}

} /* namespace FermiOwn */
