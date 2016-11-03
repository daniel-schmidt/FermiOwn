/*
 * GrossNeveuKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_GROSSNEVEUKFIELD_H_
#define SRC_DATA_GROSSNEVEUKFIELD_H_

#include <iostream>

namespace FermiOwn {

class GrossNeveuKField {
public:
	GrossNeveuKField( const size_t latticeVolume, const size_t numberOfSpin, const size_t numberOfFlavours );
	virtual ~GrossNeveuKField();

	inline void setValue( bool val, size_t x, size_t spin, size_t flavour );
	inline bool getValue( size_t x, size_t spin, size_t flavour ) const;
	inline void setRow( const RowVectorXb & newRow, size_t x );
	inline void invert( size_t x, size_t spin, size_t flavour );

	inline void Print() const;

private:
	inline size_t colIndex( size_t spin, size_t flavour) const;

	size_t V;
	size_t numSpin;
	size_t Nf;
	size_t DoFperX;

	MatrixXb data;
};

/*
 * Implementation of functions
 **********************************************************/

GrossNeveuKField::GrossNeveuKField( const size_t latticeVolume, const size_t numberOfSpin, const size_t numberOfFlavours ) :
	V( latticeVolume ),
	numSpin( numberOfSpin ),
	Nf( numberOfFlavours ),
	DoFperX( numberOfSpin*numberOfFlavours ){
	data = MatrixXb::Zero( V, DoFperX );
}

GrossNeveuKField::~GrossNeveuKField() {
}

void GrossNeveuKField::setValue( bool val, size_t x, size_t spin, size_t flavour ) {
	//TODO: define preprocessor variable to disable out-of-bounds checks
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::setValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( spin, flavour ) ) = val;
}

inline bool GrossNeveuKField::getValue( size_t x , size_t spin, size_t flavour ) const {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::getValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	return data( x, colIndex( spin, flavour ) );
}

inline void GrossNeveuKField::setRow( const RowVectorXb & newRow, size_t x ) {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::setRow: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	} else if( newRow.cols() > int( DoFperX ) ) {
		std::cerr << "Error in GrossNeveuKField::setRow: Row length " << newRow.cols() << " is longer than the number of degrees of freedom per site: " << DoFperX << std::endl;
		exit(1);
	}
	data.row(x) = newRow;
}

inline void GrossNeveuKField::invert( size_t x, size_t spin, size_t flavour ) {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::invert: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( spin, flavour ) ) = !data( x, colIndex( spin, flavour ) );
}

inline void GrossNeveuKField::Print() const {
	std::cout << data << std::endl;
}

/*
 * Private functions
 ********************************************************/

inline size_t GrossNeveuKField::colIndex( size_t spin, size_t flavour ) const {
	size_t idx = numSpin * flavour + spin;
	if( idx >= DoFperX ) {
		std::cerr << "Index " << idx << " of bound (" << DoFperX << ") in GrossNeveuKField for spin=" << spin << " and flavour=" << flavour << std::endl;
		exit(1);
	}
	return idx;
}

} /* namespace FermiOwn */

#endif /* SRC_DATA_GROSSNEVEUKFIELD_H_ */
