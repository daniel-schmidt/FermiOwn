/*
 * GrossNeveuKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_GROSSNEVEUKFIELD_H_
#define SRC_DATA_GROSSNEVEUKFIELD_H_

#include <iostream>
#include <vector>

namespace FermiOwn {

class GrossNeveuKField {
public:
	GrossNeveuKField( const size_t latticeVolume, const std::vector<size_t>& internal );
	virtual ~GrossNeveuKField();

	inline void setValue( bool val, size_t x, const std::vector<size_t>& internal );
	inline bool getValue( size_t x, const std::vector<size_t>& internal ) const;
	inline void setRow( const RowVectorXb & newRow, size_t x );
	inline void invert( size_t x, const std::vector<size_t>& internal );

	inline void Print() const;

private:
	inline size_t colIndex( const std::vector<size_t>& internal ) const;

	size_t V;
	const std::vector<size_t> internalRanges;
	size_t DoFperX;

	MatrixXb data;
};

/*
 * Implementation of functions
 **********************************************************/

GrossNeveuKField::GrossNeveuKField( const size_t latticeVolume, const std::vector<size_t>& internal ) :
			V( latticeVolume ),
			internalRanges( internal ),
			DoFperX(1)
{
	for( size_t internalN : internalRanges ) {
		DoFperX *= internalN;
	}
	data = MatrixXb::Zero( V, DoFperX );
}

GrossNeveuKField::~GrossNeveuKField() {
}

void GrossNeveuKField::setValue( bool val, size_t x, const std::vector<size_t>& internal ) {
	//TODO: define preprocessor variable to disable out-of-bounds checks
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::setValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( internal) ) = val;
}

inline bool GrossNeveuKField::getValue( size_t x, const std::vector<size_t>& internal ) const {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::getValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	return data( x, colIndex( internal ) );
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

inline void GrossNeveuKField::invert( size_t x, const std::vector<size_t>& internal) {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::invert: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( internal ) ) = !data( x, colIndex( internal ) );
}

inline void GrossNeveuKField::Print() const {
	std::cout << data << std::endl;
}

/*
 * Private functions
 ********************************************************/

inline size_t GrossNeveuKField::colIndex( const std::vector<size_t>& internal ) const {
	size_t idx = internalRanges[0] * internal[1] + internal[0];
	if( idx >= DoFperX ) {
		std::cerr << "Index " << idx << " of bound (" << DoFperX << ") in GrossNeveuKField for spin=" << internal[0] << " and flavour=" << internal[1] << std::endl;
		exit(1);
	}
	return idx;
}

} /* namespace FermiOwn */

#endif /* SRC_DATA_GROSSNEVEUKFIELD_H_ */
