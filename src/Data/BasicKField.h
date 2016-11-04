/*
 * GrossNeveuKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_BASICKFIELD_H_
#define SRC_DATA_BASICKFIELD_H_

#include <iostream>
#include <vector>
#include "Constants.h"

namespace FermiOwn {

class BasicKField {
public:
	inline BasicKField( const size_t latticeVolume, const std::vector<size_t>& internal );
	inline virtual ~BasicKField();

	virtual BasicKField* clone() const =0;

	inline void setValue( bool val, size_t x, const std::vector<size_t>& internal );
	inline bool getValue( size_t x, const std::vector<size_t>& internal ) const;
	inline void setRow( const RowVectorXb & newRow, size_t x );
	inline void setZero();
	inline void invert( size_t x, const std::vector<size_t>& internal );

	inline void Print() const;

protected:
	virtual size_t colIndex( const std::vector<size_t>& internal ) const =0;

	size_t V;
	std::vector<size_t> internalRanges;		//TODO: see, if we can make this const by providing a copy operator for the class
	size_t DoFperX;

	MatrixXb data; //TODO: make this private again.
private:
};

/*
 * Implementation of functions
 **********************************************************/

BasicKField::BasicKField( const size_t latticeVolume, const std::vector<size_t>& internal ) :
			V( latticeVolume ),
			internalRanges( internal ),
			DoFperX(1)
{
	for( size_t internalN : internalRanges ) {
		DoFperX *= internalN;
	}
	data = MatrixXb::Zero( V, DoFperX );
}

BasicKField::~BasicKField() {
}

void BasicKField::setValue( bool val, size_t x, const std::vector<size_t>& internal ) {
	//TODO: define preprocessor variable to disable out-of-bounds checks
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::setValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( internal) ) = val;
}

inline bool BasicKField::getValue( size_t x, const std::vector<size_t>& internal ) const {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::getValue: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	return data( x, colIndex( internal ) );
}

inline void BasicKField::setRow( const RowVectorXb & newRow, size_t x ) {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::setRow: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	} else if( newRow.cols() > int( DoFperX ) ) {
		std::cerr << "Error in GrossNeveuKField::setRow: Row length " << newRow.cols() << " is longer than the number of degrees of freedom per site: " << DoFperX << std::endl;
		exit(1);
	}
	data.row(x) = newRow;
}

inline void BasicKField::setZero() {
	data = MatrixXb::Zero( V, DoFperX );
}

inline void BasicKField::invert( size_t x, const std::vector<size_t>& internal) {
	if( x >= V ) {
		std::cerr << "Error in GrossNeveuKField::invert: Lattice point " << x << " is outside the volume of " << V << std::endl;
		exit(1);
	}
	data( x, colIndex( internal ) ) = !data( x, colIndex( internal ) );
}

inline void BasicKField::Print() const {
	std::cout << data << std::endl;
}


} /* namespace FermiOwn */

#endif /* SRC_DATA_BASICKFIELD_H_ */
