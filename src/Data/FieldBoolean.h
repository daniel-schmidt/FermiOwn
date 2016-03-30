/*
 * FieldBoolean.h
 *
 *  Created on: 30.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_FIELDBOOLEAN_H_
#define SRC_DATA_FIELDBOOLEAN_H_

#include "Field.h"

namespace FermiOwn {

class FieldBoolean: public Field<bool> {
public:
	FieldBoolean( const size_t latticeVolume, const size_t numberOfSpin, const size_t numberOfFlavours, std::ranlux48 * rndGen, InitType init );
	virtual ~FieldBoolean();
	bool getValue( size_t x, size_t spin, size_t flavour1, size_t flavour2 ) const;
	void setValue( bool val, size_t x, size_t spin, size_t flavour1, size_t flavour2 );
	void invert(size_t x, size_t spin, size_t flavour1, size_t flavour2);
	bool constraintViolated( size_t x) const;
	size_t sumAll() const;

private:
	size_t colIndex( size_t spin, size_t flavour1, size_t flavour2 ) const;

	size_t numSpins;
	size_t numFlavours;
};

} /* namespace FermiOwn */

#endif /* SRC_DATA_FIELDBOOLEAN_H_ */
