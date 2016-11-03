/*
 * GrossNeveuKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_GROSSNEVEUKFIELD_H_
#define SRC_DATA_GROSSNEVEUKFIELD_H_

#include <vector>
#include "BasicKField.h"

namespace FermiOwn {

class GrossNeveuKField: public BasicKField {
public:
	GrossNeveuKField( const size_t latticeVolume, const std::vector<size_t>& internal );
	virtual ~GrossNeveuKField();

protected:
	virtual inline size_t colIndex( const std::vector<size_t>& internal ) const;
};

GrossNeveuKField::GrossNeveuKField( const size_t latticeVolume, const std::vector<size_t>& internal ) :
		BasicKField( latticeVolume, internal )
{}

GrossNeveuKField::~GrossNeveuKField() {};
/*
 * Protected functions
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
