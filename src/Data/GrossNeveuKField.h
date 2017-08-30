/*
 * GrossNeveuKField.h
 *
 *  Created on: 03.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_GROSSNEVEUKFIELD_H_
#define SRC_DATA_GROSSNEVEUKFIELD_H_

#include "BasicCloneKField.h"

namespace FermiOwn {

class GrossNeveuKField: public BasicCloneKField<GrossNeveuKField> {
public:
	GrossNeveuKField( const size_t latticeVolume, const idxVec& internal ) :
			BasicCloneKField( latticeVolume, internal )
	{};
	virtual ~GrossNeveuKField() {};

	const idxMap tally() const;
	const idxMap tally( const idxSet& xValues ) const;
private:
	virtual inline size_t colIndex( const std::vector<size_t>& internal ) const;
};

/*
 * Private functions
 ********************************************************/

inline size_t GrossNeveuKField::colIndex( const idxVec& internal ) const {
	size_t idx = internalRanges[0] * internal[1] + internal[0];
	if( idx >= DoFperX ) {
		std::cerr << "Index " << idx << " of bound (" << DoFperX << ") in GrossNeveuKField for spin=" << internal[0] << " and flavour=" << internal[1] << std::endl;
		exit(1);
	}
	return idx;
}
} /* namespace FermiOwn */

#endif /* SRC_DATA_GROSSNEVEUKFIELD_H_ */
