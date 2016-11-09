/*
 * GrossNeveuKField.cpp
 *
 *  Created on: 09.11.2016
 *      Author: dschmidt
 */

#include "GrossNeveuKField.h"

namespace FermiOwn {

const idxMap GrossNeveuKField::tally() const {
	idxMap map;
	for( size_t x = 0; x < V; x++ ) {
		int kx = data.row( x ).count();
		map[kx]++;
	}
	return map;
}
} /* end namespace FermiOwn */
