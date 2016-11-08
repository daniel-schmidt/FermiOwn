/*
 * MatrixChanges.h
 *
 *  Created on: 08.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_MATRIXCHANGES_H_
#define SRC_ACTIONS_MATRIXCHANGES_H_

#include <set>
#include "ThirringKField.h"
#include "DSlashUpdater.h"

namespace FermiOwn {

class MatrixChanges {
public:
	MatrixChanges( const DSlashUpdater& dslashUp, const ThirringKField& booleanField );
	virtual ~MatrixChanges();

	const AddDelRowCol& calculateDifference( const ThirringKField& other );

	const AddDelRowCol& calculateDifference( const std::set<size_t>& changedAt, const ThirringKField& other );

	void Print() const;

	AddDelRowCol indices;
private:
	const DSlashUpdater& dslash;
	const ThirringKField& kfield;
	const idxVec& internalRanges;
	const size_t V;
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_MATRIXCHANGES_H_ */
