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

template <class FieldType> class MatrixChanges {
public:
	MatrixChanges( const DSlashUpdater& dslashUp, const FieldType& booleanField ) :
		dslash( dslashUp ),
		kfield( booleanField ),
		internalRanges( kfield.getInternalRanges() ),
		V( kfield.getVolume() )
	{}

	virtual ~MatrixChanges() {};

	/**
	 * @brief Calculate the difference of the field with the other for all lattice points
	 *
	 * The other field is assumed to be the older one, such that rows or cols
	 * that are not present there, but in the classes field, are returned in
	 * addrows/addcols.
	 *
	 * @param other is an older version of the current field to be compared
	 * @return a structure with lists of indices for rows and columns to be
	 * 		   added or deleted from the Dirac operator.
	 */
	const AddDelRowCol& calculateDifference( const FieldType& other ) {
		idxVec range( V );
		std::iota( range.begin(), range.end(), 0 );
		std::set<size_t> changedAt( range.begin(), range.end() );
		return calculateDifference( changedAt, other );
	}

	/**
	 * @brief Calculate the difference at a given number of points to the other field.
	 *
	 * @param changedAt is a set of lattice points x, where the difference to other is calculated.
	 * @param other is an older version of the current field to be compared
	 * @return a structure with lists of indices for rows and columns to be
	 * 		   added or deleted from the Dirac operator.
	 */
	const AddDelRowCol& calculateDifference( const std::set<size_t>& changedAt, const FieldType& other );

	/**
	 * @brief Prints out the currently stored indices to update a Dirac matrix.
	 */
	void Print() const {
		std::cout << "Current indices to change: " << std::endl << "delCols:";
		for( size_t col : indices.delCols ) std::cout << col << " ";
		std::cout << std::endl;
		std::cout << "delRows:";
		for( size_t row : indices.delRows ) std::cout << row << " ";
		std::cout << std::endl;
		std::cout << "addCols:";
		for( size_t col : indices.addCols ) std::cout << col << " ";
		std::cout << std::endl;
		std::cout << "addRows:";
		for( size_t row : indices.addRows ) std::cout << row << " ";
		std::cout << std::endl;
	};

	AddDelRowCol indices;

private:
	const DSlashUpdater& dslash;
	const FieldType& kfield;
	const idxVec& internalRanges;
	const size_t V;
};


template class MatrixChanges<ThirringKField>;

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_MATRIXCHANGES_H_ */
