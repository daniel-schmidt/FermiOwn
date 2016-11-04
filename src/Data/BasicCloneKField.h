/*
 * BasicCloneKField.h
 *
 *  Created on: 04.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_BASICCLONEKFIELD_H_
#define SRC_DATA_BASICCLONEKFIELD_H_

#include "BasicKField.h"

namespace FermiOwn {

/**
 * Intermediate Class providing clone implementations for all derived classes via Curiously Recurring Template Pattern
 *
 * http://stackoverflow.com/questions/12255546/c-deep-copying-a-base-class-pointer
 */

template <typename Derived> class BasicCloneKField : public BasicKField {
public:
	BasicCloneKField( const size_t latticeVolume, const std::vector<size_t>& internal ) :
		BasicKField( latticeVolume, internal )
	{};
	virtual ~BasicCloneKField() {};

	virtual BasicKField* clone() const {
		return new Derived( static_cast<Derived const&>(*this) );
	}
};

} /* namespace FermiOwn */



#endif /* SRC_DATA_BASICCLONEKFIELD_H_ */
