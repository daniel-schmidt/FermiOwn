/*
 * Constants.h
 *
 *  Created on: 14.03.2016
 *      Author: dschmidt
 */

#ifndef SRC_DATA_CONSTANTS_H_
#define SRC_DATA_CONSTANTS_H_

#include <complex>

namespace FermiOwn {

	typedef std::complex< double > Complex;
	typedef double Real;

	const double PI = 3.14159265358979;
	const Complex I = Complex(0,1);
} // namespace FermiOwn



#endif /* SRC_DATA_CONSTANTS_H_ */
