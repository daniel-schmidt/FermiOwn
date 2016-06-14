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

	// basic numerical types
	typedef std::complex< double > Complex;
	typedef double Real;

	// boolean fields from Eigen
	typedef Eigen::Matrix< bool, Eigen::Dynamic, 1 > VectorXb;
	typedef Eigen::Matrix< bool, 1, Eigen::Dynamic > RowVectorXb;
	typedef Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > MatrixXb;

	// things related to sparse matrices
	typedef std::vector< Eigen::Triplet<Complex> > MatCoeffList;
	typedef Eigen::SparseMatrix<Complex, Eigen::ColMajor> SparseMat;

	// explicit constants
	const double PI = 3.14159265358979;
	const Complex I = Complex(0.,1.);
	const double ZERO_TOL = 1e-10;
} // namespace FermiOwn



#endif /* SRC_DATA_CONSTANTS_H_ */
