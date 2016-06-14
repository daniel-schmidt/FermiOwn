/*
 * WoodburyMatrix.h
 *
 *  Created on: 14.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_MATH_WOODBURYMATRIX_H_
#define SRC_MATH_WOODBURYMATRIX_H_

#include <iostream>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include "Constants.h"

namespace FermiOwn {

typedef std::vector< Eigen::Triplet<Complex> > MatCoeffList;
typedef Eigen::SparseMatrix<Complex, Eigen::ColMajor> SparseMat;
class WoodburyMatrix {
public:
	WoodburyMatrix( const size_t matSize, const MatCoeffList & list );
	virtual ~WoodburyMatrix();

	inline Complex getDet() const;

	void updateMatrix();
	Complex updateDet();
	void updateInverse();

	inline void Print() const;
	inline void PrintInverse() const;



private:
	SparseMat mat;
	const size_t size;

	SparseMat inv;
	Complex det;

};

/*============================================================
 * Inline Functions
 *============================================================*/

inline void WoodburyMatrix::Print() const {
	std::cout << mat << std::endl;
}

inline void WoodburyMatrix::PrintInverse() const {
	std::cout << inv << std::endl;
}

inline Complex WoodburyMatrix::getDet() const {
	return det;
}

} /* namespace FermiOwn */

#endif /* SRC_MATH_WOODBURYMATRIX_H_ */
