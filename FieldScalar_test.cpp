/*
 * FieldScalar_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include "FieldScalar.h"
#include "Lattice.h"

int main() {
	size_t Nt = 3, Ns = 2, dim = 2;
	Lattice lat(Nt, Ns, dim);

	FieldScalar fs0(lat, zeroInit);
	fs0.Print();
	FieldScalar fs1(lat, oneInit);
	fs1.Print();
	FieldScalar fsrnd(lat, randomInit);
    fsrnd.Print();

}
