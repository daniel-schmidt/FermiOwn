/*
 * FieldScalar_test.cpp
 *
 *  Created on: 03.02.2016
 *      Author: dschmidt
 */

#include "FieldScalar.h"
#include "Lattice.h"

int main() {
	size_t Nt = 2, Ns = 2, dim = 3;
	Lattice lat(Nt, Ns, dim);

	FieldScalar<Real> fs0(lat, zeroInit);
	fs0.Print();
	FieldScalar<Complex> fs1(lat, oneInit);
	fs1.Print();
	FieldScalar<Complex> fsrnd(lat, randomInit);
    fsrnd.Print();

    std::cout << std::endl << "Trying +=" << std::endl;
    fs1 += fsrnd;
    fs1.Print();

    std::cout << std::endl << "Trying +" << std::endl;
    (fs1 + fsrnd).Print();

    FieldScalar<Complex> fsmult(lat, oneInit);
    std::cout << std::endl << "Trying *=Scalar" << std::endl;
    fsmult *= Complex(0.,1.);
    fsmult.Print();
    fsmult /= Complex(0.,1.);
    fsmult.Print();

    FieldScalar<Complex> fsmult2(lat, oneInit);
    fsmult2 *= 2.;
    fsmult *= 3.;
    std::cout<< "Result of scalar product is " << (fsmult*fsmult2) << " and should be (48,0)."<< std::endl;

}
