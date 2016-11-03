/*
 * FieldBoolean_test.cpp
 *
 *  Created on: 03.08.2016
 *      Author: dschmidt
 */

#include <random>
#include "ThirringKField.h"
#include "GrossNeveuKField.h"

int main() {
	using namespace FermiOwn;

	std::cout << std::endl << "Testing boolean field..." << std::endl;

	std::ranlux48 rndGen;

	ThirringKField fbool1( 2*3*3, {2, 3, 3} );
	fbool1.Print();
	ThirringKField fbool0( 2*3*3, {2, 3, 3} );
	fbool0.Print();
	std::cout << "Inverting, should be at row 3, col 6 and 10..." << std::endl;
	fbool1.setValue( false, 2, {1, 0, 1} );
	fbool1.setValue( false, 2, {0, 2, 1} );
	fbool1.Print();
	std::cout << "Inverting, should be at row 7, col 2..." << std::endl;
	fbool0.invert( 6, {0, 1, 1} );
	fbool0.Print();

//	std::cout << "Testing count, field 0 at x=6: " << std::endl << fbool0.countSummedSpin(6) << std::endl << " field 1 at x=3: " << std::endl << fbool1.countSummedSpin(3) << std::endl;
	std::cout << "Testing count offdiagonal, field 0: " << std::endl << fbool0.countOffdiagonal2() << std::endl << " field 1: " << std::endl << fbool1.countOffdiagonal2() << std::endl;

	std::cout << std::endl << "Testing copy assignment:" << std::endl;
	fbool1 = fbool0;
	fbool1.Print();

	std::cout << std::endl << "Testing change detection:" << std::endl;
	fbool1.invert( 7, {0, 1, 1} );
	fbool1.invert( 7, {1, 1, 1} );
	fbool1.invert( 8, {0, 1, 1} );
	std::cout << "Same values at " << " entries.";
	(fbool1.different(fbool0)).Print();
	std::cout << std::endl;
	fbool0.Print();
	std::cout << std::endl;
	fbool1.Print();

	std::cout << std::endl << "Testing setting whole rows: " << std::endl;
	RowVectorXb newRow(2*3*3);
	newRow << 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1;
	std::cout << newRow << std::endl << "This row should be found in the previous field in line 3 and 11:" << std::endl;
	fbool1.setRow( newRow, 2);
	fbool1.setRow( newRow, 10);

	fbool1.Print();

	std::cout << std::endl << "Testing GrossNeveu field:" << std::endl;

	std::vector<size_t> internalMax = {2, 2};
	GrossNeveuKField kxia( 2*3*3, internalMax );
	std::vector<size_t> internal = {1, 0};
	kxia.setValue( 1, 5, internal );
	RowVectorXb gnrow(4);
	gnrow << 0, 1, 1, 0;
	kxia.setRow( gnrow, 7 );
	internal = {0, 0};
	kxia.invert( 8, internal );
	kxia.invert( 8, {1, 1} );
	kxia.Print();
	internal = {1, 2};
	kxia.setValue( 1, 5, internal);
}
