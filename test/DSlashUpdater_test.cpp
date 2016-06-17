/*
 * DSlashUpdater_test.cpp
 *
 *  Created on: 16.06.2016
 *      Author: dschmidt
 */

#include "FieldBoolean.h"
#include "DSlashUpdater.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;
	FieldBoolean finit( 4, 2, 2, NULL, zeroInit );
	FieldBoolean ffinal( 4, 2, 2, NULL, zeroInit );
	DSlashUpdater dslash( 4, 1, 3, 2 );

	RowVectorXb conf80( 8 );
	conf80 << 0, 0, 1, 0, 0, 1, 0, 0;
	RowVectorXb conf96( 8 );
	conf96 << 0, 1, 0, 0, 0, 0, 1, 0;
//	finit.setRow( conf80, 0 );
//	finit.setRow( conf80, 1 );

	ffinal.setRow( conf96, 0 );
	ffinal.setRow( conf96, 1 );

	FieldBoolean diff = ffinal.different( finit );

	dslash.calculateUpdateMatrices( ffinal, diff );
}
