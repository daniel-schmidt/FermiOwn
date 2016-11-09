/*
 * FermiBoolMetropolis.cpp
 *
 *  Created on: 31.03.2016
 *      Author: dschmidt
 */

#include "FermiBoolMetropolis.h"

namespace FermiOwn {

FermiBoolMetropolis::FermiBoolMetropolis( BasicKField& boolField, const Lattice & lattice, double lambda, size_t numFlavours, std::ranlux48* randomGenerator ) :
		MetropolisStep( randomGenerator ),
		kappa( 2./lambda ),
		Nf( numFlavours ),
		kfield( boolField ),
//		oldField(kfield),
		lat(lattice),
//		uni_real_dist(),
		intV_dist( 0, lat.getVol()-1 ),
//		int2mu_dist( 0, 2*lat.getDim()-1 ),
//		intNf_dist( 0, numFlavours-1 ),
//		intSpin_dist( 0, 1 ),
//		confGen( 2, numFlavours, randomGenerator ),
		weightChange(0.,0.),
		phase(0.),
		expPhase(0.,0.),
		fWeight( "weight" + std::to_string(lambda) + ".dat" )
{
	oldField = kfield.clone();
	try {
		ThirringKField& thfield = dynamic_cast<ThirringKField&>( kfield );
		weightFun = new ThirringWeightFunction( thfield, lat.getTimeSize(), lat.getSpaceSize(), lat.getDim(), numFlavours, lambda );
		confGen = new ConfigPerPointGeneratorTh( 2, numFlavours, randomGenerator );
	} catch( const std::bad_cast& e ) {
		std::cerr << "HARHAR" << std::endl;
	}
	confGen->generateAllowedConfs();
	//	generateAllowedConfs( numFlavours );
	//	intConfIndex = std::uniform_int_distribution<int>( 0, allowedConfs.rows()-1 );
	//
	//	slac.erase( kxiab );
	//	detOld = slac.det();
	//	slac.setFull();
	//	std::string filename = "weight";
	//	filename = filename + std::to_string(lambda) + ".dat";
}

FermiBoolMetropolis::~FermiBoolMetropolis() {
	delete oldField;
	delete weightFun;
	delete confGen;
}

void FermiBoolMetropolis::initializeField() {
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		RowVectorXb newConf = confGen->getRandomConf();
		kfield.setRow( newConf, x );
	}
}

//bool FermiBoolMetropolis::updateField() {
//
//	oldField = kxiab;
//
//	weightFun.saveState();
//
//	// draw random points, which get a new configuration
//	std::set<size_t> changedPoints;
//	for( int cnt=0; cnt < 2; cnt++ ) {
//		int x = intV_dist(*rndGen);
//		RowVectorXb newConf = confGen.getRandomConf();
//		kxiab.setRow( newConf, x );
//		changedPoints.insert( size_t( x ) );
//	}
//
//	// changing two neighbouring points
////	int x = intV_dist(*rndGen);
////	kxiab.setRow( confGen.getRandomConf(), x );
////	changedPoints.insert( size_t( x ) );
////
////	int mu = int2mu_dist( *rndGen );
////	x = lat.getNeighbours( x )[mu];
////	kxiab.setRow( confGen.getRandomConf(), x );
////	changedPoints.insert( size_t( x ) );
//
//	Complex weight = weightFun.updateWeight( changedPoints );
//
//	bool accepted = accept( weight );
//
//	return accepted;
//}

//bool FermiBoolMetropolis::accept( Complex weight ) {
//	double r = uni_real_dist(*rndGen);
//	bool accepted = false;
//	if( std::fabs(weight) > r ) {
//		phase += std::arg( weight );
//		phase = fmod( phase, 2*PI );
//		accepted = true;
//		acceptanceCounter++;
//		weightFun.keep();
//	} else {
//		kxiab = oldField;
//		weightFun.reset();
//	}
//	expPhase += std::exp( I*phase );
//	fWeight << std::real( weight ) << "\t" << std::imag( weight ) << "\t" << phase << "\t" << std::real(std::exp( I*phase )) << "\t" << std::imag(std::exp( I*phase ));
//	fWeight << "\t" << accepted << std::endl;
//	return accepted;
//}


//void FermiBoolMetropolis::updateNaive( size_t x ) {
//	if( intSpin_dist(*rndGen) ) {
//
//		int spin = intSpin_dist(*rndGen);
//		int a = intNf_dist(*rndGen);
//		int b = intNf_dist(*rndGen);
//
//		kxiab.invert( x, spin, a, b );
//	}
//}

//bool FermiBoolMetropolis::updateField( size_t x, size_t spin, size_t a, size_t b ) {
//
//	kxiab.invert( x, spin, a, b );
//	kxiab.enforceConstraint( x, spin, a, b );
//
////	std::cout << std::endl << "First change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
//	if( a == b ) {
////		 updating with same flavour,
////		 choose, if we update another spin (un-/setting a kxaa=2) or another flavour (setting two kxaa=1, keeping n1 even)
//		size_t c = intNf_dist(*rndGen);
//		size_t spin2 = 0;
//		if( a==c ) {
//			// same flavour drawn again, we have to update the second spin (otherwise we would revert our initial update)
//			spin2 = 1-spin;
////			std::cout << "New c drawn." << std::endl;
//		} else {
//			// different flavour, spin can be chosen randomly
//			spin2 = intSpin_dist(*rndGen);
////			std::cout << "New spin and b drawn." << std::endl;
//		}
//		kxiab.invert( x, spin2, c, c );
//		kxiab.enforceConstraint( x, spin2, c, c );
//	} else {
//		// updating with two different flavours, enforce kxab = kxba
//		spin = intSpin_dist(*rndGen);
////		std::cout << "New spin drawn." << std::endl;
//		kxiab.invert( x, spin, b, a );
//		kxiab.enforceConstraint( x, spin, b, a );
//	}
////	std::cout << "Second change at x=" << x << " spin=" << spin << " a=" << a << " b=" << b << std::endl;
////	kxiab.Print();
//	bool success = true;
//	if( kxiab.constraintViolated( x ) ) {
////		std::cout << "x violated! Reset and continue loop..." << std::endl;
//		kxiab = oldField;
//		success = false;
//	}
//
//	return success;
//}

//Complex FermiBoolMetropolis::calculateWeight( int dk, int dntilde ) {
//				std::cout << "Matrix before erase: " << std::endl << slac.getMatrix() << std::endl << std::endl;
//	slac.erase( kxiab );
//	//				std::cout << slac.getMatrix() << std::endl << std::endl;
//	det = slac.det();
//	slac.setFull();
//
//	double factor = getHypergeometricFactor( nxNew(1), nxNew(2) )  * getHypergeometricFactor( nyNew(1), nyNew(2) );
//	factor /= getHypergeometricFactor( nxOld(1), nxOld(2) ) * getHypergeometricFactor( nyOld(1), nyOld(2) );
//
//	//					std::cout << "\tfactor1=" << factor;
//
//	//	if( dntilde != 0 ) std::cerr << "\tdntilde = " << dntilde << std::endl;
//	factor *= std::pow( 2., dntilde );
//	Complex dw = std::pow( Complex(-kappa), double(dk)/2. );
//	Complex detRatio = det/detOld;
//	fWeight << std::real(detRatio) << "\t" << std::imag(detRatio) << "\t" << dw << "\t" << std::real(factor*dw*(detRatio)) << "\t" << std::imag(factor*dw*(detRatio)) << std::endl;
////	std::cout << "reDet=" << std::real(detRatio) << "\t imDet=" << std::imag(detRatio) << "\t dk=" << dk << "\t kappa=" << kappa << "\t dw=" << dw << "\t factor=" << factor << "\t dnt=" << dntilde << "\t" << std::fabs(factor*dw*(detRatio)) << std::endl;
//	return factor*dw*(det/detOld);
//}

void FermiBoolMetropolis::propose() {
	delete oldField;
	oldField = kfield.clone();

	weightFun->saveState();

	// draw random points, which get a new configuration
	changedPoints.clear();
	for( int cnt=0; cnt < 2; cnt++ ) {
		int x = intV_dist(*rnd);
		RowVectorXb newConf = confGen->getRandomConf();
		kfield.setRow( newConf, x );
		changedPoints.insert( size_t( x ) );
	}
}

double FermiBoolMetropolis::change() {
	weightChange = weightFun->updateWeight( changedPoints );

	return std::abs( weightChange );
}

void FermiBoolMetropolis::accept() {
	weightFun->keep();

	phase += std::arg( weightChange );
	phase = fmod( phase, 2*PI );
	writeWeightFile();
}

void FermiBoolMetropolis::reject() {
	kfield = *oldField->clone();
	weightFun->reset();
	writeWeightFile();
}

void FermiBoolMetropolis::writeWeightFile() {
	Complex currExpPhase = std::exp( I*phase );
	expPhase += currExpPhase;
	fWeight << std::real( weightChange ) << "\t" << std::imag( weightChange )
			<< "\t" << phase << "\t" << std::real( currExpPhase ) << "\t" << std::imag( currExpPhase ) <<std::endl;
}

Complex FermiBoolMetropolis::getAveragePhase() {
	return expPhase/double(getStepCount());
}

} /* namespace FermiOwn */
