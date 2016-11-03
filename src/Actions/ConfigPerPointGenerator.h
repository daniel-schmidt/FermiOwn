/*
 * ConfigGenerator.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_CONFIGPERPOINTGENERATOR_H_
#define SRC_ACTIONS_CONFIGPERPOINTGENERATOR_H_

#include <iostream>
#include <Eigen/Dense>
#include "Constants.h"

namespace FermiOwn {

/**
 * @brief Creates a list of all allowed configurations at a single lattice point.
 *
 * Since all our constraints are independent of x, the whole lattice configuration
 * can be obtained by selecting one of the allowed configurations at each lattice point.
 * The only input needed are the number of spins and flavours.
 * This class does not implement the constraints on its own, but uses the methods provided
 * by the FieldBoolean.
 */
class ConfigPerPointGenerator {
public:
	/**
	 * @brief Constructor setting physical parameters.
	 *
	 * This does not create configs!
	 *
	 * @param numSpins is the number of spinor components of the model, only the value 2 is tested at the moment.
	 * @param numFlavours is the number of fermion flavours.
	 * @param randomGenerator is necessary to enable drawing random configurations
	 */
	inline ConfigPerPointGenerator( const size_t numSpins, const size_t numFlavours, std::ranlux48* randomGenerator = NULL );

	/**
	 * @brief Default destructor, doing nothing.
	 */
	inline virtual ~ConfigPerPointGenerator();

	/**
	 * @brief Start computation of all allowed configurations.
	 *
	 * This depends on the model and must be implemented in derived classes.
	 */
	virtual void generateAllowedConfs()=0;

	/**
	 * @brief Get a list of all allowed configurations
	 *
	 * @return A matrix, where each row represents an allowed configuration.
	 */
	inline MatrixXb getAllConfs();

	/**
	 * @brief Returns a random, allowed configuration
	 *
	 * This needs a random generator set in the constructor.
	 *
	 * @return a row vector of bools, representing a random, allowed configuration at a single lattice point
	 */
	inline RowVectorXb getRandomConf();

protected:
	MatrixXb allowedConfs;		///< matrix containing an allowed config in each row
	const size_t Nf;			///< the number of flavours
	const size_t dimSpinor;		///< the number of spinor components, should be always 2, since we work in the irreducible representation

	std::ranlux48* rndGen;		///< a pointer to the random number generator ( or NULL if none was passed in the constructor )
	std::uniform_int_distribution<int> int_dist;	///< a random integer distribution ranging from 0 to allowedConfs, for drawing random configurations

};

ConfigPerPointGenerator::ConfigPerPointGenerator( const size_t numSpins, const size_t numFlavours, std::ranlux48* randomGenerator ) :
	Nf( numFlavours ),
	dimSpinor( numSpins ),
	rndGen( randomGenerator )
{
	if( dimSpinor != 2 ) std::cout << "Warning: ConfigGenerator called with spinor size != 2, this may not be implemented correctly!" << std::endl;
}

ConfigPerPointGenerator::~ConfigPerPointGenerator() {}

inline MatrixXb ConfigPerPointGenerator::getAllConfs() {
	return allowedConfs;
}

inline RowVectorXb ConfigPerPointGenerator::getRandomConf() {
	if( rndGen == NULL ) {
		std::cerr << "ConfigGenerator has no random number generator! Cannot draw random configuration." << std::endl;
		exit(1);
	}
	int conf = int_dist( *rndGen );
	return allowedConfs.row( conf );
}

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_CONFIGPERPOINTGENERATOR_H_ */
