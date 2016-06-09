/*
 * ConfigGenerator.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_CONFIGGENERATOR_H_
#define SRC_ACTIONS_CONFIGGENERATOR_H_

#include <Eigen/Dense>
#include "FieldBoolean.h"

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
class ConfigGenerator {
public:
	/**
	 * @brief Constructor setting physical parameters.
	 *
	 * This does not create configs!
	 *
	 * @param numSpins is the number of spinor components of the model, only the value 2 is tested at the moment.
	 * @param numFlavours is the number of fermion flavours.
	 */
	ConfigGenerator( const size_t numSpins, const size_t numFlavours );

	/**
	 * @brief Default destructor, doing nothing.
	 */
	virtual ~ConfigGenerator();

	/**
	 * @brief Start computation of all allowed configurations.
	 *
	 * This needs currently a huge amount of memory for Nf>3 and fails.
	 */
	void generateAllowedConfs();

	/**
	 * @brief Get a list of all allowed configurations
	 *
	 * @return A matrix, where each row represents an allowed configuration.
	 */
	inline Eigen::MatrixXi getAllConfs();

private:
	Eigen::MatrixXi allowedConfs;	//TODO: should have a boolean type...
	const size_t Nf;			///< the number of flavours
	const size_t dimSpinor;	///< the number of spinor components, should be always 2, since we work in the irreducible representation
};

inline Eigen::MatrixXi ConfigGenerator::getAllConfs() {
	return allowedConfs;
}

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_CONFIGGENERATOR_H_ */
