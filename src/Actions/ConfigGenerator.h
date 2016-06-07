/*
 * ConfigGenerator.h
 *
 *  Created on: 07.06.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_CONFIGGENERATOR_H_
#define SRC_ACTIONS_CONFIGGENERATOR_H_

namespace FermiOwn {

class ConfigGenerator {
public:
	ConfigGenerator( size_t latticeVolume, size_t numSpins, size_t numFlavours );
	virtual ~ConfigGenerator();

	void generateAllowedConfs();

	inline Eigen::MatrixXi getAllConfs();

private:
	Eigen::MatrixXi allowedConfs;	//TODO: should have a boolean type...
	size_t Nf;			///< the number of flavours
	size_t dimSpinor;	///< the number of spinor components, should be always 2, since we work in the irreducible representation
	size_t volume;		///< the lattice volume
};

inline Eigen::MatrixXi ConfigGenerator::getAllConfs() {
	return allowedConfs;
}

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_CONFIGGENERATOR_H_ */
