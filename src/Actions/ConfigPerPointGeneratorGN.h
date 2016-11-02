/*
 * ConfigPerPointGeneratorGN.h
 *
 *  Created on: 02.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_CONFIGPERPOINTGENERATORGN_H_
#define SRC_ACTIONS_CONFIGPERPOINTGENERATORGN_H_

#include <bitset>
#include <climits>
#include "ConfigPerPointGenerator.h"

namespace FermiOwn {

class ConfigPerPointGeneratorGN: public ConfigPerPointGenerator {
public:
	ConfigPerPointGeneratorGN( const size_t numSpins, const size_t numFlavours, std::ranlux48* randomGenerator = NULL );
	virtual ~ConfigPerPointGeneratorGN();

	void generateAllowedConfs();
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_CONFIGPERPOINTGENERATORGN_H_ */
