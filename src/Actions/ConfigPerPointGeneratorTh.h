/*
 * ConfigPerPointGeneratorTh.h
 *
 *  Created on: 02.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_CONFIGPERPOINTGENERATORTH_H_
#define SRC_ACTIONS_CONFIGPERPOINTGENERATORTH_H_

#include "ConfigPerPointGenerator.h"
#include "ThirringKField.h"

namespace FermiOwn {

class ConfigPerPointGeneratorTh: public ConfigPerPointGenerator {
public:
	ConfigPerPointGeneratorTh( const size_t numSpins, const size_t numFlavours, std::ranlux48* randomGenerator = NULL );
	virtual ~ConfigPerPointGeneratorTh();

	/**
	 * @brief Start computation of all allowed configurations.
	 *
	 * This needs currently a huge amount of memory for Nf>3 and fails.
	 */
	void generateAllowedConfs();
};

} /* namespace FermiOwn */

#endif /* SRC_ACTIONS_CONFIGPERPOINTGENERATORTH_H_ */
