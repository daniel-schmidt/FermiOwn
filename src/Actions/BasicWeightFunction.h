/*
 * BasicWeightFunction.h
 *
 *  Created on: 08.11.2016
 *      Author: dschmidt
 */

#ifndef SRC_ACTIONS_BASICWEIGHTFUNCTION_H_
#define SRC_ACTIONS_BASICWEIGHTFUNCTION_H_

namespace FermiOwn {

	class BasicWeightFunction {
	public:

		BasicWeightFunction(){};

		virtual ~BasicWeightFunction(){};
		/**
		 * @brief Calculates the full weight
		 * @return the value of the weight
		 */
		virtual Complex calculateWeight()=0;

		/**
		 * @brief Calculates the change in weight
		 * @return returns difference
		 */
		virtual Complex updateWeight( const std::set< size_t > & changedAt )=0;

		virtual void saveState() =0;

		virtual void keep() =0;

		virtual void reset() =0;
	};
} /* end namespace FermiOwn */



#endif /* SRC_ACTIONS_BASICWEIGHTFUNCTION_H_ */
