/*
 * PlainGapPartitioner.h
 *
 *  Created on: 21.03.2017
 *      Author: moritzl
 */

#ifndef PLAINGAPPARTITIONER_H_
#define PLAINGAPPARTITIONER_H_

#include "GraphPartitioner.h"

namespace NetworKit {

class PlainGapPartitioner: public NetworKit::GraphPartitioner {
public:
	PlainGapPartitioner(const Graph& G, count numParts, double maxImbalance, std::vector<index> chargedVertices, count minGapSize = 3, Partition previous = Partition(0));
	virtual ~PlainGapPartitioner() = default;

	virtual void run() override;

protected:
	Partition previousPartition;

};

} /* namespace NetworKit */
#endif /* PLAINGAPPARTITIONER_H_ */
