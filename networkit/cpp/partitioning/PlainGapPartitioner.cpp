/*
 * PlainGapPartitioner.cpp
 *
 *  Created on: 21.03.2017
 *      Author: moritzl
 */

#include "PlainGapPartitioner.h"

namespace NetworKit {

PlainGapPartitioner::PlainGapPartitioner(const Graph& G, count numParts, double maxImbalance, std::vector<index> chargedVertices, count minGapSize, Partition previous) : GraphPartitioner(G, numParts, maxImbalance, chargedVertices, minGapSize), previousPartition(previous){

}

void PlainGapPartitioner::run() {
	//first, create the coarse graph. Run fiducciaMattheysesStep on the coarse graph. Uncoarsen. Run again, but with enforced constraints

}

} /* namespace NetworKit */
