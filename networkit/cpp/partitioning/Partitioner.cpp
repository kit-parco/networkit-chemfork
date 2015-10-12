/*
 * Partitioner.cpp
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#include "Partitioner.h"

namespace NetworKit {
namespace GraphTools {

Partitioner::Partitioner(const Graph& G) : Algorithm(), G(G), result(0) {

}

void Partitioner::run() {
	const int numParts = 10;
	std::function<Partition(Graph&)> partitionLambda = [&](Graph& g) -> Partition {
	   count n = g.numberOfNodes();

	   // coarsen recursively until graph is small enough
	   if (n <= 2 * numParts) {
		   Partition initial(n);

		   // TODO: initial partitioning with BRKGA

		   return initial;
	   }
	   else {
		   // recursive coarsening
		   LocalMaxMatcher matcher(g);
		   Matching matching = matcher.run();
		   MatchingContracter coarsener(g, matching);
		   coarsener.run();
		   Graph coarseG = coarsener.getCoarseGraph();
		   std::vector<node> fineToCoarse = coarsener.getNodeMapping();

		   // recursive call
		   Partition coarsePart = partitionLambda(coarseG);

		   // interpolation
		   ClusteringProjector projector;
		   Partition finePart = projector.projectBack(coarseG, g, fineToCoarse, coarsePart);

		   // TODO: local refinement with BRKGA

		   return finePart;
	   }
   };
}
/**
DiscreteIndividual runMultilevelPartitioning =
   */

Partition Partitioner::getPartition() {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return result;
}

std::string Partitioner::toString() const {
	return "TODO";
}


} /* namespace GraphTools */
} /* namespace NetworKit */
