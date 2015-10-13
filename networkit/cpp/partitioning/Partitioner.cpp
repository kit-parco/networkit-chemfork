/*
 * Partitioner.cpp
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#include <tuple>

#include "Partitioner.h"
#include "../auxiliary/PrioQueue.h"

using std::pair;
using std::vector;
using Aux::PrioQueue;

namespace NetworKit {
namespace GraphTools {

Partitioner::Partitioner(const Graph& G) : Algorithm(), G(G), result(0) {

}

void Partitioner::run() {
	const int numParts = 10;

	std::function<Partition(const Graph&)> partitionLambda = [&](const Graph& g) -> Partition {
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

		   fiducciaMatheyses(g, finePart);

		   // TODO: local refinement with BRKGA

		   return finePart;
	   }
   };
	result = partitionLambda(G);
	hasRun = true;
}

void Partitioner::fiducciaMatheyses(const Graph& g, Partition&  part) {
	/**
	 * necessary data structures:
	 * - integer priority queue for each partition
	 * - history of move operations and gain values
	 */

	/**
	 * allocate data structures
	 */
	const count k = part.numberOfSubsets();
	const count n = part.numberOfElements();
	vector<PrioQueue<edgeweight, index> > queues(k,0);
	vector<edgeweight> gains;
	vector<pair<index, index> > transfers;
	vector<index> transferedVertices;
	edgeweight total = g.totalEdgeWeight();

	/**
	 * fill priority queues for the first time
	 */
	for (index ID : part.getSubsetIds()) {
		for (index j : part.getMembers(ID)) {
			edgeweight maxGain = -total;
			index IdAtMax = 0;
			for (index otherSubset : part.getSubsetIds()) {
				edgeweight thisgain = calculateGain(g, part, j, otherSubset);
				if (thisgain > maxGain && otherSubset != ID) {
					IdAtMax = otherSubset;
					maxGain = thisgain;
				}
			}
			assert(ID < queues.size());
			queues[ID].insert(maxGain, j);
		}
	}

	/**
	 * do one step
	 */
	bool allQueuesEmpty = false;
	while (!allQueuesEmpty) {
		allQueuesEmpty = true;
		for (index partID : part.getSubsetIds()) {
			assert(partID < queues.size());
			index topVertex;
			double topGain;
			std::tie(topGain, topVertex) = queues[partID].extractMin();

			//now get target partition. This could be sped up by saving the target partition in a separate data structure
			edgeweight maxGain = -total;
			index IdAtMax = 0;
			for (index otherSubset : part.getSubsetIds()) {
				edgeweight thisgain = calculateGain(g, part, topVertex, otherSubset);
				if (thisgain > maxGain && otherSubset != partID) {
					IdAtMax = otherSubset;
					maxGain = thisgain;
				}
			}

			//move node there
			part.moveToSubset(IdAtMax, topVertex);

			//update history
			gains.push_back(maxGain);
			transfers.emplace_back(partID, IdAtMax);
			transferedVertices.push_back(topVertex);

			//update gains of neighbours
			g.forNeighborsOf(topVertex, [&g, &topVertex, &partID, &total, &queues, &part](index w){
				//update gain
				edgeweight newMaxGain = -total;
				for (index otherSubset : part.getSubsetIds()) {
					edgeweight thisgain = calculateGain(g, part, topVertex, otherSubset);
					if (thisgain > newMaxGain && otherSubset != partID) {
						newMaxGain = thisgain;
					}
				}

				//update prioqueue
				queues[part[w]].remove(w);
				queues[part[w]].insert(newMaxGain, w);
			});

			allQueuesEmpty = allQueuesEmpty && queues[partID].size() == 0;
		}
	}
}

edgeweight Partitioner::calculateGain(const Graph& g, const Partition& input, index v, index targetPart) {
	return 0;
}

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
