/*
 * Partitioner.h
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#ifndef PARTITIONER_H_
#define PARTITIONER_H_

#include <vector>
#include <utility>

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include "../matching/LocalMaxMatcher.h"
#include "../coarsening/ClusteringProjector.h"
#include "../coarsening/MatchingContracter.h"

namespace NetworKit {

class Partitioner : public Algorithm {
	friend class PartitionerGTest;

public:
	Partitioner(const Graph& G, count numParts = 10, double maxImbalance = 10, bool bisectRecursivelyForInitialPartitioning = false, const std::vector<index>& chargedVertices = {});

	//Partitioner(const Graph& G, const std::vector<index>& chargedVertices, double maxImbalance = 10, bool bisectRecursivelyForInitialPartitioning = true);
	virtual ~Partitioner() = default;

	/**
	 * Apply algorithm to graph
	 */
	virtual void run() override;

	/**
	 * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
	 * @return partition of the node set
	 */
	virtual Partition getPartition();

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

	static edgeweight calculateGain(const Graph& g, const Partition& input, index v, index targetPart);
	static edgeweight fiducciaMatheysesStep(const Graph& G, Partition& input);
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints);
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints, const std::vector<count>& weights, const Partition& constraint);

protected:
	static Partition partitionRecursively(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const std::vector<index>& chargedVertices);
	static Partition recursiveBisection(const Graph& g, count k);
	static void recursiveBisection(const Graph& g, count k, Partition& input, index maskID);

	static std::pair<index, index> getMaximumDistancePair(const Graph& g, const Partition& constraint, const index partition);

	const Graph& G;
	const count numParts;
	const bool charged;
	const double maxImbalance;
	const bool bisectRecursively;
	const std::vector<index> chargedNodes;
	Partition result;

};

} /* namespace NetworKit */
#endif /* PARTITIONER_H_ */
