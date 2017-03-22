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
#include "../matching/LocalMaxMatcher.h"
#include "../coarsening/ClusteringProjector.h"
#include "../coarsening/MatchingContracter.h"

#include "GraphPartitioner.h"

namespace NetworKit {

class MultiLevelPartitioner : public GraphPartitioner {
	friend class PartitionerGTest;

public:
	MultiLevelPartitioner(const Graph& G, count numParts = 10, double maxImbalance = 2, bool bisectRecursivelyForInitialPartitioning = false, const std::vector<index>& chargedVertices = {}, count minimumGapSize = 0, Partition previous = Partition(0));


	virtual ~MultiLevelPartitioner() = default;

	/**
	 * Apply algorithm to graph
	 */
	virtual void run() override;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

	static edgeweight calculateGain(const Graph& g, const Partition& input, index v, index targetPart);
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints);
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints, const std::vector<count>& weights, const Partition& constraint);

protected:
	static Partition partitionRecursively(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const std::vector<index>& chargedVertices, const Partition& previous, const std::vector<double> &nodeWeights);
	static Partition recursiveBisection(const Graph& g, count k);
	static void recursiveBisection(const Graph& g, count k, Partition& input, index maskID);
	static void repairSingleNodes(const Graph& g, Partition& intermediate);

	static index getFarthestNode(const Graph& g, std::vector<index> seedNodes);

	static std::pair<index, index> getMaximumDistancePair(const Graph& g, const Partition& constraint, const index partition);

	const bool bisectRecursively;
	Partition previousPartition;
};

} /* namespace NetworKit */
#endif /* PARTITIONER_H_ */
