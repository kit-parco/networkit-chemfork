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

class MultiLevelPartitioner : public Algorithm {
	friend class PartitionerGTest;

public:
	MultiLevelPartitioner(const Graph& G, count numParts = 10, double maxImbalance = 2, bool bisectRecursivelyForInitialPartitioning = false, const std::vector<index>& chargedVertices = {}, bool avoidSurroundedNodes = false, Partition previous = Partition(0));


	virtual ~MultiLevelPartitioner() = default;

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
	static edgeweight fiducciaMattheysesStep(const Graph& G, Partition& input, double maxImbalance, const std::vector<index> chargedVertices = {}, std::vector<double> nodeWeights = {});
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints);
	static Partition growRegions(const Graph& g, const std::vector<index>& startingPoints, const std::vector<count>& weights, const Partition& constraint);

protected:
	static void enforceBalance(const Graph& G, Partition& part, double maxImbalance, const std::vector<index>& chargedVertices, const std::vector<double>& nodeWeights);
	static Partition partitionRecursively(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const std::vector<index>& chargedVertices, const Partition& previous, const std::vector<double> &nodeWeights);
	static Partition recursiveBisection(const Graph& g, count k);
	static void recursiveBisection(const Graph& g, count k, Partition& input, index maskID);
	static void repairSingleNodes(const Graph& g, Partition& intermediate);

	static index getFarthestNode(const Graph& g, std::vector<index> seedNodes);

	static std::pair<index, index> getMaximumDistancePair(const Graph& g, const Partition& constraint, const index partition);

	static bool chargesValid(const Partition& part, const std::vector<index>& chargedVertices) {
		for (node u : chargedVertices) {
			for (node v : chargedVertices) {
				if (part[u] == part[v] && u != v) {
					return false;
				}
			}
		}
		return true;
	}

	static double getWeightedImbalance(const std::vector<double>& fragmentSizes, double sumNodeWeights, double largestNodeWeight, double n, double k) {
		const double optSize = largestNodeWeight <= 1 ? ceil(double(n) / k) : sumNodeWeights / k + largestNodeWeight;
		return ((*std::max_element(fragmentSizes.begin(), fragmentSizes.end())) / optSize) - 1;
	}

	const Graph& G;
	const count numParts;
	const double maxImbalance;
	const bool bisectRecursively;
	const std::vector<index> chargedNodes;
	const bool noSingles;
	Partition result;
	Partition previousPartition;

};

} /* namespace NetworKit */
#endif /* PARTITIONER_H_ */
