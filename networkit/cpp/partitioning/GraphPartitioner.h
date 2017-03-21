/*
 * GraphPartitioner.h
 *
 *  Created on: 21.03.2017
 *      Author: moritzl
 */

#ifndef GRAPHPARTITIONER_H_
#define GRAPHPARTITIONER_H_

#include <vector>
#include <utility>

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"

namespace NetworKit {

class GraphPartitioner : public Algorithm {
public:

	GraphPartitioner(const Graph& G, count numParts, double maxImbalance, std::vector<index> chargedVertices);

	virtual void run() = 0;

	/**
	 * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
	 * @return partition of the node set
	 */
	virtual Partition getPartition() const;

	static edgeweight fiducciaMattheysesStep(const Graph& G, Partition& input, double maxImbalance, const std::vector<index> chargedVertices = {}, std::vector<double> nodeWeights = {});

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

	static double getWeightedImbalance(const Partition&  part, const std::vector<double> &nodeWeights, count k) {

		assert(part.numberOfElements() == nodeWeights.size());

		std::vector<double> fragmentSizes(part.upperBound());
		double maxFragmentSize = 0;
		double sumNodeWeights = 0;
		double largestNodeWeight = 0;

		for (index v = 0; v < part.numberOfElements(); v++) {
			const double weight = nodeWeights[v];

			assert(weight >= 0);

			sumNodeWeights += weight;
			fragmentSizes[part[v]] += weight;

			if (weight > largestNodeWeight) {
				largestNodeWeight = weight;
			}

			if (fragmentSizes[part[v]] < maxFragmentSize) {
				maxFragmentSize = fragmentSizes[part[v]];
			}
		}
		return getWeightedImbalance(fragmentSizes, sumNodeWeights, largestNodeWeight, part.numberOfElements(), k);
	}


protected:
	const Graph& G;
	const count numParts;
	const double maxImbalance;
	const std::vector<index> chargedNodes;
	Partition result;
	Partition previousPartition;
};

} /* namespace NetworKit */
#endif /* GRAPHPARTITIONER_H_ */
