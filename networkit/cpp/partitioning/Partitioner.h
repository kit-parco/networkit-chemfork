/*
 * Partitioner.h
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#ifndef PARTITIONER_H_
#define PARTITIONER_H_

#include "../graph/Graph.h"
#include "../structures/Partition.h"
#include "../base/Algorithm.h"
#include "../matching/LocalMaxMatcher.h"
#include "../coarsening/ClusteringProjector.h"
#include "../coarsening/MatchingContracter.h"

namespace NetworKit {
namespace GraphTools {

class Partitioner : public Algorithm {
public:
	Partitioner(const Graph& G);
	virtual ~Partitioner() = default;

	/**
	 * Apply algorithm to graph
	 */
	virtual void run() = 0;

	/**
	 * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
	 * @return partition of the node set
	 */
	virtual Partition getPartition();

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;

protected:
	void fiducciaMatheyses(const Graph& G, Partition& input);
	static edgeweight calculateGain(const Graph& g, const Partition& input, index v, index targetPart);

	const Graph& G;
	Partition result;

};

} /* namespace GraphTools */
} /* namespace NetworKit */
#endif /* PARTITIONER_H_ */
