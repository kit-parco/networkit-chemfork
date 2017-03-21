/*
 * GraphPartitioner.cpp
 *
 *  Created on: 21.03.2017
 *      Author: moritzl
 */

#include "GraphPartitioner.h"

#include <tuple>
#include <queue>
#include <algorithm>

#include "../auxiliary/PrioQueue.h"

using std::vector;
using std::pair;
using std::queue;
using Aux::PrioQueue;

namespace NetworKit {

GraphPartitioner::GraphPartitioner(const Graph& G, count numParts, double maxImbalance, std::vector<index> chargedVertices) : G(G), numParts(numParts), maxImbalance(maxImbalance), chargedNodes(chargedVertices), result(0) {

}


Partition GraphPartitioner::getPartition() const {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return result;
}

//TODO: the node weights are currently copied, this could be improved. I could use a reference, but then I cannot make it const, since I resize the vector if it is empty. I allow empty nodeWeights for overloading.
edgeweight GraphPartitioner::fiducciaMattheysesStep(const Graph& g, Partition&  part, double maxImbalance, const std::vector<index> chargedVertices, std::vector<double> nodeWeights) {
	using std::vector;

	/**
	 * magic numbers
	 */
	const count KarypisKumarStoppingCriterion = 20;//currently unused
	count movesWithoutImprovement = 0;

	/**
	 * allocate data structures
	 */
	const count n = part.numberOfElements();
	const count z = g.upperNodeIdBound();
	const count k = part.numberOfSubsets();
	const count partZ = part.upperBound();

	if (k == 1) throw std::runtime_error("FM step with one partition does not really make sense.");

	if (nodeWeights.size() == 0) {
		nodeWeights.resize(n, 1);
	}
	assert(nodeWeights.size() == n);

	const auto subsetIds = part.getSubsetIds();
	vector<index> bestTargetPartition(z);
	vector<PrioQueue<edgeweight, index> > queues(part.upperBound(),n);
	vector<edgeweight> gains;
	vector<pair<index, index> > transfers;
	vector<index> transferedVertices;
	vector<double> imbalance;

	vector<double> fragmentSizes(part.upperBound());

	const edgeweight total = g.totalEdgeWeight();
	vector<bool> charged(z, false);
	vector<bool> chargedPart(part.upperBound(), false);

	double maxFragmentSize = 0;
	double sumNodeWeights = 0;
	double largestNodeWeight = 0;

	for (node v : g.nodes()) {
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

	/**
	 * The following extension of balance to non-uniform node weights is a deliberate overestimate. For arbitrary node weights, finding the optimal size amounts to solving Bin-Packing.
	 * Thus, we stick to the known definition for uniform weights and add a safety buffer otherwise.
	 */
	const double optSize = largestNodeWeight <= 1 ? ceil(double(n) / k) : double(sumNodeWeights) / k + largestNodeWeight;
	const double maxAllowablePartSize = optSize*(1+maxImbalance);


	/**
	 * if the number of nodes is not divisible by the number of partitions, a perfect balance is impossible.
	 * To avoid an endless loop, compute the theoretical minimum imbalance and adjust the parameter if necessary.
	 */
	assert(maxAllowablePartSize >= optSize);

	/**
	 * mark which partitions are charged already
	 */
	for (index c : chargedVertices) {
		charged[c] = true;
		chargedPart[part[c]] = true;
	}

	/**
	 * fill edge cut table
	 */
	vector<vector<double> > edgeCuts(n);
	g.parallelForNodes([&g, &part, &edgeCuts, partZ](index v){
		edgeCuts[v].resize(partZ, 0);
		for (index neighbor : g.neighbors(v)) {
			edgeCuts[v][part[neighbor]] += g.weight(v,neighbor);
		}
	});

	/**
	 * fill priority queues. Since we want to access the nodes with the maximum gain first but have only minQueues available, we use the negative gain as priority.
	 */
	g.forNodes([&queues, &part, &bestTargetPartition, &charged, &chargedPart, &chargedVertices, &edgeCuts, &fragmentSizes, &nodeWeights, partZ, k, total](index v){
		edgeweight maxCut = -total;
		index idAtMax = partZ;
		index ownFragment = part[v];

		for (index fragment = 0; fragment < partZ; fragment++) {
			if (fragment != ownFragment && edgeCuts[v][fragment] > maxCut && (!chargedPart[fragment] || !charged[v])) {//I could check the balance constraint even here. Should I?
				idAtMax = fragment;
				maxCut = edgeCuts[v][fragment];
			}
		}

		if (idAtMax < partZ) {
			/**
			 * usually, this should always be true. Only exception: exactly k charged nodes are given as input, and v is one of them.
			 */
			bestTargetPartition[v] = idAtMax;
			assert(ownFragment < queues.size());
			if (fragmentSizes[ownFragment] > nodeWeights[v]) {
				queues[ownFragment].insert(-(maxCut-edgeCuts[v][ownFragment]), v); //negative max gain
			}
		} else {
			assert(chargedVertices.size() == k);
		}
	});

	count queuedSum = 0;
	for (index i = 0; i < queues.size(); i++) {
		queuedSum += queues[i].size();
	}
	//assert(chargedVertices.size() == k || queuedSum == n);

	/**
	 * iterate over all vertices and move them, as long as unmoved vertices are left
	 */
	edgeweight gainsum = 0;
	bool allQueuesEmpty = false;

	vector<bool> moved(g.upperNodeIdBound(), false);

	while (!allQueuesEmpty) {
		allQueuesEmpty = true;

		//choose largest partition with non-empty queue.
		int largestMovablePart = -1;
		count largestSize = 0;

		for (index partID : subsetIds) {
			if (queues[partID].size() > 0 && fragmentSizes[partID] > largestSize) {
				largestMovablePart = partID;
				largestSize = fragmentSizes[partID];
			}
		}

		if (largestSize > 1 && largestMovablePart != -1) {
			//at least one queue is not empty
			allQueuesEmpty = false;
			index partID = largestMovablePart;

			assert(partID < queues.size());

			index topVertex;
			double topGain;
			std::tie(topGain, topVertex) = queues[partID].extractMin();
			topGain = -topGain;//invert, since the negative gain was used as priority.

			//now get target partition.
			index targetFragment = bestTargetPartition[topVertex];
			double storedGain = edgeCuts[topVertex][targetFragment] - edgeCuts[topVertex][partID];
			assert(abs(storedGain - topGain) < 0.0001);
			assert(fragmentSizes[partID] > nodeWeights[topVertex]);
			//double checkedGain = calculateGain(g, part, topVertex, targetFragment);
			//assert(abs(checkedGain - topGain) < 0.00001);

			//move node there
			TRACE("Moved node ", topVertex, " to partition ", targetFragment, " for gain of ", topGain);
			part.moveToSubset(targetFragment, topVertex);
			moved[topVertex] = true;

			//udpate size map
			fragmentSizes[partID] -= nodeWeights[topVertex];
			fragmentSizes[targetFragment] += nodeWeights[topVertex];

			//update charge state
			if (charged[topVertex]) {
				chargedPart[partID] = false;
				chargedPart[targetFragment] = true;
			}

			//update history
			gainsum += topGain;
			gains.push_back(gainsum);
			transfers.emplace_back(partID, targetFragment);
			transferedVertices.push_back(topVertex);

			imbalance.push_back(getWeightedImbalance(fragmentSizes, sumNodeWeights, largestNodeWeight, n, k));

			//update counter and possibly abort early
			if (topGain > 0) {
				movesWithoutImprovement = 0;
			} else {
				movesWithoutImprovement++;
				if (movesWithoutImprovement > KarypisKumarStoppingCriterion) {
					//break;//out of while-loop
				}
			}

			//update edge cuts of neighbours
			g.forNodes([&g, topVertex, partID, total, targetFragment, &queues, &part, &subsetIds, &moved, &bestTargetPartition, &charged, &chargedPart, &edgeCuts, &fragmentSizes, &nodeWeights, partZ](index w){
				if (!moved[w]) {
					//update gain
					edgeCuts[w][partID] -= g.weight(topVertex, w);
					edgeCuts[w][targetFragment] += g.weight(topVertex, w);

					edgeweight maxCut = -total;
					index idAtMax = partZ;
					index ownFragment = part[w];

					for (index fragment = 0; fragment < partZ; fragment++) {
						if (fragment != ownFragment && edgeCuts[w][fragment] > maxCut && (!chargedPart[fragment] || !charged[w])) {//I could check the balance constraint even here. Should I?
							idAtMax = fragment;
							maxCut = edgeCuts[w][fragment];
						}
					}

					bestTargetPartition[w] = idAtMax;

					//update prioqueue
					queues[part[w]].remove(w);
					if (idAtMax < partZ && fragmentSizes[ownFragment] > nodeWeights[w]) {
						queues[part[w]].insert(-(maxCut-edgeCuts[w][ownFragment]), w);
					}
				}
			});

		}
	}

	//g.forNodes([&moved](index v){assert(moved[v]);});

	count testedNodes = gains.size();
	if (testedNodes == 0) return 0;
	//assert(testedNodes == n);

	/**
	 * now find best partition among those tested
	 */
	int maxIndex = -1;
	edgeweight maxGain = 0;
	for (index i = 0; i < testedNodes; i++) {
		if (gains[i] > maxGain && imbalance[i] <= maxImbalance) {
			maxIndex = i;
			maxGain = gains[i];
		}
	}

	/**
	 * apply partition modifications in reverse until best is recovered
	 */
	for (int i = testedNodes-1; i > maxIndex; i--) {
		assert(part[transferedVertices[i]] == transfers[i].second);
		TRACE("Reversing move of ", transferedVertices[i], " from ", transfers[i].second, " back to ", transfers[i].first);
		part.moveToSubset(transfers[i].first, transferedVertices[i]);
	}
	return maxGain;
}


} /* namespace NetworKit */
