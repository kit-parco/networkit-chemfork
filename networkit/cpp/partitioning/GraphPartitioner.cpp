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

GraphPartitioner::GraphPartitioner(const Graph& G, count numParts, double maxImbalance, std::vector<index> chargedVertices, count minGapSize) : G(G), numParts(numParts), maxImbalance(maxImbalance), chargedNodes(chargedVertices), minGapSize(minGapSize), result(0) {
	if (G.numberOfSelfLoops() > 0) throw std::runtime_error("Graph must not have self-loops.");
	if (chargedNodes.size() > numParts) throw std::runtime_error("Cannot have more charged nodes than partitions.");
	std::set<index> chargedVerticesSet(chargedVertices.begin(), chargedVertices.end());
	if (chargedVerticesSet.size() < chargedVertices.size()) {
		throw std::runtime_error("Charged vertices are not unique.");
	}

	for (index i : chargedVertices) {
		if (!G.hasNode(i)) throw std::runtime_error("At least one of the charged nodes is missing from the graph.");
	}

	ParallelConnectedComponents comp(G);
	comp.run();
	if (comp.numberOfComponents() > 1) {
		throw std::runtime_error("Graph is unconnected..");
	}

	if (minGapSize > 3) {
		throw std::logic_error("Forbidden gap sizes of > 3 are not yet implemented.");
	}
}


Partition GraphPartitioner::getPartition() const {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return result;
}

//TODO: the node weights are currently copied, this could be improved. I could use a reference, but then I cannot make it const, since I resize the vector if it is empty. I allow empty nodeWeights for overloading.
edgeweight GraphPartitioner::fiducciaMattheysesStep(const Graph& g, Partition&  part, double maxImbalance, const std::vector<index> &chargedVertices, std::vector<double> nodeWeights, count minGapSize) {
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
	if (*std::min_element(nodeWeights.begin(), nodeWeights.end()) < 1) {
		throw std::runtime_error("Node weights must be at least 1.");
	}


	const auto subsetIds = part.getSubsetIds();
	vector<index> bestTargetBlock(z);
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
	assert(maxAllowablePartSize >= optSize);

	/**
	 * mark which partitions are charged already
	 */
	for (index c : chargedVertices) {
		charged[c] = true;
		chargedPart[part[c]] = true;
	}

	auto moveLegal = [&](index node, index targetBlock) {
		if (part[node] == targetBlock) return false;
		if (chargedPart[targetBlock] && charged[node]) return false;
		if (fragmentSizes[targetBlock] + nodeWeights[node] > maxAllowablePartSize) return false;

		if (node -2 >= 0 && part[node-2] == targetBlock && part[node-1] != targetBlock && nodeWeights[node-1] < minGapSize) return false;
		if (node > 0 && node < n - 1 && part[node-1] == part[node+1] && targetBlock != part[node-1] && nodeWeights[node] < minGapSize) return false;
		if (node + 2 < n && part[node+2] == targetBlock && part[node+1] != targetBlock && nodeWeights[node+1] < minGapSize) return false;

		if (minGapSize == 3) {
			bool noGaps = true;
			if (node - 3 >= 0 && part[node-3] == targetBlock && part[node-2] != targetBlock && nodeWeights[node-2] + nodeWeights[node-1] < minGapSize) noGaps = false;
			if (node - 2 >= 0 && node + 1 < n && part[node-2] == part[node+1] && targetBlock != part[node+1] && nodeWeights[node] + nodeWeights[node-1] < minGapSize) noGaps = false;
			if (node - 1 >= 0 && node + 2 < n && part[node-1] == part[node+2] && targetBlock != part[node-1] && nodeWeights[node] + nodeWeights[node+1] < minGapSize) noGaps = false;
			if (node + 3 < n  && part[node+3] == targetBlock && part[node+2] != targetBlock && nodeWeights[node+2] + nodeWeights[node+1] < minGapSize) noGaps = false;
			if (!noGaps) {
				return false;
			}
		}

		return true;
	};

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
	g.forNodes([&](index v){
		edgeweight maxCut = -total;
		index idAtMax = partZ;
		index ownFragment = part[v];

		for (index fragment = 0; fragment < partZ; fragment++) {
			if (moveLegal(v, fragment) && edgeCuts[v][fragment] > maxCut) {
				idAtMax = fragment;
				maxCut = edgeCuts[v][fragment];
			}
		}

		if (idAtMax < partZ) {
			/**
			 * usually, this should always be true. Only exception: exactly k charged nodes are given as input, and v is one of them.
			 */
			bestTargetBlock[v] = idAtMax;
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
			index targetFragment = bestTargetBlock[topVertex];
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
			g.forNodes([&](index w){
				if (!moved[w]) {
					//update gain
					edgeCuts[w][partID] -= g.weight(topVertex, w);
					edgeCuts[w][targetFragment] += g.weight(topVertex, w);

					edgeweight maxCut = -total;
					index idAtMax = partZ;
					index ownFragment = part[w];

					for (index fragment = 0; fragment < partZ; fragment++) {
						if (moveLegal(w, fragment) && edgeCuts[w][fragment] > maxCut) {
							idAtMax = fragment;
							maxCut = edgeCuts[w][fragment];
						}
					}

					bestTargetBlock[w] = idAtMax;

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

void GraphPartitioner::enforceBalance(const Graph& g, Partition& part, double maxImbalance, const vector<index>& chargedVertices, const std::vector<double>& nodeWeights) {
	const count n = g.numberOfNodes();
	const count z = g.upperNodeIdBound();
	const count k = part.numberOfSubsets();
	if (part.upperBound() > k) {
		part.compact();
	}
	const count partZ = part.upperBound();
	assert(part.numberOfElements() == n);

	/**
	 * allocate data structures similar to FM
	 */
	std::vector<PrioQueue<edgeweight, index> > queues(part.upperBound(),n);
	const auto subsetIds = part.getSubsetIds();
	vector<index> bestTargetPartition(z, 0);
	vector<double> fragmentSizes(part.upperBound());

	const edgeweight total = g.totalEdgeWeight();
	vector<bool> charged(z, false);
	vector<bool> chargedPart(part.upperBound(), false);
	vector<bool> moved(g.upperNodeIdBound(), false);

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

		if (fragmentSizes[part[v]] > maxFragmentSize) {
			maxFragmentSize = fragmentSizes[part[v]];
		}
	}

	/**
	 * The following extension of balance to non-uniform node weights is a deliberate overestimate. For arbitrary node weights, finding the optimal size amounts to solving Bin-Packing.
	 * Thus, we stick to the known definition for uniform weights and add a safety buffer otherwise.
	 */
	const double optSize = largestNodeWeight <= 1 ? ceil(double(n) / k) : double(sumNodeWeights) / k + largestNodeWeight;
	//const double maxAllowablePartSize = optSize*(1+maxImbalance);

	/**
	 * mark which partitions are charged already
	 */
	for (index c : chargedVertices) {
		charged[c] = true;
		chargedPart[part[c]] = true;
	}

	/**
	 *
	 */


	/**
	 * save values before rebalancing to compare
	 */
	double currentImbalance = getWeightedImbalance(fragmentSizes, sumNodeWeights, largestNodeWeight, n, k);

	/**
	 * repeat until partitions are balanced
	 */
	while(currentImbalance > maxImbalance) {
		/**
		 * reset priority queues
		 */
		for (auto queue : queues) {
			queue.clear();
		}

		/**
		 * reset moved status
		 */
		moved.clear();
		moved.resize(z, false);

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
		g.forNodes([&queues, &part, &bestTargetPartition, &charged, &chargedPart, &edgeCuts, &fragmentSizes, &nodeWeights, partZ, total](index v){
			edgeweight maxCut = -total;
			index idAtMax = partZ;
			index ownFragment = part[v];
			for (index fragment = 0; fragment < partZ; fragment++) {
				if (fragment != ownFragment && edgeCuts[v][fragment] > maxCut && (!chargedPart[fragment] || !charged[v]) && fragmentSizes[fragment] + nodeWeights[v] < fragmentSizes[ownFragment]) {//I could check the balance constraint even here. Should I?
					idAtMax = fragment;
					maxCut = edgeCuts[v][fragment];
				}
			}

			assert(ownFragment < queues.size());
			bestTargetPartition[v] = idAtMax;
			if (idAtMax < partZ) {
				//if movable, put into queue
				double key = -(maxCut-edgeCuts[v][ownFragment]);
				queues[ownFragment].insert(key, v); //negative max gain
			}

		});

		/**
		 * main movement loop
		 */
		bool allQueuesEmpty = false;
		while(!allQueuesEmpty) {
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

			if (largestMovablePart != -1) {
				index p = largestMovablePart;
				allQueuesEmpty = false;
				while(queues[p].size() > 0) {

					/**
					 * the loop content is copied almost verbatim from the FM-algorithm
					 */

					index topVertex;
					double topGain;
					std::tie(topGain, topVertex) = queues[p].extractMin();
					assert(!moved[topVertex]);
					topGain = -topGain;//invert, since the negative gain was used as priority.

					if (topGain <= -total) continue;//could not be moved anywhere

					//now get target partition.
					index targetFragment = bestTargetPartition[topVertex];
					double storedGain = edgeCuts[topVertex][targetFragment] - edgeCuts[topVertex][p];
					assert(abs(storedGain - topGain) < 0.0001);
					//double checkedGain = calculateGain(g, part, topVertex, targetFragment);
					//assert(abs(checkedGain - topGain) < 0.0001);

					count ownSize = fragmentSizes[p];
					count targetSize = fragmentSizes[targetFragment];

					if (targetSize >= ownSize) {
						//out of date, update gain

						edgeweight maxCut = -total;
						index idAtMax = partZ;
						for (index fragment = 0; fragment < partZ; fragment++) {
							if (fragment != p && edgeCuts[topVertex][fragment] > maxCut && (!chargedPart[fragment] || !charged[topVertex])  && fragmentSizes[fragment] + nodeWeights[topVertex] < fragmentSizes[p]) {
								idAtMax = fragment;
								maxCut = edgeCuts[topVertex][fragment];
							}
						}
						if (idAtMax < partZ) {
							/**
							 * the node can be moved, but to another fragment
							 * reinsert into queue
							 */
							double newGain = maxCut - edgeCuts[topVertex][p];
							bestTargetPartition[topVertex] = idAtMax;
							queues[p].insert(-newGain, topVertex);
						}

						continue;
					}

					part.moveToSubset(targetFragment, topVertex);
					moved[topVertex] = true;

					DEBUG("Node ", topVertex, " of weight ", nodeWeights[topVertex], " moved, Imbalance now ", getWeightedImbalance(fragmentSizes, sumNodeWeights, largestNodeWeight,  n,  k));

					//udpate size map
					fragmentSizes[p] -= nodeWeights[topVertex];
					fragmentSizes[targetFragment] += nodeWeights[topVertex];

					//update charge state
					if (charged[topVertex]) {
						chargedPart[p] = false;
						chargedPart[targetFragment] = true;
					}

					//update edge cuts of neighbours
					g.forNeighborsOf(topVertex, [&g, topVertex, total, p, targetFragment, &queues, &part, &subsetIds, &moved, &bestTargetPartition, &charged, &chargedPart, &edgeCuts, partZ, &nodeWeights, &fragmentSizes](index w){
						if (!moved[w]) {
							//update gain
							edgeCuts[w][p] -= g.weight(topVertex, w);
							edgeCuts[w][targetFragment] += g.weight(topVertex, w);

							edgeweight maxCut = -total;
							index idAtMax = partZ;
							index ownFragment = part[w];

							for (index fragment = 0; fragment < partZ; fragment++) {
								if (fragment != ownFragment && edgeCuts[w][fragment] > maxCut && (!chargedPart[fragment] || !charged[w])  && fragmentSizes[fragment] + nodeWeights[w] < fragmentSizes[part[w]]) {
									idAtMax = fragment;
									maxCut = edgeCuts[w][fragment];
								}
							}

							bestTargetPartition[w] = idAtMax;

							//update prioqueue
							queues[part[w]].remove(w);
							if (idAtMax < partZ) {
								queues[part[w]].insert(-(maxCut-edgeCuts[w][ownFragment]), w);
							}

						}
					});
				}
			}
		}
		double oldImbalance = currentImbalance;
		currentImbalance = getWeightedImbalance(fragmentSizes, largestNodeWeight, sumNodeWeights, n,  k);
		DEBUG("Imbalance reduced ", oldImbalance, "->", currentImbalance);
		assert(currentImbalance < oldImbalance);
	}
}



} /* namespace NetworKit */
