/*
 * Partitioner.cpp
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#include <tuple>
#include <queue>
#include <algorithm>

#include "Partitioner.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Random.h"

#include "../graph/GraphDistance.h"

using std::vector;
using std::pair;
using std::queue;
using Aux::PrioQueue;

namespace NetworKit {

Partitioner::Partitioner(const Graph& G, count numParts) : Algorithm(), G(G), numParts(numParts), result(0) {
	if (G.numberOfSelfLoops() > 0) throw std::runtime_error("Graph must not have self-loops.");
}

void Partitioner::run() {

	std::function<Partition(const Graph&)> partitionLambda = [&](const Graph& g) -> Partition {
	   count n = g.numberOfNodes();
	   DEBUG("Partitioning graph with ", n, " nodes into ", numParts, " parts.");

	   // coarsen recursively until graph is small enough
	   if (n <= 2 * numParts) {
		   Partition initial = recursiveBisection(g, numParts);

		   return initial;
	   }
	   else {
		   // recursive coarsening
		   LocalMaxMatcher matcher(g);
		   Matching matching = matcher.run();
		   assert(matching.isProper(g));
		   MatchingContracter coarsener(g, matching);
		   coarsener.run();
		   Graph coarseG = coarsener.getCoarseGraph();
		   std::vector<node> fineToCoarse = coarsener.getNodeMapping();

		   // recursive call
		   Partition coarsePart = partitionLambda(coarseG);

		   // interpolation
		   ClusteringProjector projector;
		   Partition finePart = projector.projectBack(coarseG, g, fineToCoarse, coarsePart);

		   // local refinement with Fiduccia-Matheyses
		   edgeweight gain;
		   do {
			    gain = fiducciaMatheysesStep(g, finePart);
		   } while (gain > 0);

		   return finePart;
	   }
   };
	result = partitionLambda(G);
	hasRun = true;
}

edgeweight Partitioner::fiducciaMatheysesStep(const Graph& g, Partition&  part) {
	/**
	 * allocate data structures
	 */
	const count k = part.numberOfSubsets();
	const count n = part.numberOfElements();
	const auto subsetIds = part.getSubsetIds();
	vector<index> bestTargetPartition(g.upperNodeIdBound());
	vector<PrioQueue<edgeweight, index> > queues(part.upperBound(),n);
	vector<edgeweight> gains;
	vector<pair<index, index> > transfers;
	vector<index> transferedVertices;
	edgeweight total = g.totalEdgeWeight();

	/**
	 * fill priority queues
	 */
	part.forEntries([total, &g, &part, &queues, &subsetIds, &bestTargetPartition](index node, index clusterID){
		assert(g.hasNode(node));
		edgeweight maxGain = -total;
		index IdAtMax = 0;
		for (index otherSubset : subsetIds) {
			edgeweight thisgain = calculateGain(g, part, node, otherSubset);
			if (thisgain > maxGain && otherSubset != clusterID) {
				IdAtMax = otherSubset;
				maxGain = thisgain;
			}
		}
		bestTargetPartition[node] = IdAtMax;
		assert(clusterID < queues.size());
		queues[clusterID].insert(maxGain, node);
	});

	count queuedSum = 0;
	for (index i = 0; i < queues.size(); i++) {
		queuedSum += queues[i].size();
	}
	assert(queuedSum == n);

	/**
	 * iterate over all vertices,
	 */
	edgeweight gainsum = 0;
	bool allQueuesEmpty = false;

	vector<bool> moved(g.upperNodeIdBound(), false);


	while (!allQueuesEmpty) {
		allQueuesEmpty = true;
		for (index partID : subsetIds) {
			assert(partID < queues.size());
			if (queues[partID].size() == 0) continue; //nothing to move here, queue is empty
			index topVertex;
			double topGain;
			std::tie(topGain, topVertex) = queues[partID].extractMin();

			index IdAtMax = bestTargetPartition[topVertex];
			//now get target partition.

			//move node there
			TRACE("Moved node ", topVertex, " to partition ", IdAtMax, " for gain of ", topGain);
			part.moveToSubset(IdAtMax, topVertex);
			moved[topVertex] = true;

			//update history
			gainsum += topGain;
			gains.push_back(gainsum);
			transfers.emplace_back(partID, IdAtMax);
			transferedVertices.push_back(topVertex);

			//update gains of neighbours
			g.forNeighborsOf(topVertex, [&g, topVertex, partID, total, &queues, &part, &subsetIds, &moved, &bestTargetPartition](index w){
				if (!moved[w]) {
					//update gain
					edgeweight newMaxGain = -total;
					index IdAtMax = 0;
					for (index otherSubset : subsetIds) {
						edgeweight thisgain = calculateGain(g, part, w, otherSubset);
						if (thisgain > newMaxGain && otherSubset != part[w]) {
							newMaxGain = thisgain;
							IdAtMax = otherSubset;
						}
					}
					bestTargetPartition[w] = IdAtMax;

					//update prioqueue
					queues[part[w]].remove(w);
					queues[part[w]].insert(newMaxGain, w);

				}
			});

			allQueuesEmpty = false;
		}
	}

	g.forNodes([&moved](index v){assert(moved[v]);});

	count testedNodes = gains.size();
	assert(testedNodes == n);



	/**
	 * now find best partition among those tested
	 */
	int maxIndex = 0;
	for (index i = 0; i < n; i++) {
		if (gains[i] > gains[maxIndex]) maxIndex = i;
	}

	/**
	 * apply partition modifications in reverse until best is recovered
	 */
	for (int i = n-1; i > maxIndex; i--) {
		assert(part[transferedVertices[i]] == transfers[i].second);
		TRACE("Reversing move of ", transferedVertices[i], " from ", transfers[i].second, " back to ", transfers[i].first);
		part.moveToSubset(transfers[i].first, transferedVertices[i]);
	}
	return gains[maxIndex];
}

edgeweight Partitioner::calculateGain(const Graph& g, const Partition& input, index u, index targetPart) {
	assert(input.numberOfElements() >= g.numberOfNodes());
	assert(g.hasNode(u));
	assert(input.contains(u));

	edgeweight extDegreeNow = 0;
	edgeweight extDegreeAfterMove = 0;

	const index oldPartition = input[u];

	g.forNeighborsOf(u, [&extDegreeNow, &extDegreeAfterMove, &u, &targetPart, &input, &g, &oldPartition](index v){
		if (input[v] != oldPartition) extDegreeNow += g.weight(u, v);
		if (input[v] != targetPart) extDegreeAfterMove  += g.weight(u, v);
	});

	return extDegreeNow - extDegreeAfterMove;
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

Partition Partitioner::recursiveBisection(const Graph& g, count k) {
	assert(k <= g.numberOfNodes());
	Partition trivialConstraint(g.upperNodeIdBound());
	trivialConstraint.allToOnePartition();
	recursiveBisection(g, k, trivialConstraint, trivialConstraint[0]);
	return trivialConstraint;
}

void Partitioner::recursiveBisection(const Graph& g, count k, Partition& mask, index maskID) {
	if (k == 1) return;
	auto beforeMap = mask.subsetSizeMap();
	count nodes = beforeMap.at(maskID);
	if (nodes < 2) return; //this risk returning less partitions then demanded.
	//assert(k <= nodes);

	index a, b;
	std::tie(a,b) = getMaximumDistancePair(g, mask, maskID);
	assert(a != b);
	const count firstWeight = k/2 + (k % 2 != 0 && Aux::Random::real() < 0.5) ? 1 : 0; //uneven numbers are randomly rounded to the first or the second half
	const count secondWeight = k - firstWeight;

	vector<index> points(2);
	vector<count> weights(2);
	points[0] = a;
	points[1] = b;
	weights[0] = firstWeight;
	weights[1] = secondWeight;

	/**
	 * TODO: we need a region growing implementation with custom weights
	 */
	mask = growRegions(g,  points, weights, mask);
	auto map = mask.subsetSizeMap();
	count firstRegionSize = map.at(a);
	count secondRegionSize = map.at(b);

	assert(firstRegionSize + secondRegionSize == nodes);

	//assert(firstRegionSize >= firstWeight);
	//assert(secondRegionSize >= secondWeight);

	recursiveBisection(g, firstWeight, mask, a);
	recursiveBisection(g, secondWeight, mask, b);
}

Partition Partitioner::growRegions(const Graph& g, const vector<index>& startingPoints) {
	Partition constraint(g.numberOfNodes());
	constraint.allToOnePartition();
	vector<count> weights(startingPoints.size(), 1);
	return growRegions(g, startingPoints, weights, constraint);
}

Partition Partitioner::growRegions(const Graph& g, const vector<index>& startingPoints, const vector<count>& weights, const Partition& constraint) {
	/**
	 * validate input
	 */
	const count n = g.numberOfNodes();
	const count z = g.upperNodeIdBound();
	assert(startingPoints.size() <= n);
	assert(startingPoints.size() == weights.size());
	assert(constraint.numberOfElements() == n);
	for (index point : startingPoints) assert(g.hasNode(point));

	/**
	 * allocate data structures
	 */
	vector<bool> visited(z, false);
	vector<queue<index>> bfsQueues(startingPoints.size());
	Partition result(constraint);
	result.setUpperBound(std::max(*std::max_element(startingPoints.begin(), startingPoints.end())+1, constraint.upperBound()));

	//TODO: shuffle partitions around randomly to avoid giving an advantage to the one with the lower index

	/**
	 * fill BFS queues with starting points
	 */
	for (index p = 0; p < startingPoints.size(); p++) {
		index point = startingPoints[p];
		bfsQueues[p].push(point);
		visited[point] = true;
		result.moveToSubset(point, point);
	}

	/**
	 * run BFS from sources and assign partitions
	 */
	bool allQueuesEmpty = false;
	while (!allQueuesEmpty) {
		allQueuesEmpty = true;
		for (index p = 0; p < bfsQueues.size(); p++) {
			if (bfsQueues[p].empty()) continue;
			allQueuesEmpty = false;

			index nextNode;

			do {
				nextNode = bfsQueues[p].front();
				bfsQueues[p].pop();
				assert(g.hasNode(nextNode));
			} while (!bfsQueues[p].empty() && visited[nextNode]  && startingPoints[p] != nextNode);

			if (visited[nextNode] && !(startingPoints[p] == nextNode)) continue;
			//if (!visited[nextNode]) {
				result.moveToSubset(startingPoints[p], nextNode);
				visited[nextNode] = true;
			//}

			for (index neighbor : g.neighbors(nextNode)) {
				if (visited[neighbor] || !constraint.inSameSubset(nextNode, neighbor)) continue;

				//if not visited and in same partition
				bfsQueues[p].push(neighbor);
			}
		}
	}

	//TODO: check that all nodes have been visited
	return result;
}

/**
 * optimization: give one node as parameter, there should be one lying around somewhere
 */
std::pair<index, index> Partitioner::getMaximumDistancePair(const Graph& g, const Partition& constraint, const index partition) {
	assert(partition < constraint.upperBound());
	assert(constraint.subsetSizeMap().at(partition) >= 2);
	assert(constraint.numberOfElements() == g.numberOfNodes());

	index z = g.upperNodeIdBound();
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	vector<edgeweight> distances(z, infDist);

	/**
	 * get node in correct partition
	 */
	index startingNode;
	g.forNodes([&constraint, partition, &startingNode](index u){
		if (constraint[u] == partition)
			startingNode = u;
			//this is wasteful, since it keeps iterating over the whole graph even if a node has been found
	});

	/**
	 * run breadth-first-searches constrained to the partition until distance no longer grows
	 */
	edgeweight maxDistance = 0;
	edgeweight lastDistance = -1;
	index a = startingNode;
	index b = startingNode;

	bool checkedConnectedness = false;

	while (maxDistance > lastDistance) {
		a = b;
		lastDistance = maxDistance;
		queue<index> bfsQueue;
		vector<bool> visited(z, false);
		distances.clear();
		distances.resize(z, infDist);

		bfsQueue.push(a);
		visited[a] = true;
		distances[a] = 0;

		while (!bfsQueue.empty()) {
			index currentNode = bfsQueue.front();
			bfsQueue.pop();

			for (index neighbor : g.neighbors(currentNode)) {
				if (visited[neighbor] || !constraint.inSameSubset(currentNode, neighbor)) continue;

				assert(distances[neighbor] == infDist);
				distances[neighbor] = distances[currentNode] + 1;

				//if not visited and in same partition
				visited[neighbor] = true;
				bfsQueue.push(neighbor);
			}
			b = currentNode;
			if (maxDistance < distances[currentNode]) {
				maxDistance = distances[currentNode];
			}
		}

		if (!checkedConnectedness) {
			for (index v : constraint.getMembers(constraint[startingNode])) {
				if (!visited[v]) {
					//partition is disconnected!
					return {a, v};
				}
			}
			checkedConnectedness = true;
		}

	}

	//assert(GraphDistance().unweightedDistance(g, a,b) == maxDistance); does not need to be true, there might be a shortcut through another partitition
	return {a, b};
}



} /* namespace NetworKit */
