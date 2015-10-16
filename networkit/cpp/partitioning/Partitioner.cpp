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
	   if (n <= 4 * numParts) {
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
	 * fill priority queues
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
	 * iterate over all vertices,
	 */
	edgeweight gainsum = 0;
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
			gainsum += maxGain;
			gains.push_back(gainsum);
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

	assert(gains.size() == n);

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
		part.moveToSubset(transferedVertices[i], transfers[i].first);
	}
	return gains[maxIndex];
}

edgeweight Partitioner::calculateGain(const Graph& g, const Partition& input, index u, index targetPart) {
	assert(input.numberOfElements() >= g.numberOfNodes());
	assert(g.hasNode(u));
	assert(input.contains(u));

	auto subsetIDs = input.getSubsetIds();
	assert(subsetIDs.count(targetPart) == 1);

	edgeweight extDegreeNow = 0;
	edgeweight extDegreeAfterMove = 0;

	const index oldPartition = input[u];

	g.forNeighborsOf(u, [&extDegreeNow, &extDegreeAfterMove, &u, &targetPart, &input, &g, &oldPartition](index v){
		if (input[v] != oldPartition) extDegreeNow += g.weight(u, v);
		if (input[v] != targetPart) extDegreeAfterMove  += g.weight(u, v);
	});

	return extDegreeAfterMove - extDegreeNow;
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

	index a, b;
	std::tie(a,b) = getMaximumDistancePair(g, mask, maskID);
	const count firstWeight = k/2;
	const count secondWeight = k - firstWeight;

	vector<index> points(2);
	vector<double> weights(2);
	points[0] = a;
	points[1] = b;
	weights[0] = firstWeight;
	weights[1] = secondWeight;

	/**
	 * here, we assume that the region growing overwrites the partitions of the input.
	 * Additionally, we need a region growing implementation with custom weights
	 */
	mask = growRegions(g,  points, mask);
	auto map = mask.subsetSizeMap();
	count firstRegionSize = map.at(a);
	count secondRegionSize = map.at(b);
	assert(firstRegionSize >= firstWeight);
	assert(secondRegionSize >= secondWeight);

	recursiveBisection(g, firstWeight, mask, a);
	recursiveBisection(g, secondWeight, mask, b);
}

Partition Partitioner::growRegions(const Graph& g, const vector<index>& startingPoints) {
	Partition constraint(g.numberOfNodes());
	constraint.allToOnePartition();
	return growRegions(g, startingPoints, constraint);
}

Partition Partitioner::growRegions(const Graph& g, const vector<index>& startingPoints, const Partition& constraint) {
	/**
	 * validate input
	 */
	const count n = g.numberOfNodes();
	const count z = g.upperNodeIdBound();
	assert(startingPoints.size() <= n);
	assert(constraint.numberOfElements() == n);
	for (index point : startingPoints) assert(g.hasNode(point));

	/**
	 * allocate data structures
	 */
	vector<bool> visited(z, false);
	queue<index> bfsQueue;
	Partition result(constraint);
	result.setUpperBound(std::max(*std::max_element(startingPoints.begin(), startingPoints.end())+1, constraint.upperBound()));

	/**
	 * fill BFS queue with starting points
	 */
	for (index point : startingPoints) {
		bfsQueue.push(point);
		visited[point] = true;
		result.moveToSubset(point, point);
	}

	/**
	 * run BFS from sources and assign partitions
	 */
	while (!bfsQueue.empty()) {
		index nextNode = bfsQueue.front();
		bfsQueue.pop();
		for (index neighbor : g.neighbors(nextNode)) {
			if (visited[neighbor] || !constraint.inSameSubset(nextNode, neighbor)) continue;

			//if not visited and in same partition
			visited[neighbor] = true;
			result.moveToSubset(result[nextNode], neighbor);
			bfsQueue.push(neighbor);
		}
	}

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
	}

	return {a, b};
}



} /* namespace NetworKit */
