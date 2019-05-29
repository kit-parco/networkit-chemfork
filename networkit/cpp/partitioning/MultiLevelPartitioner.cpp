/*
 * Partitioner.cpp
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#include <tuple>
#include <queue>
#include <algorithm>
#include <set>

#include "MultiLevelPartitioner.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Random.h"

#include "../components/ParallelConnectedComponents.h"

#include "../graph/GraphDistance.h"
#include "../community/ClusteringGenerator.h"

using std::vector;
using std::pair;
using std::queue;
using Aux::PrioQueue;

namespace NetworKit {

MultiLevelPartitioner::MultiLevelPartitioner(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const vector<index>& chargedVertices, count minGapSize, Partition previous) : GraphPartitioner(G, numParts, maxImbalance, chargedVertices, minGapSize), bisectRecursively(bisectRecursively) {
	if (bisectRecursively && chargedNodes.size() > 0) throw std::runtime_error("If using charged nodes, use region growing for the initial graph.");

	count n = G.numberOfNodes();
	if (previous.numberOfElements() == 0) {
		previousPartition = Partition(n);
		previousPartition.allToOnePartition();
		assert(previousPartition.numberOfSubsets() == 1);
	} else {
		if (previous.numberOfElements() != n) {
			throw std::runtime_error("Previous partition given, but of wrong size.");
		}
		if (previous.numberOfSubsets() > numParts) {
			throw std::runtime_error("Previous partition given, but with too many blocks.");
		}
		previousPartition = previous;
	}
}

void MultiLevelPartitioner::run() {
	count n = G.numberOfNodes();

	if (hasRun) {
		previousPartition = result;
	}

	std::vector<double> dummyWeights(n, 1);
	result = partitionRecursively(G, numParts, maxImbalance, bisectRecursively, chargedNodes, previousPartition, dummyWeights, minGapSize);
	if (result.numberOfSubsets() != numParts) {
		throw std::runtime_error("After recursive partitioning, got " + std::to_string(result.numberOfSubsets()) + " parts instead of " + std::to_string(numParts) + ".");
	}

	INFO("Cut after recursive Call: ", result.calculateCutWeight(G));

	/**
	 * make sure that the partition is balanced. Only necessary if the balance constraint was relaxed during the multi-level-partitioning
	 */
	enforceBalance(G, result, maxImbalance, chargedNodes, dummyWeights);
	fiducciaMattheysesStep(G, result, maxImbalance, chargedNodes, dummyWeights, minGapSize);
	INFO("Cut after rebalancing and final FM Step: ", result.calculateCutWeight(G));

	hasRun = true;
}

Partition MultiLevelPartitioner::partitionRecursively(const Graph& G, const count numParts, double maxImbalance, bool bisectRecursively,
		const std::vector<index>& chargedVertices, const Partition& previous, const std::vector<double> &nodeWeights, const count minGapSize) {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	assert(previous.numberOfElements() == n);
	assert(previous.numberOfSubsets() <= numParts);
	assert(nodeWeights.size() == n);

	DEBUG("Partitioning graph with ", n, " nodes, ", m, " edges and total edge weight ",  G.totalEdgeWeight() , " into ", numParts, " parts.");

	bool coarsestLevelReached = n <= 2 * numParts;

	//get cut edges
	std::vector<std::pair<node, node> > forbiddenEdges;

	if (!coarsestLevelReached) {
		G.forEdges([&](node u, node v, edgeweight w) {
			if (previous[u] != previous[v]) {
				forbiddenEdges.push_back(std::make_pair(u,v));
			}
		});
		if (forbiddenEdges.size() == G.numberOfEdges()) {
			//no more edge to contract, use initial partitioning step
			coarsestLevelReached = true;
		}
	}

	if (!coarsestLevelReached) {
		LocalMaxMatcher matcher(G, chargedVertices, forbiddenEdges, nodeWeights, true);
		matcher.run();
		if (matcher.getMatching().size(G) == 0) {
			coarsestLevelReached = true;
		}
	}

	if (coarsestLevelReached) {
	   Partition initial;
	   if (bisectRecursively) {
		   initial = recursiveBisection(G, numParts);
		   assert(initial.numberOfSubsets() == numParts);
	   } else {
		   vector<index> startingPoints(chargedVertices);

		   /**
		    * fill up starting points with other points
		    */
		   for (index i = chargedVertices.size(); i < numParts; i++) {
			   index farthestNode = getFarthestNode(G, startingPoints);
			   startingPoints.push_back(farthestNode);
			}
		   initial = growRegions(G, startingPoints);//TODO: adapt region growing to accept node weights
		   assert(initial.numberOfSubsets() == numParts);
	   }

	   ClusteringGenerator gen;
	   Partition naiveInitial = gen.makeContinuousBalancedClustering(G, numParts);

	   if (chargesValid(naiveInitial, chargedVertices) && naiveInitial.calculateCutWeight(G) < initial.calculateCutWeight(G)
			   && getWeightedImbalance(initial, nodeWeights, numParts) <= maxImbalance && naiveInitial.numberOfSubsets() == numParts) {
		   initial = naiveInitial;
		   DEBUG("Replaced initial partition with naive solution, since it was better.");
	   }

	   bool previousValid = previous.numberOfElements() == n && getWeightedImbalance(previous, nodeWeights, numParts) <= maxImbalance && chargesValid(previous, chargedVertices)
	   	   && previous.numberOfSubsets() == numParts;
	   if (previousValid && previous.calculateCutWeight(G) < initial.calculateCutWeight(G)) {
		   initial = previous;
	   }

	   count initialK = initial.numberOfSubsets();
	   DEBUG("Initial solution has ", initialK, " blocks, a cut of ", initial.calculateCutWeight(G), " and an imbalance of ", initial.getImbalance(numParts));
	   assert(initialK == numParts);
	   return initial;
	}
	else {
		// recursive coarsening
	   LocalMaxMatcher matcher(G, chargedVertices, forbiddenEdges, nodeWeights, true);
	   matcher.run();
	   Matching matching = matcher.getMatching();
	   assert(matching.isProper(G));
	   MatchingContracter coarsener(G, matching);
	   coarsener.run();
	   Graph coarseG = coarsener.getCoarseGraph();
	   assert(coarseG.numberOfNodes() < G.numberOfNodes());
	   if (coarseG.numberOfNodes() >= G.numberOfNodes()) {
			throw std::runtime_error("Graph not smaller after coarsening.");
	   }

	   std::vector<node> fineToCoarse = coarsener.getFineToCoarseNodeMapping();

	   Partition coarsePrevious(coarseG.numberOfNodes());
	   coarsePrevious.allToSingletons();

	   //map node weights and previous partition
	   vector<double> coarseWeights(coarseG.numberOfNodes(), 0);

	   for (node v : G.nodes()) {
		   node coarseNode = fineToCoarse[v];
		   index coarsePart = previous[v];
		   coarsePrevious.moveToSubset(coarsePart, coarseNode);

		   coarseWeights[coarseNode] += nodeWeights[v];
	   }
	   coarsePrevious.compact();
	   assert(coarsePrevious.numberOfSubsets() == previous.numberOfSubsets());

	   // map charged vertices to coarser nodes
	   vector<index> coarseCharged;
	   for (node v : chargedVertices) {
		   coarseCharged.push_back(fineToCoarse[v]);
	   }

	   // recursive call
	   Partition coarsePart = partitionRecursively(coarseG, numParts, maxImbalance, bisectRecursively, coarseCharged, coarsePrevious, coarseWeights, minGapSize);
		if (coarsePart.numberOfSubsets() != numParts) {
			throw std::runtime_error("Coarse partition has " + std::to_string(coarsePart.numberOfSubsets()) + " parts instead of " + std::to_string(numParts) + ".");
		}

	   // interpolation
	   ClusteringProjector projector;
	   Partition finePart = projector.projectBack(coarseG, G, fineToCoarse, coarsePart);

	   edgeweight preBalancingCut = finePart.calculateCutWeight(G);
	   double preBalancingImbalance = finePart.getImbalance(numParts);
	   count preBalancingK = finePart.numberOfSubsets();

	   enforceBalance(G, finePart, maxImbalance*2, chargedVertices, nodeWeights);

	   edgeweight preRefinementCut = finePart.calculateCutWeight(G);
	   double preRefinementImbalance = finePart.getImbalance(numParts);
	   count preRefinementK = finePart.numberOfSubsets();
	   assert(preBalancingK == preRefinementK);

		DEBUG("Rebalancing, n: ", n, ", cut: ", preBalancingCut, "->", preRefinementCut, ", imbalance:", preBalancingImbalance, "->", preRefinementImbalance);

	   // local refinement with Fiduccia-Matheyses
	   edgeweight gain;
	   do {
			gain = fiducciaMattheysesStep(G, finePart, maxImbalance, chargedVertices, nodeWeights, minGapSize);
			assert(gain == gain);
			TRACE("Found gain ", gain, " in FM-step with ", G.numberOfNodes(), " nodes and ", finePart.numberOfSubsets(), " partitions.");
	   } while (gain > 1e-10);
	   assert(gain == 0);
	   edgeweight postRefinementCut = finePart.calculateCutWeight(G);
	   assert(postRefinementCut <= preRefinementCut);

	   INFO("Refinement, n: ", G.numberOfNodes(), " k: ", preRefinementK, "->", finePart.numberOfSubsets(), ", cut: ", preRefinementCut, "->", postRefinementCut, ", imbalance:", preRefinementImbalance, "->", finePart.getImbalance(numParts));

	   assert(finePart.numberOfSubsets() == numParts);
	   return finePart;
	}
}
edgeweight MultiLevelPartitioner::calculateGain(const Graph& g, const Partition& input, index u, index targetPart) {
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

std::string MultiLevelPartitioner::toString() const {
	return "TODO";
}

Partition MultiLevelPartitioner::recursiveBisection(const Graph& g, count k) {
	assert(k <= g.numberOfNodes());
	Partition trivialConstraint(g.upperNodeIdBound());
	trivialConstraint.allToOnePartition();
	recursiveBisection(g, k, trivialConstraint, trivialConstraint[0]);
	return trivialConstraint;
}

void MultiLevelPartitioner::recursiveBisection(const Graph& g, count k, Partition& mask, index maskID) {
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

Partition MultiLevelPartitioner::growRegions(const Graph& g, const vector<index>& startingPoints) {
	Partition constraint(g.numberOfNodes());
	constraint.allToOnePartition();
	vector<count> weights(startingPoints.size(), 1);
	return growRegions(g, startingPoints, weights, constraint);
}

Partition MultiLevelPartitioner::growRegions(const Graph& g, const vector<index>& startingPoints, const vector<count>& weights, const Partition& constraint) {
	/**
	 * validate input
	 */
	const count n = g.numberOfNodes();
	const count z = g.upperNodeIdBound();
	const count k = startingPoints.size();
	assert(startingPoints.size() <= n);
	assert(startingPoints.size() == weights.size());
	assert(constraint.numberOfElements() == n);
	for (index point : startingPoints) assert(g.hasNode(point));

	/**
	 * make sure starting points are unique
	 */
	for (index i = 0; i < k; i++) {
		for (index j = 0; j < i; j++) {
			assert(startingPoints[i] != startingPoints[j]);
		}
	}

	/**
	 * allocate data structures
	 */
	vector<bool> visited(z, false);
	vector<queue<index>> bfsQueues(k);
	vector<count> partitionSizes(k, 0);

	//partitions already present in the input are copied
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
			if (bfsQueues[p].empty()) continue; //here one could also check whether the partition has already reached its alloted size
			allQueuesEmpty = false;

			index currentNode;

			do {
				currentNode = bfsQueues[p].front();
				bfsQueues[p].pop();
				assert(g.hasNode(currentNode));
			} while (!bfsQueues[p].empty() && visited[currentNode]  && startingPoints[p] != currentNode);

			if (visited[currentNode] && !(startingPoints[p] == currentNode)) continue;

			result.moveToSubset(startingPoints[p], currentNode);
			visited[currentNode] = true;
			partitionSizes[p]++;

			for (index neighbor : g.neighbors(currentNode)) {
				if (visited[neighbor] || !constraint.inSameSubset(currentNode, neighbor)) continue;

				//if not visited and in same partition
				bfsQueues[p].push(neighbor);
			}
		}
	}
	result.compact();
	assert(result.numberOfSubsets() == startingPoints.size());
	//TODO: make sure that all nodes in this subpartition have been visited
	//g.forNodes([&visited](index v){assert(visited[v]);});
	return result;
}

/**
 * optimization: give one node as parameter, there should be one lying around somewhere
 */
std::pair<index, index> MultiLevelPartitioner::getMaximumDistancePair(const Graph& g, const Partition& constraint, const index partition) {
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
					/**
					 * partition is disconnected!
					 * The distance between a and all nodes in the other component is infinity, we can just return any of them.
					 */
					return {a, v};
				}
			}
			checkedConnectedness = true;
		}

	}

	//assert(GraphDistance().unweightedDistance(g, a,b) == maxDistance); does not need to be true, there might be a shortcut through another partition
	return {a, b};
}

index MultiLevelPartitioner::getFarthestNode(const Graph& G, std::vector<index> seedNodes) {
	/**
	 * Yet another BFS. This currently has problems with unconnected graphs.
	 */
	const count z = G.upperNodeIdBound();

	if (seedNodes.size() == 0) return G.randomNode();

	vector<bool> visited(z, false);
	queue<index> bfsQueue;

	for (index seed : seedNodes) {
		bfsQueue.push(seed);
		assert(G.hasNode(seed));
		visited[seed] = true;
	}

	index nextNode = G.randomNode(); //will be overwritten anyway
	while (bfsQueue.size() > 0) {
		nextNode = bfsQueue.front();
		bfsQueue.pop();
		visited[nextNode] = true;
		G.forNeighborsOf(nextNode, [&visited, &bfsQueue](index v){if (!visited[v]) bfsQueue.push(v);});
	}

	//if nodes are unvisited, the graph is unconnected and the unvisited nodes are in fact the farthest
	for (node v : G.nodes()) {
		if (!visited[v]) nextNode = v;
		break;
	}

	return nextNode;
}

void MultiLevelPartitioner::repairSingleNodes(const Graph& G, Partition& intermediate) {
	const count n = G.numberOfNodes();
	assert(intermediate.numberOfElements() == n);

	//TODO: this has problems with deleted nodes
	for (index i = 1; i < n-1; i++) {
		index lastPart = intermediate[i-1];
		index currentPart = intermediate[i];
		index nextPart = intermediate[i+1];
		if (lastPart == nextPart && lastPart != currentPart) {
			/**
			 * we have a single node wedged into a wrong partition. We now have three choices:
			 * 1. Move surrounded node to enclosing partition
			 * 2. Move left node to some other partition
			 * 3. Move right node to some other partition
			 */

			intermediate.moveToSubset(lastPart, i);
			DEBUG("Moved node ", i, " to subset ", lastPart, " as it was surrounded by it.");
		}
	}
}


} /* namespace NetworKit */
