/*
 * Partitioner.cpp
 *
 *  Created on: 25.09.2015
 *      Author: moritzl
 */

#include <tuple>
#include <queue>
#include <algorithm>

#include "MultiLevelPartitioner.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Random.h"

#include "../graph/GraphDistance.h"
#include "../community/ClusteringGenerator.h"

using std::vector;
using std::pair;
using std::queue;
using Aux::PrioQueue;

namespace NetworKit {

MultiLevelPartitioner::MultiLevelPartitioner(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const vector<index>& chargedVertices, bool avoidSingleNodes, Partition previous) : Algorithm(), G(G), numParts(numParts), maxImbalance(maxImbalance), bisectRecursively(bisectRecursively), chargedNodes(chargedVertices), noSingles(avoidSingleNodes), result(0) {
	if (G.numberOfSelfLoops() > 0) throw std::runtime_error("Graph must not have self-loops.");
	if (chargedNodes.size() > numParts) throw std::runtime_error("Cannot have more charged nodes than partitions.");
	if (bisectRecursively && chargedNodes.size() > 0) throw std::runtime_error("If using charged nodes, use region growing for the initial graph.");
	for (index i : chargedVertices) {
		if (!G.hasNode(i)) throw std::runtime_error("At least one of the charged nodes is missing from the graph.");
	}
	count n = G.numberOfNodes();
	if (previous.numberOfElements() == 0) {
		previousPartition = Partition(n);
		previousPartition.allToOnePartition();
		assert(previousPartition.numberOfSubsets() == 1);
	} else {
		if (previous.numberOfElements() != n) {
			throw std::runtime_error("Previous partition given, but of wrong size.");
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
	result = partitionRecursively(G, numParts, maxImbalance, bisectRecursively, chargedNodes, previousPartition, dummyWeights);

	/**
	 * make sure that the partition is balanced. Only necessary if the balance constraint was relaxed during the multi-level-partitioning
	 */
	enforceBalance(G, result, maxImbalance, chargedNodes, dummyWeights);
	fiducciaMattheysesStep(G, result, maxImbalance, chargedNodes, dummyWeights);

	if (noSingles) repairSingleNodes(G, result);

	hasRun = true;
}

Partition MultiLevelPartitioner::partitionRecursively(const Graph& G, count numParts, double maxImbalance, bool bisectRecursively, const std::vector<index>& chargedVertices, const Partition& previous, const std::vector<double> &nodeWeights) {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	assert(previous.numberOfElements() == n);
	assert(nodeWeights.size() == n);

	DEBUG("Partitioning graph with ", n, " nodes, ", m, " edges and total edge weight ",  G.totalEdgeWeight() , " into ", numParts, " parts.");

	// coarsen recursively until graph is small enough
	if (n <= 2 * numParts) {
	   Partition initial;
	   if (bisectRecursively) {
		   initial = recursiveBisection(G, numParts);
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
	   }

	   ClusteringGenerator gen;
	   Partition naiveInitial = gen.makeContinuousBalancedClustering(G, numParts);

	   if (chargesValid(naiveInitial, chargedVertices) && naiveInitial.calculateCutWeight(G) < initial.calculateCutWeight(G) && getWeightedImbalance(initial, nodeWeights, numParts) <= maxImbalance) {
		   initial = naiveInitial;
		   DEBUG("Replaced initial partition with naive solution, since it was better.");
	   }

	   bool previousValid = previous.numberOfElements() == n && getWeightedImbalance(previous, nodeWeights, numParts) <= maxImbalance && chargesValid(previous, chargedVertices);
	   if (previousValid && previous.calculateCutWeight(G) < initial.calculateCutWeight(G)) {
		   initial = previous;
	   }

	   count initialK = initial.numberOfSubsets();
	   DEBUG("Initial solution has ", initialK, " partitions, a cut of ", initial.calculateCutWeight(G), " and an imbalance of ", initial.getImbalance(numParts));
	   return initial;
	}
	else {

		//get cut edges
		std::vector<std::pair<node, node> > forbiddenEdges;
		G.forEdges([&](node u, node v, edgeweight w) {
			if (previous[u] != previous[v]) {
				forbiddenEdges.push_back(std::make_pair(u,v));
			}
		});
		assert(forbiddenEdges.size() < G.numberOfEdges());

		// recursive coarsening
	   LocalMaxMatcher matcher(G, chargedVertices, forbiddenEdges, nodeWeights, true);
	   matcher.run();
	   Matching matching = matcher.getMatching();
	   assert(matching.isProper(G));
	   MatchingContracter coarsener(G, matching);
	   coarsener.run();
	   Graph coarseG = coarsener.getCoarseGraph();
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
	   Partition coarsePart = partitionRecursively(coarseG, numParts, maxImbalance, bisectRecursively, coarseCharged, coarsePrevious, coarseWeights);

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
			gain = fiducciaMattheysesStep(G, finePart, maxImbalance, chargedVertices, nodeWeights);
			assert(gain == gain);
			TRACE("Found gain ", gain, " in FM-step with ", G.numberOfNodes(), " nodes and ", finePart.numberOfSubsets(), " partitions.");
	   } while (gain > 0);
	   assert(gain == 0);
	   edgeweight postRefinementCut = finePart.calculateCutWeight(G);
	   assert(postRefinementCut <= preRefinementCut);

	   DEBUG("Refinement, n: ", G.numberOfNodes(), " k: ", preRefinementK, "->", finePart.numberOfSubsets(), ", cut: ", preRefinementCut, "->", postRefinementCut, ", imbalance:", preRefinementImbalance, "->", finePart.getImbalance(numParts));

	   return finePart;
	}
}
//TODO: the node weights are currently copied, this could be improved. I could use a reference, but then I cannot make it const, since I resize the vector if it is empty. I allow empty nodeWeights for overloading.
edgeweight MultiLevelPartitioner::fiducciaMattheysesStep(const Graph& g, Partition&  part, double maxImbalance, const std::vector<index> chargedVertices, std::vector<double> nodeWeights) {
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
	g.forNodes([&queues, &part, &bestTargetPartition, &charged, &chargedPart, &chargedVertices, &edgeCuts, partZ, k, total](index v){
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
			queues[ownFragment].insert(-(maxCut-edgeCuts[v][ownFragment]), v); //negative max gain
		} else {
			assert(chargedVertices.size() == k);
		}
	});

	count queuedSum = 0;
	for (index i = 0; i < queues.size(); i++) {
		queuedSum += queues[i].size();
	}
	assert(chargedVertices.size() == k || queuedSum == n);

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
			g.forNeighborsOf(topVertex, [&g, topVertex, partID, total, targetFragment, &queues, &part, &subsetIds, &moved, &bestTargetPartition, &charged, &chargedPart, &edgeCuts, partZ](index w){
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
					if (idAtMax < partZ) {
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

Partition MultiLevelPartitioner::getPartition() {
	if(!hasRun) {
		throw std::runtime_error("Call run()-function first.");
	}
	return result;
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

	//assert(GraphDistance().unweightedDistance(g, a,b) == maxDistance); does not need to be true, there might be a shortcut through another partitition
	return {a, b};
}

index MultiLevelPartitioner::getFarthestNode(const Graph& G, std::vector<index> seedNodes) {
	/**
	 * Yet another BFS.
	 */
	const count z = G.upperNodeIdBound();

	vector<bool> visited(z, false);
	queue<index> bfsQueue;

	for (index seed : seedNodes) {
		bfsQueue.push(seed);
		assert(G.hasNode(seed));
	}

	index nextNode = G.randomNode(); //will be overwritten anyway, except if seedNodes was empty
	while (bfsQueue.size() > 0) {
		nextNode = bfsQueue.front();
		bfsQueue.pop();
		visited[nextNode] = true;
		G.forNeighborsOf(nextNode, [&visited, &bfsQueue](index v){if (!visited[v]) bfsQueue.push(v);});
	}
	return nextNode;
}

void MultiLevelPartitioner::enforceBalance(const Graph& g, Partition& part, double maxImbalance, const vector<index>& chargedVertices, const std::vector<double>& nodeWeights) {
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
	const double maxAllowablePartSize = optSize*(1+maxImbalance);

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
