/*
 * PartitionerGTest.cpp
 *
 *  Created on: 15.10.2015
 *      Author: moritzl
 */

#include "PartitionerGTest.h"

#include "../Partitioner.h"
#include "../../graph/Graph.h"
#include "../../Globals.h"
#include "../../generators/BarabasiAlbertGenerator.h"
#include "../../generators/HyperbolicGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../../community/ClusteringGenerator.h"

namespace NetworKit {

TEST_F(PartitionerGTest, testGainUnweighted) {
	Graph G(6);
	G.addEdge(0,1);
	G.addEdge(0,2);
	G.addEdge(1,3);
	G.addEdge(2,3);
	G.addEdge(1,4);
	G.addEdge(3,5);
	G.addEdge(4,5);
	G.addEdge(3,4);

	Partition part(6);
	part.allToOnePartition();
	part.setUpperBound(2);
	index secondPartition = 1;
	part.moveToSubset(secondPartition, 4);
	part.moveToSubset(secondPartition, 5);

	auto map = part.subsetSizeMap();
	EXPECT_EQ(4, map.at(0));
	EXPECT_EQ(2, map.at(1));

	EXPECT_EQ(3, part.calculateCutWeight(G));

	EXPECT_EQ(0, Partitioner::calculateGain(G, part, 3, secondPartition));
	EXPECT_EQ(-1, Partitioner::calculateGain(G, part, 1, secondPartition));
}

TEST_F(PartitionerGTest, testGainRandomGraph) {
	const count k = 10;
	Graph H = HyperbolicGenerator(1000, 6).generate();
	ClusteringGenerator clusterGen;
	Partition part = clusterGen.makeRandomClustering(H, k);
	edgeweight cut = part.calculateCutWeight(H);

	for (index i = 0; i < 100; i++) {
		index v = H.randomNode();
		index targetPart = Aux::Random::index(k);
		edgeweight gain = Partitioner::calculateGain(H, part, v, targetPart);
		part.moveToSubset(targetPart, v);
		edgeweight newcut = part.calculateCutWeight(H);
		EXPECT_EQ(cut, newcut + gain);
		cut = newcut;
	}
}

TEST_F(PartitionerGTest, testPartitionerOnRealGraphWithChargedNodes) {
	/**
	 * read in Graph
	 */

	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count n = G.numberOfNodes();
	const count targetK = 10;
	const count numCharged = 5;

	ASSERT_GE(n, targetK);
	ASSERT_GE(targetK, numCharged);

	/**
	 * prepare random charged nodes
	 */
	std::vector<index> charged(numCharged);
	for (index i = 0; i < numCharged; i++) {
		bool present;
		index chIndex;
		do {
			/**
			 * sample random index. If already present, sample again.
			 */
			chIndex = Aux::Random::index(n);
			present = false;
			for (index j = 0; j < i; j++) {
				if (charged[j] == chIndex) {
					present = true;
				}
			}
		} while (present);
		charged[i] = chIndex;
	}

	/**
	 * call partitioner
	 */
	Partitioner part(G, targetK, 10, false, charged);
	part.run();
	Partition result = part.getPartition();

	/**
	 * check if charged nodes are in different partitions
	 */
	for (index i = 0; i < numCharged; i++) {
		for (index j = 0; j < i; j++) {
			EXPECT_NE(result[i], result[j]);
		}
	}
}

TEST_F(PartitionerGTest, testPartitionerOnBarabasiAlbert) {
	BarabasiAlbertGenerator gen(10, 5000, 10);
	Graph G =  gen.generate();

	const count targetK = 10;

	count n = G.numberOfNodes();
	Partitioner part(G, targetK);
	part.run();
	Partition result = part.getPartition();
	DEBUG("Resulted in ", result.numberOfSubsets(), " partitions, with a cut of weight ", result.calculateCutWeight(G));
	auto map = result.subsetSizeMap();

	double averageSize = n / targetK;
	double variance = 0;
	double maxImbalance = 0;

	for (auto entry : map) {
		DEBUG("Partition ", entry.first, " of size ", entry.second);
		variance += std::abs(entry.second - averageSize);
		double imbalance = std::abs(entry.second - averageSize) / averageSize;
		if (maxImbalance < imbalance) maxImbalance = imbalance;
	}
	variance += (targetK - result.numberOfSubsets())*averageSize; //missing partitions have a cardinality of zero
	if (targetK > result.numberOfSubsets() && maxImbalance < 1) maxImbalance = 1;
	DEBUG("Variance of sizes is ", variance / targetK, " imbalance is ", maxImbalance);

}

TEST_F(PartitionerGTest, testPartitionerOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count targetK = 10;

	count n = G.numberOfNodes();
	Partitioner part(G, targetK);
	part.run();
	Partition result = part.getPartition();
	DEBUG("Resulted in ", result.numberOfSubsets(), " partitions, with a cut of weight ", result.calculateCutWeight(G));
	auto map = result.subsetSizeMap();

	double averageSize = n / targetK;
	double variance = 0;
	double maxImbalance = 0;

	for (auto entry : map) {
		DEBUG("Partition ", entry.first, " of size ", entry.second);
		variance += std::abs(entry.second - averageSize);
		double imbalance = std::abs(entry.second - averageSize) / averageSize;
		if (maxImbalance < imbalance) maxImbalance = imbalance;
	}
	variance += (targetK - result.numberOfSubsets())*averageSize; //missing partitions have a cardinality of zero
	if (targetK > result.numberOfSubsets() && maxImbalance < 1) maxImbalance = 1;
	DEBUG("Variance of sizes is ", variance / targetK, " imbalance is ", maxImbalance);
}

}//end namespace NetworKit
