/*
 * MultiLevelPartitionerGTest.cpp
 *
 *  Created on: 15.10.2015
 *      Author: moritzl
 */

#include "MultiLevelPartitionerGTest.h"

#include "../MultiLevelPartitioner.h"
#include "../../graph/Graph.h"
#include "../../Globals.h"
#include "../../generators/BarabasiAlbertGenerator.h"
#include "../../generators/HyperbolicGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../../io/PartitionReader.h"
#include "../../community/ClusteringGenerator.h"

namespace NetworKit {

TEST_F(MultiLevelPartitionerGTest, testGainUnweighted) {
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

	EXPECT_EQ(0, MultiLevelPartitioner::calculateGain(G, part, 3, secondPartition));
	EXPECT_EQ(-1, MultiLevelPartitioner::calculateGain(G, part, 1, secondPartition));
}

TEST_F(MultiLevelPartitionerGTest, testGainRandomGraph) {
	const count k = 10;
	Graph H = HyperbolicGenerator(1000, 6).generate();
	ClusteringGenerator clusterGen;
	Partition part = clusterGen.makeRandomClustering(H, k);
	edgeweight cut = part.calculateCutWeight(H);

	for (index i = 0; i < 100; i++) {
		index v = H.randomNode();
		index targetPart = Aux::Random::index(k);
		edgeweight gain = MultiLevelPartitioner::calculateGain(H, part, v, targetPart);
		part.moveToSubset(targetPart, v);
		edgeweight newcut = part.calculateCutWeight(H);
		EXPECT_EQ(cut, newcut + gain);
		cut = newcut;
	}
}

TEST_F(MultiLevelPartitionerGTest, testFMOnContractedBarabasiAlbert) {
	METISGraphReader reader;
	Graph G = reader.read("input/intermediate.graph");
	const count n = G.numberOfNodes();

	PartitionReader partReader;
	Partition part = partReader.read("input/intermediate.part");
	const count k = part.numberOfSubsets();
	ASSERT_EQ(n, part.numberOfElements());

	edgeweight gain;
	do {
		gain = MultiLevelPartitioner::fiducciaMatheysesStep(G, part);
		DEBUG("Found gain ", gain, " in FM-step with ", G.numberOfNodes(), " nodes and ", part.numberOfSubsets(), " partitions.");
	} while (gain > 0);
	EXPECT_EQ(k, part.numberOfSubsets());
}

TEST_F(MultiLevelPartitionerGTest, testFMOnRealGraphAndRandomPartition) {
	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count k = 10;

	ClusteringGenerator clusterGen;
	Partition part = clusterGen.makeRandomClustering(G, k);
	edgeweight cut = part.calculateCutWeight(G);

	edgeweight gainsum = 0;
	edgeweight gain;
	do {
		gain = MultiLevelPartitioner::fiducciaMatheysesStep(G, part);
		gainsum += gain;
		DEBUG("Found gain ", gain, " in FM-step with ", G.numberOfNodes(), " nodes and ", part.numberOfSubsets(), " partitions.");
	} while (gain > 0);
	EXPECT_EQ(0, gain);
	edgeweight newcut = part.calculateCutWeight(G);
	EXPECT_NEAR(cut, newcut + gainsum, 0.001);
	DEBUG("Final cut weight: ", newcut);
	EXPECT_EQ(k, part.numberOfSubsets());
}

TEST_F(MultiLevelPartitionerGTest, testRegionGrowingOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count k = 10;
	const count n = G.numberOfNodes();

	/**
	 * prepare random starting nodes
	 */
	std::vector<index> startingNodes(k);
	for (index i = 0; i < k; i++) {
		bool present;
		index chIndex;
		do {
			/**
			 * sample random index. If already present, sample again.
			 */
			chIndex = Aux::Random::index(n);
			present = false;
			for (index j = 0; j < i; j++) {
				if (startingNodes[j] == chIndex) {
					present = true;
				}
			}
		} while (present);
		startingNodes[i] = chIndex;
	}

	Partition part = MultiLevelPartitioner::growRegions(G, startingNodes);
	EXPECT_EQ(k, part.numberOfSubsets());
}

TEST_F(MultiLevelPartitionerGTest, testPartitionerOnRealGraphWithChargedNodes) {
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
	MultiLevelPartitioner part(G, targetK, 10, false, charged);
	part.run();
	Partition result = part.getPartition();

	EXPECT_EQ(targetK, result.numberOfSubsets());

	/**
	 * check if charged nodes are in different partitions
	 */
	for (index i = 0; i < numCharged; i++) {
		for (index j = 0; j < i; j++) {
			EXPECT_NE(result[charged[i]], result[charged[j]]);
		}
	}
}

TEST_F(MultiLevelPartitionerGTest, testPartitionerOnBarabasiAlbert) {
	/**
	 * generate graph
	 */
	BarabasiAlbertGenerator gen(10, 2000, 10);
	Graph G =  gen.generate();

	const count targetK = 10;
	const double maxImbalance = 0.1;

	/**
	 * call partitioner
	 */
	MultiLevelPartitioner part(G, targetK, maxImbalance);
	part.run();
	Partition result = part.getPartition();

	/**
	 * check results
	 */
	DEBUG("Resulted in ", result.numberOfSubsets(), " partitions, with a cut of weight ", result.calculateCutWeight(G));
	EXPECT_EQ(targetK, result.numberOfSubsets());

	auto map = result.subsetSizeMap();
	double measuredImbalance = result.getImbalance(targetK);

	EXPECT_LE(measuredImbalance, maxImbalance);

}

TEST_F(MultiLevelPartitionerGTest, testPartitionerOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count targetK = 10;

	const count n = G.numberOfNodes();
	const double maxImbalance = 0.1;

	MultiLevelPartitioner part(G, targetK, maxImbalance);
	part.run();
	Partition result = part.getPartition();
	EXPECT_EQ(n, result.numberOfElements());
	edgeweight cutWeight = result.calculateCutWeight(G);
	DEBUG("Resulted in ", result.numberOfSubsets(), " partitions, with a cut of weight ", cutWeight);
	EXPECT_EQ(targetK, result.numberOfSubsets());

	double measuredImbalance = result.getImbalance(targetK);

	EXPECT_LE(measuredImbalance, maxImbalance);

	/**
	 * compare with naive partition
	 */
	ClusteringGenerator gen;
	Partition naive = gen.makeContinuousBalancedClustering(G, targetK);

	ASSERT_EQ(targetK, naive.numberOfSubsets());
	edgeweight naivecut = naive.calculateCutWeight(G);
	EXPECT_LT(cutWeight, naivecut);
}

TEST_F(MultiLevelPartitionerGTest, testPartitionerNaiveComparisonRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/bacteriorhodopsin-10-2.5.graph");

	const count targetK = 10;

	const count n = G.numberOfNodes();
	const double maxImbalance = 0.1;

	MultiLevelPartitioner part(G, targetK, maxImbalance);
	part.run();
	Partition result = part.getPartition();
	edgeweight cutWeight = result.calculateCutWeight(G);
	DEBUG("Resulted in ", result.numberOfSubsets(), " partitions, with a cut of weight ", cutWeight);
	EXPECT_EQ(targetK, result.numberOfSubsets());

	/**
	 * compare with naive partition and refinement step
	 */
	ClusteringGenerator gen;
	Partition naive = gen.makeContinuousBalancedClustering(G, targetK);

	edgeweight gain;
	do {
		gain = MultiLevelPartitioner::fiducciaMatheysesStep(G, naive);
		assert(gain == gain);
		DEBUG("Found gain ", gain, " in FM-step with ", G.numberOfNodes(), " nodes and ", naive.numberOfSubsets(), " partitions.");
	} while (gain > 0);

	edgeweight cutNaiveFM = naive.calculateCutWeight(G);
	EXPECT_LE(cutWeight, cutNaiveFM);
}

}//end namespace NetworKit
