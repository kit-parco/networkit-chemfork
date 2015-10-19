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
#include "../../io/METISGraphReader.h"

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

TEST_F(PartitionerGTest, testRegionGrowing) {

}

TEST_F(PartitionerGTest, testRecursiveBisection) {

}

TEST_F(PartitionerGTest, testChargedNodes) {

}

TEST_F(PartitionerGTest, testPartitioner) {
	//const count n = 10000;
	BarabasiAlbertGenerator gen(10, 10000, 10);
	Graph G =  gen.generate();

	const count targetK = 10;


	//Graph G = METISGraphReader().read("input/caidaRouterLevel.graph");
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
