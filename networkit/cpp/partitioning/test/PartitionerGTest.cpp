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

TEST_F(PartitionerGTest, testPartitioner) {
	BarabasiAlbertGenerator gen(10, 1000, 10);

	Graph G = gen.generate();
	Partitioner part(G, 10);
	part.run();
}

}//end namespace NetworKit
