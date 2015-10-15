/*
 * PartitionerGTest.cpp
 *
 *  Created on: 15.10.2015
 *      Author: moritzl
 */

#include "PartitionerGTest.h"

#include "../Partitioner.h"
#include "../../graph/Graph.h"

using namespace NetworKit;

TEST_F(PartitionerGTest, testGain) {
	Graph G(6);
	G.addEdge(1,2);
	G.addEdge(1,2);
	G.addEdge(1,2);
	G.addEdge(1,2);
	G.addEdge(1,2);
	G.addEdge(1,2);

}

