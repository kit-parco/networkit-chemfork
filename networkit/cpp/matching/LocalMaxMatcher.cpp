/*
 * LocalMaxMatcher.cpp
 *
 *  Created on: 05.12.2012
 */

#include "LocalMaxMatcher.h"

namespace NetworKit {

LocalMaxMatcher::LocalMaxMatcher(const Graph& graph, const std::vector<index>& chargedVertices, const std::vector<std::pair<index,index> >& forbiddenEdges, const std::vector<double>& nodeWeights, bool useConductance): Matcher(graph), chargedVertices(chargedVertices), forbiddenEdges(forbiddenEdges), nodeWeights(nodeWeights), useConductance(useConductance)
{
	if (graph.isDirected()) throw std::runtime_error("Matcher only defined for undirected graphs");
	if (this->nodeWeights.size() == 0) {
		this->nodeWeights.resize(graph.numberOfNodes(), 1);
	}
}

// TODO: update to new edge attribute system
// TODO: make local max matching parallel


void LocalMaxMatcher::run() {
	int64_t z = G.upperNodeIdBound();
	count E = G.numberOfEdges();

	// put edges into array of triples
	struct MyEdge {
		node s; // source
		node t; // target
		edgeweight w; // weight
	};

	std::vector<MyEdge> edges(E);
	index e = 0;
	G.forEdges([&](node u, node v, edgeweight w) {
		//we remove forbidden edge, since they should not be matched.
		if (std::none_of(forbiddenEdges.begin(), forbiddenEdges.end(), [u,v](std::pair<index,index> edge){return (edge.first == u && edge.second == v) || (edge.first == v && edge.second == u);})) {
			edges[e].s = u;
			edges[e].t = v;
			edges[e].w = w + Aux::Random::real(1e-6);
			++e;
		} else {
			DEBUG("Edge (", u, ",", v, ") was forbidden, not matched.");
		}
	});

	edges.resize(e);

	// candidates records mating candidates
	std::vector<MyEdge> candidates(z);
	G.parallelForNodes([&](node u) {
		candidates[u].w = (edgeweight) 0;
		candidates[u].t = u; // itself as mating partner => unmatched
	});

	//note charged nodes
	std::vector<bool> charged(z, false);
	for (index c : chargedVertices) {
		charged[c] = true;
	}

	auto conductance = [&](MyEdge edge) -> double {return edge.w / (nodeWeights[edge.s] * nodeWeights[edge.t]);};

	while (E > 0) {
		// for each edge find out if it is locally maximum and allowed
		for (auto edge: edges) {
			const double edgeRating = useConductance ? conductance(edge) : edge.w;
			if (edgeRating > candidates[edge.s].w && edgeRating > candidates[edge.t].w && edge.s != edge.t && (!charged[edge.s] || !charged[edge.t]) ) {
				candidates[edge.s].t = edge.t;
				candidates[edge.s].w = edgeRating;
				candidates[edge.t].t = edge.s;
				candidates[edge.t].w = edgeRating;
			}
		}

		// check if candidates agree to match; if so, then match them
		for (auto edge: edges) {
			node u = edge.s;
			node v = edge.t;
			if (candidates[u].t == v && candidates[v].t == u && u != v && (!charged[u] || !charged[v])) {
				// both nodes agree
				M.match(u, v);
			}
		}

		// create remaining "graph" by selecting remaining edges (as triples)
		// adjust candidates
		std::vector<MyEdge> newEdges;
		for (auto edge: edges) {
			if (! M.isMatched(edge.s) && ! M.isMatched(edge.t) && edge.s != edge.t && (!charged[edge.s] || !charged[edge.t])) {
				newEdges.push_back(edge);
				candidates[edge.s].w = (edgeweight) 0;
				candidates[edge.t].w = (edgeweight) 0;
			}
		}
		edges = newEdges;
		E = edges.size();
	}

}

} /* namespace NetworKit */
