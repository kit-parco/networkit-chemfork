
#include "UnionMaximumSpanningForest.h"
#include "../auxiliary/SignalHandling.h"

namespace NetworKit {

UnionMaximumSpanningForest::UnionMaximumSpanningForest(const Graph &G) : G(G), hasWeightedEdges(false), hasUMSF(false), hasAttribute(false) { };

template <typename A>
UnionMaximumSpanningForest::UnionMaximumSpanningForest(const Graph &G, const std::vector< A > &attribute) : G(G), hasWeightedEdges(false), hasUMSF(false), hasAttribute(false) {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("Error: Edges of G must be indexed for using edge attributes");
	}

	weightedEdges.reserve(G.numberOfEdges());

	G.forEdges([&](node u, node v, edgeid eid) {
		weightedEdges.emplace_back(u, v, attribute[eid], eid);
	});

	INFO(weightedEdges.size(), " weighted edges saved");

	hasWeightedEdges = true;
}

// instantiate for count and edgeweight
template UnionMaximumSpanningForest::UnionMaximumSpanningForest<edgeweight>(const Graph &G, const std::vector<edgeweight>&);
template UnionMaximumSpanningForest::UnionMaximumSpanningForest<count>(const Graph &G, const std::vector<count>&);

void UnionMaximumSpanningForest::run() {
	hasRun = false;
	hasUMSF = false;
	hasAttribute= false;

	Aux::SignalHandler handler;

	umsf = G.copyNodes();

	handler.assureRunning();

	bool useEdgeWeights = false;

	if (!hasWeightedEdges) {
		weightedEdges.reserve(G.numberOfEdges());

		G.forEdges([&](node u, node v, edgeweight weight, edgeid eid) {
			weightedEdges.emplace_back(u, v, weight, eid);
		});

		hasWeightedEdges = true;
		useEdgeWeights = true;
	}

	handler.assureRunning();

	bool calculateAttribute = false;

	if (G.hasEdgeIds()) {
		umsfAttribute.clear();
		umsfAttribute.resize(G.upperEdgeIdBound(), false);
		calculateAttribute = true;
	}

	std::sort(weightedEdges.begin(), weightedEdges.end(), std::greater<weightedEdge>());

	handler.assureRunning();

	edgeweight currentAttribute = std::numeric_limits<edgeweight>::max();

	std::vector<std::pair<node, node> > nodesToMerge;
	UnionFind uf(G.upperNodeIdBound());

	for (weightedEdge e : weightedEdges) {
		if (e.attribute != currentAttribute) {
			for (auto candidate : nodesToMerge) {
				uf.merge(candidate.first, candidate.second);
			}

			nodesToMerge.clear();
			currentAttribute = e.attribute;
		}


		if (uf.find(e.u) != uf.find(e.v)) {
			if (useEdgeWeights) {
				umsf.addEdge(e.u, e.v, e.attribute);
			} else {
				umsf.addEdge(e.u, e.v);
			}

			if (calculateAttribute) {
				umsfAttribute[e.eid] = true;
			}

			nodesToMerge.emplace_back(e.u, e.v);

		}
	}

	handler.assureRunning();

	hasUMSF = true;
	hasAttribute = calculateAttribute;
	hasRun = true;
}

bool UnionMaximumSpanningForest::inUMSF(edgeid eid) const {
	if (!hasAttribute) throw std::runtime_error("Error: Either the attribute hasn't be calculated yet or the graph has no edge ids.");

	return umsfAttribute[eid];
}

bool UnionMaximumSpanningForest::inUMSF(node u, node v) const {
	if (hasUMSF) {
		return umsf.hasEdge(u, v);
	} else if (hasAttribute) {
		return umsfAttribute[G.edgeId(u, v)];
	} else {
		throw std::runtime_error("Error: The run() method must be executed first");
	}
}

std::vector< bool > UnionMaximumSpanningForest::getAttribute(bool move) {
	std::vector<bool> result;

	if (!hasAttribute) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(umsfAttribute);
		hasAttribute = false;
	} else {
		result = umsfAttribute;
	}

	return result;
}

Graph UnionMaximumSpanningForest::getUMSF(bool move) {
	Graph result;

	if (!hasUMSF) throw std::runtime_error("Error: The run() method must be executed first");

	if (move) {
		result = std::move(umsf);
		hasUMSF = false;
	} else {
		result = umsf;
	}

	return result;
}

std::string UnionMaximumSpanningForest::toString() const {
	return "Union maximum-weight spanning forest";
}

bool UnionMaximumSpanningForest::isParallel() const {
	return false;
}


}