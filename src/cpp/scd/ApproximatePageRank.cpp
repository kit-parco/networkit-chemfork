/*
 * ApproximatePageRank.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "ApproximatePageRank.h"

namespace NetworKit {

ApproximatePageRank::ApproximatePageRank(Graph& g, double alpha_, double epsilon):
		G(g), alpha(alpha_), oneMinusAlphaOver2((1.0 - alpha) * 0.5), eps(epsilon)
{

}

void ApproximatePageRank::push(node u, node seed, std::vector<double>& pr, std::vector<double>& residual, std::queue<node>& activeNodes)
{
	double res = residual[u];
	pr[u] = pr[u] + alpha * res;
//	TRACE("normalizedResid[", u, "]: ", normalizedResid[u]);
	double mass = oneMinusAlphaOver2 * res / G.degree(u);

	G.forNeighborsOf(u, [&](node v) {
		residual[v] = residual[v] + mass;
		normalizedResid[v] = residual[v] / G.degree(v);
		if (normalizedResid[v] >= eps) {
			activeNodes.push(v);
		}
//		TRACE("normalizedResid[", v, "]: ", normalizedResid[v]);
	});

	residual[u] = oneMinusAlphaOver2 * res;
	normalizedResid[u] = residual[u] / G.degree(u);
	if (normalizedResid[u] >= eps) {
		activeNodes.push(u);
	}

//	TRACE("residual[", u, "]: ", residual[u]);
}

ApproximatePageRank::~ApproximatePageRank() {

}

std::vector<double> ApproximatePageRank::run(node seed) {
	// initialize vectors: pr = 0, residual = characteristic(seed)
	count n = G.numberOfNodes();
	std::vector<double> pr(n, 0.0);
	normalizedResid = pr;
	std::vector<double> residual = pr;
	residual[seed] = 1.0;
	normalizedResid[seed] = 1.0 / (double) G.degree(seed);
	std::queue<node> activeNodes;
	activeNodes.push(seed);


	while (activeNodes.size() > 0) {
		node v = activeNodes.front();
		activeNodes.pop();
		push(v, seed, pr, residual, activeNodes);
	}

//	auto converged([&](node& argmax) {
//		// TODO: accelerate
//		std::vector<double>::iterator max_elem = std::max_element(normalizedResid.begin(), normalizedResid.end());
//		argmax = std::distance(normalizedResid.begin(), max_elem);
////		TRACE("argmax: ", argmax, ", max: ", (* max_elem), ", eps: ", eps);
//		return ((* max_elem) < eps);
//	});
//
//
//	node argMax = seed;
//	while (! converged(argMax)) {
//		push(argMax, seed, pr, residual);
//	}

	return pr;
}

} /* namespace NetworKit */
