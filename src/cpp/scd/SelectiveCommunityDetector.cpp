/*
 * SelectiveCommunityDetector.cpp
 *
 *  Created on: 15.05.2013
 *      Author: cls, Yassine Marrakchi
 */

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

/**
 * @param[in] G	pointer to the considered graph
 */
SelectiveCommunityDetector::SelectiveCommunityDetector(const Graph& G) : G(G) {
}

SelectiveCommunityDetector::~SelectiveCommunityDetector() {
	// TODO Auto-generated destructor stub
}

std::map<node, std::set<node> > SelectiveCommunityDetector::getResult() {
	return result;

}

std::map<node, double> SelectiveCommunityDetector::getTimings() {
	return timings;
}

} /* namespace NetworKit */
