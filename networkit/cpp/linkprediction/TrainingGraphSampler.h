/*
 * TrainingGraphSampler.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef TRAININGGRAPHSAMPLER_H_
#define TRAININGGRAPHSAMPLER_H_

#include <utility>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Samples a given amount of links from a given graph in order to generate a training graph.
 */
class TrainingGraphSampler {
public:
  /**
   * Returns a graph that contains @a trainPercentage percent of links form the given graph @a G.
   * The links are randomly selected from G until the given percentage is reached.
   * @param G The graph to construct the training graph from
   * @param trainPercentage Percentage of links regarding the number of links in the
   * given graph that should be in the returned graph
   * @return a graph that contains the given percentage of links from @a G
   */
  static Graph byPercentage(const Graph& G, double trainPercentage);

  /**
   * Returns a graph that contains @a numTrainEdges links from the given graph @a G.
   * The links are randomly selected from G until the given count is reached.
   * @param G The graph to construct the training graph from
   * @param numTrainLinks Number of links the returned graph should consist of
   * @return a graph that contains the given number of links from @a G
   */
  static Graph byCount(const Graph& G, count numTrainLinks);

};

} // namespace NetworKit

#endif /* TRAININGGRAPHSAMPLER_H_ */