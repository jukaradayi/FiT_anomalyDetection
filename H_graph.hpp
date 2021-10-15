/*
 *  G_graph.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
 */

#ifndef STREAMGRAPHS_H_GRAPH_HPP_
#define STREAMGRAPHS_H_GRAPH_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>

#include <limits>
#include <utility>
#include <cstdint>
#include <map>
#include <string>
#include "history_graph.hpp"

using namespace NetworKit;
//using namespace StreamGraphs;
namespace StreamGraphs {

class HGraph final : public HistoryGraph {


    // constructor & destructor

public:
    HGraph(bool use_projection, bool use_unpacked, Bound main_bound, Bound proj_bound, Count N, Bound window); 


    ~HGraph();

    //void updateGraph(Interaction i);

    void trimQueue(Time t);

    inline Bound getWindow() {return window;};
//private:

    // window size
    Bound window;

    
};
}
#endif
