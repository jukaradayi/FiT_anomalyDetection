/*
 *  G_graph.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
 */

#ifndef STREAMGRAPHS_G_GRAPH_HPP_
#define STREAMGRAPHS_G_GRAPH_HPP_

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

class GGraph final : public HistoryGraph {


    // constructor & destructor

public:
    GGraph(bool use_projection, bool use_unpacked, Bound main_bound, Bound proj_bound, Count N, Bound window); 


    ~GGraph();

    //void updateGraph(Interaction i);
    
    void trimQueue(Time t);
    std::vector<std::function<void()>> nctions;

    inline Bound getWindow(){return window;};
//private:

    // window size
    Bound window;


    
};
}
#endif
