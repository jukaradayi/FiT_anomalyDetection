/*
 *  history_graph.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
 */

#ifndef STREAMGRAPHS_HISTORY_GRAPH_HPP_
#define STREAMGRAPHS_HISTORY_GRAPH_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>

#include "globals.hpp"
#include <limits>
//#include <boost>
#include <utility>
#include <cstdint>
#include <map>
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
//#include <string_view>

using namespace NetworKit;
namespace StreamGraphs {
//using Time = uint64_t;
//using Bound = uint64_t;
//using Count = uint64_t;
//
//constexpr int none = std::numeric_limits<uint64_t>::max();
//// (t,u,v) triplet
//struct Interaction {
//    Time t;
//    node u, v;
//    
//    Interaction() : t(none), u(none), v(none) {}
//
//    Interaction(Time _t, node _u, node _v) {
//        t = _t;
//        //u = std::min(_u, _v);
//        //v = std::max(_u, _v);
//        u = _u;
//        v = _v;
//    }
//};
//struct Edge {
//    node u, v;
//
//    Edge() : u(none), v(none) {}
//
//    Edge(node _u, node _v, bool sorted = false) {
//        u = sorted ? std::min(_u, _v) : _u;
//        v = sorted ? std::max(_u, _v) : _v;
//    }
//    bool operator <( const Edge &rhs ) const
//    {
//       if (u < rhs.u) {
//            return true;
//       } else if (u == rhs.u && v < rhs.v) {
//           return true;
//       } else {
//           return false;
//       }
//    }
//};

class HistoryGraph { // TODO check who is private or public
public:
    // main graph
    NetworKit::Graph main_graph = NetworKit::Graph(0,true,false);

    // projection graph
    NetworKit::Graph top_graph = NetworKit::Graph(0,true,false);
    NetworKit::Graph bot_graph = NetworKit::Graph(0,true,false);


    // unpacked graph
    NetworKit::Graph unpk_graph = NetworKit::Graph(0,true,false);

    // degree limits
    Bound main_bound;
    Bound proj_bound;

    // number of nodes
    Count N;

    // enable or disable some graph usage
    bool use_projection;
    bool use_unpacked;

    // node names mappings
    std::vector<uint64_t> node2main; // main graph node mapping
    std::vector<uint64_t> main2node;
    std::vector<uint64_t> node2top; // projection graph node mapping
    std::vector<uint64_t> top2node;
    std::vector<uint64_t> node2bot; // projection graph node mapping
    std::vector<uint64_t> bot2node;
    std::vector<uint64_t> node2unpk; // unpacked graph node mapping
    std::vector<uint64_t> unpk2node;

    // queue and storage
    std::vector<Interaction> queue;
    std::map<Edge, uint64_t> counter;
    std::vector<node> removed_main; // nodes from main graph to be restored
    std::vector<node> removed_top; // nodes from projection graph to be restored
    std::vector<node> removed_bot; // nodes from projection graph to be restored
    std::vector<node> removed_unpk; // nodes from unpacked graph to be restored

    // degree distribution
    std::vector<Count> main_degree_distribution;
    std::vector<Count> top_degree_distribution;
    std::vector<Count> bot_degree_distribution;
    std::vector<Count> unpk_degree_distribution;

    // weighted degree sequence
    std::vector<Count> main_weightedDegree_sequence;
    std::vector<Count> top_weightedDegree_sequence;
    std::vector<Count> bot_weightedDegree_sequence;
    std::vector<Count> unpk_weightedDegree_sequence;


    HistoryGraph(bool use_projection, bool use_unpacked, Bound main_bound, Bound proj_bound, Count N);
    ~HistoryGraph();
    //void HistoryGraph::restoreNode(node u);

    void updateGraph(Interaction i);

    node addNode(node _u, bool is_top);

    void removeNode(node _u, bool is_top);

    void removeEdgeProjection(NetworKit::Graph& proj_graph, node node_main, node proj_node, bool is_top);


    void addEdgeProjection(node node_main, node proj_node, bool is_top);

    inline virtual void trimQueue(Time t){};

    inline virtual Bound getWindow(){};
};

}
#endif
