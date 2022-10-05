/**
 *  Implementation of dynamic history graphs.
 *  Given a "link stream" in input, with the following format:
 *    
 *    t_0 u_0 v_0
 *    t_1 u_1 v_1
 *    ...
 *    t_n u_n v_n
 *
 *  a history graph H_s(t_n) (or G_d(t_n)) is defined as a static graph
 *  containing link (u_n, v_n) and the links preceding this one.
 *  The conditions to include or exclude links from H_s(t_n) and G_d(t_n)
 *  are defined in H_graph.cpp and G_graph.cpp.
 *
 *
 *
 */

#ifndef STREAMGRAPHS_HISTORY_GRAPH_HPP_
#define STREAMGRAPHS_HISTORY_GRAPH_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>

#include "SortedCounter.hpp"
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

class HistoryGraph {
public:

    NetworKit::Graph& main_graph; ///< the main graph, stored as a NetworKit::Graph

    // projection graph
    NetworKit::Graph& top_graph ; ///< the top projection graph, stored as a NetworKit::Graph
    NetworKit::Graph& bot_graph ; ///< the bottom projection graph, stored as a NetworKit::Graph

    // unpacked graph
    //NetworKit::Graph& unpk_graph;// = NetworKit::Graph(0,true,false);

    // degree limits
    const Bound main_bound; ///< the degree bound to be used on the main graph
    const Bound proj_bound; ///< the degree bound to be used in the projections graph

    // number of nodes
    Count N; ///< the number of nodes in the main graph

    // enable or disable some graph usage
    const bool use_projection; ///< a boolean, used to enable/disable the use of projection graphs
    const bool use_unpacked; ///< unused TODO REMOVE
    const bool is_bipartite; ///< a boolean, enable if the input graph is bipartite

    // node names mappings // TODO might define them as arrays
    std::vector<uint64_t> node2main; ///< mapping from input data node names to main_graph node names
    std::vector<uint64_t> main2node; ///< mapping from main_graph node names to input data node names
    std::vector<uint64_t> node2top; ///< mapping from input data node names to top_graph node names
    std::vector<uint64_t> top2node; ///< mapping from top_graph node names to input data node names
    std::vector<uint64_t> node2bot; ///< mapping from input data node names to bot_graph node names
    std::vector<uint64_t> bot2node; ///< mapping from bot_graph node names to input data node names
    std::vector<uint64_t> node2unpk; 
    std::vector<uint64_t> unpk2node;

    //uint64_t* node2main; // main graph node mapping
    //uint64_t* main2node;
    //uint64_t* node2top; // projection graph node mapping
    //uint64_t* top2node;
    //uint64_t* node2bot; // projection graph node mapping
    //uint64_t* bot2node;
    //uint64_t* node2unpk; // unpacked graph node mapping
    //uint64_t* unpk2node;


    // queue and storage
    //std::vector<Interaction> queue;
    std::queue<Interaction> queue; ///< queue storing the interactions (t,u,v) used to build main_graph
    std::map<Edge, uint64_t> counter; ///< mapping from an edge to its weight in main_graph 
    std::vector<node> removed_main; ///< vector storing the nodes that have been used and were removed from main_graph. When adding a new node, we "resurect" those removed node in priority, insted creating a new node
    std::vector<node> removed_top; ///< same as removed_main for top_graph
    std::vector<node> removed_bot; ///< same as removed_main for bot_graph
    std::vector<node> removed_unpk; //
    
    SortedCounters<node> degree_counter;
    SortedCounters<node> weightedDegree_counter;
    SortedCounters<Edge> weight_counter;

    // degree distribution
    std::vector<Count> main_degree_distribution; ///< the degree distribution of the main graph
    std::vector<Count> top_degree_distribution; ///< the degree distribution of the top graph
    std::vector<Count> bot_degree_distribution; ///< the degree distribution of the bottom graph
    std::vector<Count> unpk_degree_distribution;

    // weighted degree sequence
    std::vector<Count> main_weightedDegree_sequence; ///< the weighted degree sequence of the main graph
    std::vector<Count> top_weightedDegree_sequence; ///< the weighted degree sequence of the top graph
    std::vector<Count> bot_weightedDegree_sequence; ///< the weighted degree sequence of the bottom graph
    std::vector<Count> unpk_weightedDegree_sequence;

    /** Maintain a history graph and its top and bottom projections (if enabled)
     *
     * @brief history graph constructor
     *
     * @param main_graph (NetworKit::Graph) The main graph. 
     * @param top_graph (NetworKit::Graph) The top projection graph
     * @param bot_graph (NetworKit::Graph) the bottom projection graph
     * @param use_projection (bool) Set to True when you want to use the top and bottom projection graphs
     * @param is_bipartite (bool)  Set to True when the main graph is bipartite. When set to False and use_projection == True, only one projection graph is used (top_graph).
     * @param main_bound (int) The degree bound to be used on the main_graph. When accessing the neighborhood of a node, it will only iterate on nodes with degree lesser than main_bound
     * @param proj_bound (int) The degree bound to be used on the projections.
     * @param N (int) The number of nodes in the link stream.
     *
     */
    HistoryGraph(NetworKit::Graph& main_graph, NetworKit::Graph& top_graph, NetworKit::Graph& bot_graph, const bool use_projection, const bool use_unpacked, const bool is_bipartite, const Bound main_bound, const Bound proj_bound, Count N);
    ~HistoryGraph();

    /**
     * Given a new interaction (t,u,v), add it to the queue, and add the 
     * edge (u,v) to the main_graph, and update the top and bottom projections accordingly.
     *
     * Note: when updating the projection graph top_graph (resp. bot_graph), if the degree of 
     * the node v_0 in main_graph (resp. u_0) was previously d(v_0) == main_bound but is now 
     * d(v_0) > main_bound (resp d(u_0) ), the influence of v_0 on top_graph 
     * (resp u_0 on bot_graph) was taken into account but now is not anymore,
     * so loop through each neighbor n of v_0 (resp u_0) in main_graph and 
     * remove the edge (u_0, n) in top_graph (resp. (v_0, n) in bot_graph).
     *

     * @brief Add a new interaction to the main graph and its projections.
     * 
     *
     * @param i (Interaction) The interaction (t, u, v) to add to the main graph.
     */
    void updateGraph(const Interaction i);

    /**
     * Add a new node to the main_graph and the projections.
     * If removed_main is not empty, the name of the added node will be n=removed_main.pop_back(), else it will be n=N+1 where N is the biggest node ID in the graph.
     * @brief Add a new node to the graph.
     * @see removed_main
     * @param _u (node) The node to add to the main graph and its projections
     * @param is_top (bool) Set to true if the node should be added to the top projection, false if it should be added to the bottom projection
     *
     */
    node addNode(const node _u, const bool is_top);

    /**
     * Remove a node from the main_graph and the projections. When removing a node from the main_graph, add its name to removed_main.
     * @brief Remove a node from the graph.
     * @see removed_main
     * @param _u (node) The node to remove from the main graph and its projections
     * @param is_top (bool) Set to true if the node should be removed from the top projection, false if it should be removed from the bottom projection
     *
     */
    void removeNode(const node _u, const bool is_top);

    /**
     * Remove edges from one of the projection graph
     * @param proj_graph (NetworKit::Graph) The projection from which the edge should be removed
     * @param proj_node (node) The first node of the edges to be removed
     * @param node_main (node) For each neighbor n of node_main in main_graph, remove edge (proj_node, n) from proj_graph
     *
     */
    void removeEdgeProjection(NetworKit::Graph& proj_graph,const node node_main,const node proj_node,const bool is_top);

    /**
     * Add edges from one of the projection graph
     * @param proj_graph (NetworKit::Graph) the projection graph to which we should add edges
     * @param proj_node (node) The first node of the edges to be added
     * @param node_main (node) For each neighbor n of node_main in main_graph, add edge (proj_node, n) to proj_graph
     *
     */ 
    void addEdgeProjection(const node node_main,const node proj_node,const bool is_top);

    /**
     * Virtual function
     * @brief Remove edges from main_graph, top_graph and bottom_graph that don't respect the time condition anymore
     * See H_graph.cpp and G_graph.ccp for implementations
     */
    inline virtual void trimQueue(const Time t){};

    /**
     * @brief Return the window size.
     */
    inline virtual Bound getWindow(){};

    /**
     * Increase the total weight of main_graph by 1
     */
    inline void increaseTotalWeight(){total_weight+=1;};

    /**
     * Decrease the total weight of main_graph by 1
     */
    inline void decreaseTotalWeight(){total_weight-=1;};

    /**
     * Return the total weight of main_graph
     */
    inline Count getTotalWeight(){return total_weight;};

    /**
     * Template. Iterate on the neighbors n of u in graph, and apply the function handle to
     * n if the degree of n is less than main_bound
     * @param graph (NetworKit::Graph) The graph on which we iterate
     * @param u (node) Iterate on the neighborhood of u
     * @param handle The function to be applied to the neighbors n of u
     */
    template <typename L>
    void forBoundedNeighbors(NetworKit::Graph& graph, node u, L handle) const;

private:
    Count total_weight = 0; ///< Integer, the total weight of the main_graph
    Count redundancy_top; ///< Integer, the redundancy of the top graph projection
    Count redundancy_bot; ///< Integer, the redundancy of the bottom graph projection


};

template <typename L>
void HistoryGraph::forBoundedNeighbors(NetworKit::Graph& graph, node u, L handle) const {

    /**
     * Iterate on the neighbors of u, and if the degree is lesser than main_bound
     *  apply handle to it.
     */
    for (NetworKit::Graph::NeighborIterator bN = graph.neighborRange(u).begin();
                bN != graph.neighborRange(u).end(); ++bN) {
        if (graph.degree(*bN) < main_bound) {
            handle(*bN);
        }


    }
}
}
#endif
