/*
 * Implementation of dynamic history graphs.
 * Given a number of edges s, HGraph is defined as the aggregated graph
 * of the input link-stream, with the s last edges.
 *  Given a "link stream" L in input, with the following format, with 
 *  set of nodes V and set of edges E:
 *    
 *    t_0 u_0 v_0
 *    t_1 u_1 v_1
 *    ...
 *    t_n u_n v_n
 *
 * we note Q_s(t_p) = [(t_(p-s), u_(p-s), v_(p-s)), ... (t_p, u_p, v_p)].
 * the HGraph, for p in [0,n], H_d(t_p) = (E_H_d = {(u_i, v_i) for i in [p-s,p]} , V_H_d = {u such that u in V and there exist v in V such that (u,v) in E_H_d}).
 *
 * The HGraph class inherits from HistoryGraph and implements the "trimQueue" method.
 * Note : the HGraph are valued graphs, the weight of an edge is the number of occurence of this edge in the window.

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
     * @param window (int) When a new interaction (t_n, u_n, v_n) is added to main_graph, an interaction (t_0, u_0, v_0) is removed from the queue and the edge (u_0, v_0) has its weight reduced in main_graph if (t_n - t_0) > window.
     */
    HGraph(NetworKit::Graph& main_graph, NetworKit::Graph& top_graph, NetworKit::Graph& bot_graph, const bool use_projection, const bool use_unpacked, const bool is_bipartite, const Bound main_bound, const Bound proj_bound, Count N, const Bound window); 


    ~HGraph();

    //void updateGraph(Interaction i);
    /** Remove edges from main_graph and its projections graph if they don't satisfy the window
     * condition anymore.
     * When a new interaction (t_n, u_n, v_n) is added to main_graph, an interaction (t_d, u_d, v_d) is removed from the queue and the edge (u_d, v_d) has its weight reduced in main_graph if (n -d  > window.
     * If the projections are enabled, when removing edge (u_0, v_0) from the graph, 
     * for each neighbor n of v_0 in main_graph, decrement the weight of edge (u_0, n) in top_graph, and for 
     * each neighbor n' of u_0 in main_graph, decrement the weight of edge (v_0, n') from bot_graph.
     * Note: when updating the projection graph top_graph (resp. bot_graph), if the degree of 
     * the node v_0 in main_graph (resp. u_0) was previously d(v_0) > main_bound but is now 
     * d(v_0) == main_bound (resp d(u_0) ), the influence of v_0 on top_graph (resp u_0 on bot_graph) was not taken into account but now is, so loop through each neighbor n of v_0 (resp u_0) in main_graph and add the edge (u_0, n) in top_graph (resp. (v_0, n) in bot_graph).
     *
     *
     * @brief Remove edges from graph when they don't satisfy the window condition.
     *
     */ 
    void trimQueue(Time t);

    /**
     *  Return the window.
     */
    inline Bound getWindow() {return window;};
//private:

    // window size
    Bound window;

    
};
}
#endif
