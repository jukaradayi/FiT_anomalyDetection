// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "G_graph.hpp"
#include <networkit/graph/Graph.hpp>
#include <iostream>

namespace StreamGraphs {
using namespace StreamGraphs;
/** constructors **/
GGraph::GGraph(NetworKit::Graph& main_graph, NetworKit::Graph& top_graph, NetworKit::Graph& bot_graph, bool use_projection, bool use_unpacked, bool is_bipartite, Bound main_bound, Bound proj_bound, Count N, Bound window) : HistoryGraph(main_graph, top_graph, bot_graph,  use_projection, use_unpacked, is_bipartite, main_bound, proj_bound, N), window(window) {}

GGraph::~GGraph() {}

void GGraph::trimQueue(Time t) {

    auto decreaseMainDegree = [&](node u, node v) {
        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = main_graph.degree(u);
        Count v_degree = main_graph.degree(v);

        if (u_degree > 0) {
            main_degree_distribution[u_degree -1] += 1;
        }
        if (v_degree > 0) {
            main_degree_distribution[v_degree -1] += 1;
        }
        main_degree_distribution[u_degree] -= 1;
        main_degree_distribution[v_degree] -= 1;

    };

    while (t - queue.front().t > window) {
        // TODO change to .back()
        Interaction i0 = queue.front();
        queue.pop();

        node u0_main = node2main[i0.u];
        node v0_main = node2main[i0.v];
       
        // decrease edge counter
        Edge e0(i0.u, i0.v);
        Edge e0_main(u0_main, v0_main);

        --counter[e0];

        weightedDegree_counter.decrease_counter(u0_main);
        weightedDegree_counter.decrease_counter(v0_main);
        weight_counter.decrease_counter(e0_main);

        main_graph.setWeight(u0_main, v0_main, counter[e0]);
        decreaseTotalWeight();
        main_weightedDegree_sequence[u0_main] -= 1;
        main_weightedDegree_sequence[v0_main] -= 1;

        // remove edge when needed
        if (counter[e0] == 0 && main_graph.hasEdge(u0_main, v0_main)) {

            // remove edge in projection // TODO modulariser ce bout de code
            if (use_projection) {

                node u0_top = node2top[i0.u];
                node v0_bot = (is_bipartite) ? node2bot[i0.v] : node2top[i0.v];
                //if (is_bipartite) v0_bot = node2bot[i0.v];
                //removeEdgeProjection(top_graph, v0_main, u0_top, true);
                if (main_graph.degree(v0_main) < proj_bound) {
                    removeEdgeProjection(top_graph, v0_main, u0_top, true);
                } else if (main_graph.degree(v0_main) == proj_bound) {
                    forBoundedNeighbors(main_graph, v0_main, [&] (node n_main) {
                        node n = main2node[n_main];
                        node n_proj = node2top[n];
                        addEdgeProjection(v0_main, n_proj, true);
                       });
                }

                if (is_bipartite) {
                    //removeEdgeProjection(bot_graph, u0_main, v0_bot, false);
                    if (main_graph.degree(u0_main) < proj_bound) {
                        forBoundedNeighbors(main_graph, u0_main, [&] (node n_main) {
                            node n = main2node[n_main];
                            node n_proj = node2bot[n];
                            addEdgeProjection(u0_main, n_proj, false);
                           });

                    } else if (main_graph.degree(u0_main) ==  proj_bound){
                        removeEdgeProjection(bot_graph, u0_main, v0_bot, false);
                    }

                } else {
                    if (main_graph.degree(u0_main) < proj_bound) {
                        forBoundedNeighbors(main_graph, u0_main, [&] (node n_main) {
                            node n = main2node[n_main];
                            node n_proj = node2top[n];
                            addEdgeProjection(u0_main, n_proj, true);
                           });

                    } else if (main_graph.degree(u0_main) ==  proj_bound){
                        removeEdgeProjection(top_graph, u0_main, v0_bot, false);
                    }

                    //removeEdgeProjection(top_graph, u0_main, v0_bot, false);
                }

            }
            degree_counter.decrease_counter(u0_main);
            degree_counter.decrease_counter(v0_main);
            weight_counter.remove_counter(e0_main);
            main_graph.removeEdge(u0_main, v0_main);
            counter.erase(e0);

            decreaseMainDegree(u0_main, v0_main);

            if (main_graph.degree(u0_main) == 0) {
                removeNode(i0.u, true);
            }

            if (main_graph.degree(v0_main) == 0) {
                removeNode(i0.v, false);
            }
        }
    }
};
}





