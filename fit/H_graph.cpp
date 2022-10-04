// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "H_graph.hpp"
#include <networkit/graph/Graph.hpp>
#include <iostream>

namespace StreamGraphs {
using namespace StreamGraphs;
/** constructors **/
HGraph::HGraph(NetworKit::Graph& main_graph, NetworKit::Graph& top_graph, NetworKit::Graph& bot_graph, const bool use_projection, const bool use_unpacked, const bool is_bipartite, const Bound main_bound, const Bound proj_bound, Count N, const Bound window) : HistoryGraph(main_graph, top_graph, bot_graph, use_projection, use_unpacked, is_bipartite, main_bound, proj_bound, N), window(window) {}

HGraph::~HGraph() {}

void HGraph::trimQueue(Time t) {
    auto decreaseMainDegree = [&](node u, node v) {
        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = main_graph.degree(u);
        Count v_degree = main_graph.degree(v);
        main_degree_distribution[u_degree] -= 1;
        main_degree_distribution[v_degree] -= 1;
        if (u_degree > 0) {
            main_degree_distribution[u_degree -1] += 1;
        }
        if (v_degree > 0) {
            main_degree_distribution[v_degree -1] += 1;
        }
    };
    while (queue.size() > window) { // could be an if

        Interaction i0 = queue.front();

        queue.pop();

        node u0_main = node2main[i0.u];
        node v0_main = node2main[i0.v];
        
        // decrease edge counter
        Edge e0(i0.u, i0.v);
        --counter[e0];

        main_graph.setWeight(u0_main, v0_main, counter[e0]);
        decreaseTotalWeight();
        main_weightedDegree_sequence[u0_main] -= 1;
        main_weightedDegree_sequence[v0_main] -= 1;

        // remove edge when needed
        if (counter[e0] == 0 && main_graph.hasEdge(u0_main, v0_main)) { // TODO virer check hasEdge

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
                    //removeEdgeProjection(top_graph, u0_main, v0_bot, false);
                    if (main_graph.degree(u0_main) < proj_bound) {
                        forBoundedNeighbors(main_graph, u0_main, [&] (node n_main) {
                            node n = main2node[n_main];
                            node n_proj = node2top[n];
                            addEdgeProjection(u0_main, n_proj, true);
                           });

                    } else if (main_graph.degree(u0_main) ==  proj_bound){
                        removeEdgeProjection(top_graph, u0_main, v0_bot, false);
                    }


                }
            }
            std::cout << " decreasing \n";
            degree_counter.decrease_counter(u0_main);
            degree_counter.decrease_counter(v0_main);
            main_graph.removeEdge(u0_main, v0_main);
            //std::cout << "node u_main "<<u0_main<<" main " << main_graph.degree(u0_main) << " counter " << degree_counter.get_value(u0_main) <<" \n";

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





