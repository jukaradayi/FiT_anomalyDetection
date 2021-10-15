// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "H_graph.hpp"
#include <networkit/graph/Graph.hpp>

namespace StreamGraphs {
using namespace StreamGraphs;
/** constructors **/
HGraph::HGraph(bool use_projection, bool use_unpacked, Bound main_bound, Bound proj_bound, Count N, Bound window) : HistoryGraph( use_projection, use_unpacked, main_bound, proj_bound, N), window(window) {}

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
    while (queue.size() > window) {

        Interaction i0 = queue.front();
        queue.erase(queue.begin()); 

        node u0_main = node2main[i0.u];
        node v0_main = node2main[i0.v];
        
        // decrease edge counter
        Edge e0(i0.u, i0.v);
        --counter[e0];

        main_graph.setWeight(u0_main, v0_main, counter[e0]);
        main_weightedDegree_sequence[i0.u] -= 1;
        main_weightedDegree_sequence[i0.v] -= 1;

        // remove edge when needed
        if (counter[e0] == 0 && main_graph.hasEdge(u0_main, v0_main)) {

            // remove edge in projection // TODO modulariser ce bout de code
            if (use_projection) {
                node u0_top = node2top[i0.u];
                node v0_bot = node2bot[i0.v];

                removeEdgeProjection(top_graph, v0_main, u0_top, true);
                removeEdgeProjection(bot_graph, u0_main, v0_bot, false);
            }
            main_graph.removeEdge(u0_main, v0_main);
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





