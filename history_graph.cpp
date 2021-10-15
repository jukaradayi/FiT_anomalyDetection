// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include <networkit/graph/Graph.hpp>
#include <unordered_set>
#include <boost/iterator/filter_iterator.hpp>
namespace StreamGraphs {

/** constructors **/
HistoryGraph::HistoryGraph(bool use_projection, bool use_unpacked, Bound main_bound, Bound proj_bound, Count N) : 
    main_bound(main_bound), 
    proj_bound(proj_bound), 
    use_projection(use_projection), 
    use_unpacked(use_unpacked), 
    queue(0, Interaction(0,0,0)),
    removed_main(0,0),
    removed_top(0,0),
    removed_bot(0,0),
    removed_unpk(0,0),
    node2main(N,none),
    main2node(N,none),
    node2top(N,none), 
    top2node(N,none), 
    node2bot(N,none), 
    bot2node(N,none), 
    node2unpk(N,none),
    unpk2node(N,none),
    main_degree_distribution(N, 0),
    top_degree_distribution(N, 0),
    bot_degree_distribution(N, 0),
    unpk_degree_distribution(N, 0),
    main_weightedDegree_sequence(N, 0),
    top_weightedDegree_sequence(N, 0),
    bot_weightedDegree_sequence(N, 0),
    unpk_weightedDegree_sequence(N, 0),
    N(N)
      {
          std::map<Edge, uint64_t> mapping;
          this->counter = mapping;
      }

HistoryGraph::~HistoryGraph() {}

node HistoryGraph::addNode(node _u, bool is_top) {
    // TODO _addNode static function
    // update main graph
    node u;
    if (removed_main.size() == 0) {
        // new node
        u = main_graph.addNode();
    } else {
        // restore node
        u = removed_main.back();
        removed_main.pop_back();
        main_graph.restoreNode(u);
    }
    node2main[_u] = u;
    main2node[u] = _u;

    // update projection
    if (is_top && removed_top.size() > 0) {
        node u_top = removed_top.back();
        removed_top.pop_back();
        top_graph.restoreNode(u_top);
        node2top[_u] = u_top;
        top2node[u_top] = _u;
    } else if (is_top && removed_top.size() == 0) {
        node u_top = top_graph.addNode();
        node2top[_u] = u_top;
        top2node[u_top] = _u;

    } else if (!is_top && removed_bot.size() > 0) {
        node u_bot = removed_bot.back();
        removed_bot.pop_back();
        bot_graph.restoreNode(u_bot);
        node2bot[_u] = u_bot;
        bot2node[u_bot] = _u;

    } else if (!is_top && removed_bot.size() == 0) {
        node u_bot = bot_graph.addNode();
        node2bot[_u] = u_bot;
        bot2node[u_bot] = _u;
    }

    //// update unpacked
    //if (use_unpacked && removed_unpk.size() > 0) {
    //    node u_unpk = removed_unpk.back();
    //    removed_unpk.pop_back();
    //    unpk_graph.restoreNode(u_unpk);
    //    node2unpk[_u] = u_unpk;
    //    unpk2node[u_unpk] = _u;
    //    unpk_graph.restoreNode(u_unpk + 1);
    //    unpk2node[u_unpk + 1] = _u;

    //} else if (use_unpacked && removed_unpk.size() == 0) {
    //    node u_unpk = unpk_graph.addNode();
    //    node u_unpk_bis = unpk_graph.addNode();
    //    node2unpk[_u] = u_unpk;
    //    unpk2node[u_unpk] = _u;
    //    unpk2node[u_unpk + 1] = _u;

    //}
    return u;
}

void HistoryGraph::removeNode(node _u, bool is_top) {
   node u = node2main[_u];
   main_graph.removeNode(u);
   removed_main.push_back(u);
   node2main[_u] = none;
   main2node[u] = none;

   // in projection : remove node and its neighbors
   if (use_projection && is_top) {
       // remove node from top graph
       node u_top = node2top[_u];
       top_graph.removeNode(u_top);
       removed_top.push_back(u_top);

       node2top[_u] = none;
       top2node[u_top] = none;
   } else if (use_projection && !is_top) {

       // Remove node from bottom graph
       node u_bot = node2bot[_u];
       bot_graph.removeNode(u_bot);
       removed_bot.push_back(u_bot);

       node2bot[_u] = none;
       bot2node[u_bot] = none;

   }

   //if (use_unpacked) {;
   //}

}

void HistoryGraph::addEdgeProjection(node node_main, node proj_node, bool is_top) {
    auto increaseTopDegree = [&](node u, node v) {
        // to be called ~AFTER~ increasing edge weight
        //
        Count u_degree = top_graph.degree(u);
        Count v_degree = top_graph.degree(v);

        top_degree_distribution[u_degree -1] += 1;
        top_degree_distribution[v_degree -1] += 1;
        top_weightedDegree_sequence[top2node[u]] += 1;
        top_weightedDegree_sequence[top2node[v]] += 1;


        if (u_degree > 1) {
            top_degree_distribution[u_degree -2] -= 1;
        }
        if (v_degree > 1) {
            top_degree_distribution[v_degree -2] -= 1;
        }
    };
    auto increaseBotDegree = [&](node u, node v) {

        // to be called ~AFTER~ increasing edge weight
        Count u_degree = bot_graph.degree(u);
        Count v_degree = bot_graph.degree(v);

        bot_degree_distribution[u_degree -1] += 1;
        bot_degree_distribution[v_degree -1] += 1;
        bot_weightedDegree_sequence[bot2node[u]] += 1;
        bot_weightedDegree_sequence[bot2node[v]] += 1;

        if (u_degree > 1) {
            bot_degree_distribution[u_degree -2] -= 1;
        }
        if (v_degree > 1) {
            bot_degree_distribution[v_degree -2] -= 1;
        }
    };

   
    for (BoundedNeighborIterator bN = BoundedNeighborRange(main_graph, node_main, proj_bound).begin(); 
            bN != BoundedNeighborRange(main_graph, node_main, proj_bound).end(); ++bN) {
    node n = main2node[*bN];
    // for debugging purposes 
    node n_proj;
    if (is_top) {
        n_proj = node2top[n];
    } else {
        n_proj = node2bot[n];
    }

    if (n_proj == none) {
        std::cout << " n_proj is none\n";
    }
    if (n_proj == proj_node) {
        // no self loop
        continue;
    }
    if (is_top) {
        top_graph.increaseWeight(proj_node, n_proj, 1); 
        increaseTopDegree(proj_node, n_proj);
    } else {
        bot_graph.increaseWeight(proj_node, n_proj, 1); 
        increaseBotDegree(proj_node, n_proj);
    }

    }

};

void HistoryGraph::removeEdgeProjection(NetworKit::Graph& proj_graph, node node_main, node proj_node, bool is_top) {
    auto decreaseTopDegree = [&](node u, node v) {
        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = top_graph.degree(u);
        Count v_degree = top_graph.degree(v);

        if (u_degree > 0) {
            top_degree_distribution[u_degree -1] += 1;
        }
        if (v_degree > 0) {
            top_degree_distribution[v_degree -1] += 1;
        }
        top_degree_distribution[u_degree] -= 1;
        top_degree_distribution[v_degree] -= 1;

    };
    auto decreaseBotDegree = [&](node u, node v) {

        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = bot_graph.degree(u);
        Count v_degree = bot_graph.degree(v);

        if (u_degree > 0) {
            bot_degree_distribution[u_degree -1] += 1;
        }
        if (v_degree > 0) {
            bot_degree_distribution[v_degree -1] += 1;
        }
        bot_degree_distribution[u_degree] -= 1;
        bot_degree_distribution[v_degree] -= 1;

    };


        //for (BoundedNeighborIterator bN = BoundedNeighborRange(main_graph, node_main, proj_bound).begin();
        //        bN != BoundedNeighborRange(main_graph, node_main, proj_bound).end(); ++bN) {
        //node n = main2node[*bN];
        for (NetworKit::Graph::NeighborIterator N_it = main_graph.neighborRange(node_main).begin();
                N_it != main_graph.neighborRange(node_main).end(); ++N_it) {
        node n = main2node[*N_it];
        

        // for debugging purposes 
        node n_proj;
        //NetworKit::Graph* proj_graph;
        if (is_top) {
            n_proj = node2top[n];
            //proj_graph = &top_graph;
        } else {
            n_proj = node2bot[n];
            //proj_graph = &bot_graph;
        }
        if (n_proj == none) {
        }
        if (proj_graph.hasEdge(proj_node, n_proj)) {
            Count w = proj_graph.weight(proj_node, n_proj);
            if (w > 1) {
                proj_graph.setWeight(proj_node, n_proj, w-1); 
            } else {
                proj_graph.removeEdge(proj_node, n_proj); 
            }

            
        }
        if (is_top) {

            decreaseTopDegree(proj_node, n_proj);
        }else {

            decreaseBotDegree(proj_node, n_proj);
        }
        }
    };


void HistoryGraph::updateGraph(Interaction i){
    auto increaseMainDegree = [&](node u, node v) {
        // to be called ~AFTER~ increasing edge weight
        Count u_degree = main_graph.degree(u);
        Count v_degree = main_graph.degree(v);
        main_degree_distribution[u_degree -1] += 1;
        main_degree_distribution[v_degree -1] += 1;
        main_weightedDegree_sequence[u] += 1;
        main_weightedDegree_sequence[v] += 1;

        if (u_degree > 1) {
            main_degree_distribution[u_degree -2] -= 1;
        }
        if (v_degree > 1) {
            main_degree_distribution[v_degree -2] -= 1;
        }

    };
    Time t = i.t;
    Edge e(i.u, i.v);

    // if node doesn't exist, add node, else increase weight
    node u_main;
    if (node2main[i.u] == none){
        u_main = addNode(i.u, true);
    } else {
        u_main = node2main[i.u]; 
    }

    node v_main;
    if (node2main[i.v] == none){
        v_main = addNode(i.v, false);
    } else {
        v_main = node2main[i.v];
    }

    // update weight counter
    ++counter[e];

    // update queue and main graph
    queue.push_back(i);
    main_graph.increaseWeight(u_main, v_main, 1);
    increaseMainDegree(u_main, v_main); 

    // update projection graph
    if (use_projection) {
        node u_top = node2top[i.u];
        node v_bot = node2bot[i.v];

        // update top 
        
        // update bottom 

        // check if either node has passed the degree limit
        if (main_graph.degree(v_main) >= proj_bound) {
            removeEdgeProjection(top_graph, v_main, u_top, true);
        } else {
            addEdgeProjection(v_main, u_top, true);


        }


        if (main_graph.degree(u_main) >= proj_bound) {

            removeEdgeProjection(bot_graph, u_main, v_bot, false);

        } else {

            addEdgeProjection(u_main, v_bot, false);
        }

    }
    node u_top = node2top[i.u];
    node v_bot = node2bot[i.v];

    // remove links to fit window
    trimQueue(t);

}

};


