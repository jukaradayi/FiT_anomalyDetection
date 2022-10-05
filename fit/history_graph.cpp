// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

//#include <catch2/catch_test_macros.hpp> // TODO Make it optional
#include "csv.hpp"

#include "history_graph.hpp"
#include <networkit/graph/Graph.hpp>
#include <unordered_set>
//#include <boost/iterator/filter_iterator.hpp>
namespace StreamGraphs {

// TODO TODO si pas bipartite, bot = top* ? 

/** constructors **/
HistoryGraph::HistoryGraph(NetworKit::Graph& main_graph, NetworKit::Graph& top_graph, NetworKit::Graph& bot_graph, const bool use_projection, const bool use_unpacked, const bool is_bipartite, const Bound main_bound, const Bound proj_bound, Count N) : 
    main_graph(main_graph),
    top_graph(top_graph),
    bot_graph(bot_graph),
    main_bound(main_bound), 
    proj_bound(proj_bound), 
    node2main(N,none),
    main2node(N,none),
    node2top(N,none), 
    top2node(N,none), 
    node2bot(N,none), 
    bot2node(N,none), 
    node2unpk(N,none),
    unpk2node(N,none),
    use_projection(use_projection), 
    use_unpacked(use_unpacked), 
    is_bipartite(is_bipartite),
    removed_main(0,0),
    removed_top(0,0),
    removed_bot(0,0),
    removed_unpk(0,0),
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
          //this->node2main = new uint64_t[N] {none};
          //this->main2node = new uint64_t[N] {none};
          //this->node2top = new uint64_t[N] {none}; 
          //this->top2node = new uint64_t[N] {none};
          //this->node2bot = new uint64_t[N] {none}; 
          //this->bot2node = new uint64_t[N] {none};
          //this->node2unpk = new uint64_t[N] {none};
          //this->unpk2node = new uint64_t[N] {none};

          //if (is_bipartite) {
          //    this->bot_graph = top_graph;
          //} else {
          //    this->bot_graph =  bot_graph;
          //}

          std::map<Edge, uint64_t> mapping;
          this->counter = mapping;

      }

HistoryGraph::~HistoryGraph() {}

node HistoryGraph::addNode(const node _u, const bool is_top) {
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
    weightedDegree_counter.add_counter(u);
    degree_counter.add_counter(u);
    if (use_projection) {
        // update projection -- Not update bot graph when main graph is not bipartite
        if ((is_top && removed_top.size() > 0) || (!is_top && !is_bipartite && removed_top.size() >0)) {
            node u_top = removed_top.back();
            removed_top.pop_back();
            top_graph.restoreNode(u_top);
            node2top[_u] = u_top;
            top2node[u_top] = _u;
        } else if ((is_top && removed_top.size() == 0)|| (!is_top && !is_bipartite && removed_top.size() == 0)) {
            node u_top = top_graph.addNode();
            node2top[_u] = u_top;
            top2node[u_top] = _u;

        } else if (!is_top && is_bipartite && removed_bot.size() > 0) {
            node u_bot = removed_bot.back();
            removed_bot.pop_back();
            bot_graph.restoreNode(u_bot);
            node2bot[_u] = u_bot;
            bot2node[u_bot] = _u;

        } else if (!is_top && is_bipartite && removed_bot.size() == 0) {
            node u_bot = bot_graph.addNode();
            node2bot[_u] = u_bot;
            bot2node[u_bot] = _u;
        }
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

void HistoryGraph::removeNode(const node _u, const bool is_top) {
   node u = node2main[_u];
   int bbb = main_graph.numberOfNodes();
   main_graph.removeNode(u);
   removed_main.push_back(u);
   node2main[_u] = none;

   main2node[u] = none;
   degree_counter.remove_counter(u);
   weightedDegree_counter.remove_counter(u);

   // in projection : remove node and its neighbors
   if (use_projection && (is_top || !is_top && !is_bipartite)) {
       // remove node from top graph
       node u_top = node2top[_u];
       top_graph.removeNode(u_top);
       removed_top.push_back(u_top);

       node2top[_u] = none;
       top2node[u_top] = none;
   } else if (use_projection && is_bipartite && !is_top) {

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

void HistoryGraph::addEdgeProjection(const node node_main, const node proj_node,const bool is_top) {
    auto increaseTopDegree = [&](node u, node v, int64_t top_weight) {
        // to be called ~AFTER~ increasing edge weight
        //
        Count u_degree = top_graph.degree(u);
        Count v_degree = top_graph.degree(v);
        if (top_weight == 1) {
            top_degree_distribution[u_degree -1] += 1;
            top_degree_distribution[v_degree -1] += 1;
            if (u_degree > 1) {
                top_degree_distribution[u_degree -2] -= 1;
            }
            if (v_degree > 1) {
                top_degree_distribution[v_degree -2] -= 1;
            }
        }

        top_weightedDegree_sequence[top2node[u]] += 1;
        top_weightedDegree_sequence[top2node[v]] += 1;


    };
    auto increaseBotDegree = [&](node u, node v, int64_t bot_weight) {

        // to be called ~AFTER~ increasing edge weight
        Count u_degree = bot_graph.degree(u);
        Count v_degree = bot_graph.degree(v);

        if (bot_weight == 1) {
            bot_degree_distribution[u_degree -1] += 1;
            bot_degree_distribution[v_degree -1] += 1;
            if (u_degree > 1) {
                bot_degree_distribution[u_degree -2] -= 1;
            }
            if (v_degree > 1) {
                bot_degree_distribution[v_degree -2] -= 1;
            }
        }
        bot_weightedDegree_sequence[bot2node[u]] += 1;
        bot_weightedDegree_sequence[bot2node[v]] += 1;

    };

    // iterate over bounded neighbors to add nodes 
    //for (BoundedNeighborIterator bN = BoundedNeighborRange(main_graph, node_main, proj_bound).begin(); 
    //        bN != BoundedNeighborRange(main_graph, node_main, proj_bound).end(); ++bN) {
    //for (NetworKit::Graph::NeighborIterator bN = main_graph.neighborRange(node_main).begin();
    //        bN != main_graph.neighborRange(node_main).end(); ++bN) {
    forBoundedNeighbors(main_graph, node_main, [&] (node n_main) {


        node n = main2node[n_main];
        ////node n = main2node[*N_it];
        //if (main_graph.degree(*bN) >= main_bound) {
        //     continue;
        //}

        node n_proj;
        if (is_top || (!is_top && !is_bipartite)) {
            n_proj = node2top[n];
        } else if (!is_top && is_bipartite) {
            n_proj = node2bot[n];
        }

        if (n_proj == none) { // should never happen
            std::cout << " n_proj is none\n";
        }
        if (n_proj == proj_node) {
            // no self loop
            return;
        } // TODO CHANGE IS TOP 
        if (is_top || (!is_top && !is_bipartite)) {
            // redundancy_top = (top_proj.hasEdge(proj_node, n_proj)) ? redundancy_top : redundancy_top + 1;
            top_graph.increaseWeight(proj_node, n_proj, 1); 
            Count w = top_graph.weight(proj_node, n_proj);
            increaseTopDegree(proj_node, n_proj, w);
        } else if (!is_top && is_bipartite){
            bot_graph.increaseWeight(proj_node, n_proj, 1); 
            Count w = bot_graph.weight(proj_node, n_proj);
            increaseBotDegree(proj_node, n_proj, w);
        }

    });

};

void HistoryGraph::removeEdgeProjection(NetworKit::Graph& proj_graph,const node node_main, const node proj_node, const bool is_top) {
    auto decreaseTopDegree = [&](node u, node v, uint64_t weight) {
        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = top_graph.degree(u);
        Count v_degree = top_graph.degree(v);

        if (weight == 1) {
            if (u_degree > 0) {
                top_degree_distribution[u_degree -1] += 1;
            }
            if (v_degree > 0) {
                top_degree_distribution[v_degree -1] += 1;
            }
            top_degree_distribution[u_degree] -= 1;
            top_degree_distribution[v_degree] -= 1;
        }
        top_weightedDegree_sequence[u] -= 1;
        top_weightedDegree_sequence[v] -= 1;


    };
    auto decreaseBotDegree = [&](node u, node v, uint64_t weight) {

        // to be called ~AFTER~ decreasing edge weight
        Count u_degree = bot_graph.degree(u);
        Count v_degree = bot_graph.degree(v);
        
        if (weight == 1) {
            if (u_degree > 0) {
                bot_degree_distribution[u_degree -1] += 1;
            }
            if (v_degree > 0) {
                bot_degree_distribution[v_degree -1] += 1;
            }
            bot_degree_distribution[u_degree] -= 1;
            bot_degree_distribution[v_degree] -= 1;
        }
        bot_weightedDegree_sequence[u] -= 1;
        bot_weightedDegree_sequence[v] -= 1;

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
            if (is_top || (!is_top && !is_bipartite)) {
                n_proj = node2top[n];
                //proj_graph = &top_graph;
            } else  if (!is_top && is_bipartite){
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
                // decrease degree
                if (is_top || (!is_top && !is_bipartite)) {

                    decreaseTopDegree(proj_node, n_proj, w);
                }else if (!is_top && is_bipartite){

                    decreaseBotDegree(proj_node, n_proj, w);
                }


                
            }
        }
    };


void HistoryGraph::updateGraph(const Interaction i){
    auto increaseMainDegree = [&](node u, node v, int64_t weight) {
        // to be called ~AFTER~ increasing edge weight
        // increase degree only if weight is equal to 1 (i.e. if edges was not present before)
        Count u_degree = main_graph.degree(u);
        Count v_degree = main_graph.degree(v);
        if (weight == 1) {
            main_degree_distribution[u_degree -1] += 1;
            main_degree_distribution[v_degree -1] += 1;
            if (u_degree > 1) {
                main_degree_distribution[u_degree -2] -= 1;
            }
            if (v_degree > 1) {
                main_degree_distribution[v_degree -2] -= 1;
            }
        }
        main_weightedDegree_sequence[u] += 1;
        main_weightedDegree_sequence[v] += 1;

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
    queue.push(i);

    // update counters
    if (!main_graph.hasEdge(u_main, v_main)) {
        degree_counter.increase_counter(u_main);
        degree_counter.increase_counter(v_main);
    }
    weightedDegree_counter.increase_counter(u_main);
    weightedDegree_counter.increase_counter(v_main);
   
    main_graph.increaseWeight(u_main, v_main, 1);
    //main_graph.max_weighted_degree(u_main, v_main);

    increaseTotalWeight();
    increaseMainDegree(u_main, v_main, counter[e]);
    
    // update projection graph
    if (use_projection) {
        node u_top = node2top[i.u];
        node v_bot = (is_bipartite) ? node2bot[i.v] : node2top[i.v];
        //if (is_bipartite) v_bot = node2bot[i.v];

        // update top 
        // check if either node has passed the degree limit
        if (main_graph.degree(v_main) >= proj_bound) {
            //removeEdgeProjection(top_graph, v_main, u_top, true);
            forBoundedNeighbors(main_graph, v_main, [&] (node n_main) {
                node n = main2node[n_main];
                node n_proj = node2top[n];
                removeEdgeProjection(top_graph, v_main, n_proj, true);
            });

        } else {
            addEdgeProjection(v_main, u_top, true);
        }

        if (main_graph.degree(u_main) >= proj_bound && is_bipartite) {
            //removeEdgeProjection(bot_graph, u_main, v_bot, false);
            forBoundedNeighbors(main_graph, u_main, [&] (node n_main) {
                node n = main2node[n_main];
                node n_proj = node2bot[n];
                removeEdgeProjection(bot_graph, u_main, n_proj, true);
            });

        } else if (is_bipartite) {
            addEdgeProjection(u_main, v_bot, false);
        }

        // update bottom 
        if (main_graph.degree(u_main) >= proj_bound && !is_bipartite) {
            //removeEdgeProjection(top_graph, u_main, v_bot, false);
            forBoundedNeighbors(main_graph, u_main, [&] (node n_main) {
                node n = main2node[n_main];
                node n_proj = node2top[n];
                removeEdgeProjection(top_graph, u_main, n_proj, true);
            });

        } else if (!is_bipartite) {
            addEdgeProjection(u_main, v_bot, false);
        }

    }
    // remove links to fit window
    trimQueue(t);
}

};

