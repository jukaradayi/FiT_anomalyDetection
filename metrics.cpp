// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree
// TODO slow:  egonets, BFS

#include "G_graph.hpp"
#include "history_graph.hpp"
#include "metrics.hpp"
#include "PLM.hpp"
#include <omp.h>
#include <chrono>
#include <ctime>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ClusteringProjector.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <iostream>
#include <fstream>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/centrality/LocalClusteringCoefficient.hpp>
#include <networkit/linkprediction/JaccardIndex.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/centrality/PageRank.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/distance/BFS.hpp>

namespace StreamGraphs {
//using namespace StreamGraphs;
using namespace NetworKit;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


Metrics::Metrics(HistoryGraph& history,const bool use_basic,const bool use_local,const bool use_global,const bool use_nonLinear): history(history), use_basic(use_basic), use_local(use_local), use_global(use_global), use_nonLinear(use_nonLinear), header_done(false) {}

Metrics::~Metrics() {}

//double Metrics::localClustering(const node u, const node v) {
//    return 1.0;
//}

std::pair<double, double> Metrics::localClustering(const NetworKit::Graph& G, const node u, const node v) {
    count z = G.upperNodeIdBound();
    std::vector<double> scoreData;
    scoreData.clear();
    scoreData.resize(z); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$
    std::vector<index> inBegin;
    std::vector<node> inEdges;

    //if (turbo) {
    auto isOutEdge = [&](node u, node v) {
        return G.degree(u) > G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
    };

    inBegin.resize(G.upperNodeIdBound() + 1);
    inEdges.resize(G.numberOfEdges());
    index pos = 0;

    for (index u = 0; u < G.upperNodeIdBound(); ++u) {
        inBegin[u] = pos;
        if (G.hasNode(u)) {
            G.forEdgesOf(u, [&](node v) {
                if (isOutEdge(v, u)) {
                    inEdges[pos++] = v;
                }
            });
        }
    }

    inBegin[G.upperNodeIdBound()] = pos;
    //}

    std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());

    for (auto & nm : nodeMarker) {
        nm.resize(z, false);
    }


    G.balancedParallelForNodes([&](node u_node) {
    //auto computeClustering = [&](node u) {
    //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_main, history.main_bound).begin();
    //    bN != BoundedNeighborRange(history.main_graph, u_main, history.main_bound).end(); ++bN) {
    //    node n = *bN;
        if (u_node != u && u_node != v) {
            return;
        }
        count d = G.degree(u_node);

        if (d < 2) {

            scoreData[u_node] = 0.0;

        } else {

            size_t tid = omp_get_thread_num();
            count triangles = 0;

            //G.forEdgesOf(n, [&](node v) {

        //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_node, history.main_bound).begin();
        //    bN != BoundedNeighborRange(history.main_graph, u_node, history.main_bound).end(); ++bN) {
        for (NetworKit::Graph::NeighborIterator bN = history.main_graph.neighborRange(u_node).begin();
             bN != history.main_graph.neighborRange(u_node).end(); ++bN) {

                //node nn_main = *bNN;
                node n = *bN;
                if (history.main_graph.degree(n) >= history.main_bound) {
                    continue;
                }

                //node n = *bN;

                nodeMarker[tid][n] = true;
            }

            //G.forEdgesOf(u, [&](node, node v) {

        //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_node, history.main_bound).begin();
        //    bN != BoundedNeighborRange(history.main_graph, u_node, history.main_bound).end(); ++bN) {
        //        node n = *bN;
         for (NetworKit::Graph::NeighborIterator bN = history.main_graph.neighborRange(u_node).begin();
             bN != history.main_graph.neighborRange(u_node).end(); ++bN) {

                //node nn_main = *bNN;
                node n = *bN;
                if (history.main_graph.degree(n) >= history.main_bound) {
                    continue;
                }

                //if (turbo) {
                for (index i = inBegin[n]; i < inBegin[n+1]; ++i) {

                    node w = inEdges[i];

                    if (nodeMarker[tid][w]) {
                        triangles += 1;
                    }
                }

                //} else {
                //    G.forEdgesOf(v, [&](node, node w) {
                //        if (nodeMarker[tid][w]) {
                //            triangles += 1;
                //        }
                //    });
                //}
            }

            //G.forEdgesOf(u, [&](node, node v) {
            //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_node, history.main_bound).begin();
            //    bN != BoundedNeighborRange(history.main_graph, u_node, history.main_bound).end(); ++bN) {
            for (NetworKit::Graph::NeighborIterator bN = history.main_graph.neighborRange(u_node).begin();
                 bN != history.main_graph.neighborRange(u_node).end(); ++bN) {

                    //node nn_main = *bNN;
                    node n = *bN;
                    if (history.main_graph.degree(n) >= history.main_bound) {
                        continue;
                    }

                    //node n = *bN;
                    nodeMarker[tid][n] = false;
            }

            scoreData[u_node] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
            //if (turbo) 
            scoreData[u_node] *= 2; // in turbo mode, we count each triangle only once

        }

    });
    //computeClustering(u);
    //std::cout << "compute u" << scoreData[u] << " \n";


    //computeClustering(v);
    //std::cout << "compute v " << scoreData[v] << " \n";

    //std::cout << "make pair \n";

    std::pair<double, double> res = std::make_pair(scoreData[u], scoreData[v]);
    //hasRun = true;

    return res;


}

std::string Metrics::run(const node u,const node v, const int interaction_id) {
    // initialisations
    Edge e(u, v);

    node u_main = history.node2main[u];
    node v_main = history.node2main[v];
    
    node u_top;
    node v_bot;
    if (history.use_projection) {
        u_top = history.node2top[u];
        v_bot = (history.is_bipartite) ? history.node2bot[v] : history.node2top[v];
    }
    int N_plm = 5;

    //if (history.is_bipartite) v_bot = history.node2bot[v];
    integerResults["interaction_id"] = interaction_id;

    //if (use_basic) {
    //    /***************/
    //    /*Basic metrics*/
    //    /***************/
    if (use_basic) {

        integerResults["number of nodes"] = history.main_graph.numberOfNodes(); 
        integerResults["number of links"] = history.main_graph.numberOfEdges();
        integerResults["count nodes degree 1"] = history.main_degree_distribution[0];
        integerResults["count nodes degree 2"] = history.main_degree_distribution[1];

        integerResults["degree u"] = degree(u);
        integerResults["degree v"] = degree(v);

        // TODO En O(degré(u))
        integerResults["weighted degree u"] = weighted_degree(u);
        integerResults["weighted degree v"] = weighted_degree(v);

        // TODO en O(n) - devrait se changer facilement
        //integerResults["max degree"] = max_degree();

        //integerResults["max weighted degree"] = max_weighted_degree();

        integerResults["degree absolute difference"] = degree_absolute_difference(u, v);

        integerResults["weighted degree absolute difference"] = weighted_degree_absolute_difference(u, v);

        integerResults["weight"] = history.counter[e];

        integerResults["total weight"] = total_weight();

        /***********************/
        /*Basic derived metrics*/
        /***********************/

        doubleResults["average weight"] = integerResults["total weight"]/integerResults["number of links"];

        doubleResults["average degree"] = 2 * integerResults["number of links"] / integerResults["number of nodes"];

        doubleResults["average weighted degree"] = integerResults["total weight"]/(2*integerResults["number of nodes"]);

        doubleResults["density"] = 2 * integerResults["number of links"]/(integerResults["number of nodes"] *(integerResults["number of nodes"] - 1)) ;



        /********************/
        /*Projection metrics*/
        /********************/
        if (history.use_projection) {
            //TODO keep top & bot nodes directly from main graph
            integerResults["number of top nodes"] = history.top_graph.numberOfNodes();
            if (history.is_bipartite) integerResults["number of bot nodes"] = history.bot_graph.numberOfNodes();
            integerResults["number of top links"] = history.top_graph.numberOfEdges();
            if (history.is_bipartite) integerResults["number of bot links"] = history.bot_graph.numberOfEdges();
            integerResults["top degree"] = top_degree(u);
            if (history.is_bipartite) {
                integerResults["bot degree"] = bot_degree(v);
            } else {
                integerResults["top degree v"] = top_degree(v);
            }

            integerResults["top weighted degree"] = top_weighted_degree(u);
            if (history.is_bipartite) {
                integerResults["bot weighted degree"] = bot_weighted_degree(v);
            } else {
                integerResults["top weighted degree"] = top_weighted_degree(v);
            }

            //integerResults["top max weighted degree"] = top_max_weighted_degree();
            //if (history.is_bipartite) integerResults["bot max weighted degree"] = bot_max_weighted_degree();
        }
    };
    
    if (use_local) {
        //    //// sets
        //    // get distance 1 and 2 neighbor sets, basic sets

        std::unordered_set<node> bNNu;
        std::unordered_set<node> bNNv;
        std::unordered_set<node> bNu;
        std::unordered_set<node> bNv;
        std::unordered_set<node> bNu_n_bNNv;
        std::unordered_set<node> bNNu_n_bNv;
        std::unordered_set<node> bNNu_u_bNv;
        std::unordered_set<node> bNu_u_bNNv;

        std::unordered_set<node> bNNuNv_u_bNuNNv;
        std::unordered_set<node> bNNuNv_n_bNuNNv;
        std::unordered_set<node> bNu_u_bNNu;
        std::unordered_set<node> bNv_u_bNNv;
        std::unordered_set<node> bNu_u_bNv;
        std::unordered_set<node> bNu_n_bNv;

        //// fill sets
        //std::cout << "u " << u_main << " v " << v_main << "\n";
        //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_main, history.main_bound).begin();
        //    bN != BoundedNeighborRange(history.main_graph, u_main, history.main_bound).end(); ++bN) {
        for (NetworKit::Graph::NeighborIterator bN = history.main_graph.neighborRange(u_main).begin();
                bN != history.main_graph.neighborRange(u_main).end(); ++bN) { // boundedIncrementor(history.main_graph, history.main_bound, bN, history.main_graph.neighborRange(u_main).end())) {
            node n_main = *bN;
            //node n_main = *N_it;

            if (history.main_graph.degree(n_main) >= history.main_bound) {
                continue;
            }
            //std::cout << n_main << " is a neighbor of u with degree bipbip " << history.main_graph.degree(n_main) << " \n";

            bNu.insert(n_main);
            bNu_u_bNNv.insert(n_main);

            //for (BoundedNeighborIterator bNN = BoundedNeighborRange(history.main_graph, n_main, history.main_bound).begin();
            //    bNN != BoundedNeighborRange(history.main_graph, n_main, history.main_bound).end(); ++bNN) {
            for (NetworKit::Graph::NeighborIterator bNN = history.main_graph.neighborRange(n_main).begin();
                 bNN != history.main_graph.neighborRange(n_main).end(); ++bNN) {

                node nn_main = *bNN;
                //node nn_main = *NN_it;
                if (history.main_graph.degree(nn_main) >= history.main_bound) {
                    continue;
                }

                bNNu.insert(nn_main);
                bNNu_u_bNv.insert(nn_main);
            }
        }

        //for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, v_main, history.main_bound).begin();
        //    bN != BoundedNeighborRange(history.main_graph, v_main, history.main_bound).end(); ++bN) {
        for (NetworKit::Graph::NeighborIterator bN = history.main_graph.neighborRange(v_main).begin();
                bN != history.main_graph.neighborRange(v_main).end(); ++bN) {

            node n_main = *bN;


        //for (NetworKit::Graph::NeighborIterator N_it = history.main_graph.neighborRange(v_main).begin();
        //        N_it != history.main_graph.neighborRange(v_main).end(); ++N_it) {
            //node n_main = *bN;
            //node n_main = *N_it;

            if (history.main_graph.degree(n_main) >= history.main_bound) {
                continue;
            }
            //std::cout << n_main << " is a neighbor of v with degree " << history.main_graph.degree(n_main) << " \n";

            bNv.insert(n_main);
            //std::pair<std::unordered_set<node>,bool> insertion = bNNu_u_bNv.insert(n_main);
            auto insertion = bNNu_u_bNv.insert(n_main);


            //if (bNNu.find(n_main) != bNNu.end()) {
            // if insertion.second returns false, it means the element already exists in bNNu_u_bNv 
            if (insertion.second == false){
                bNNu_n_bNv.insert(n_main);
            }
            //for (BoundedNeighborIterator bNN = BoundedNeighborRange(history.main_graph, n_main, history.main_bound).begin();
            //    bNN != BoundedNeighborRange(history.main_graph, n_main, history.main_bound).end(); ++bNN) {
            //    node nn_main = *bNN;
            for (NetworKit::Graph::NeighborIterator bNN = history.main_graph.neighborRange(n_main).begin();
                 bNN != history.main_graph.neighborRange(n_main).end(); ++bNN) {

                //node nn_main = *bNN;
                node nn_main = *bNN;
                if (history.main_graph.degree(nn_main) >= history.main_bound) {
                    continue;
                }

                bNNv.insert(nn_main);
                //std::pair<std::unordered_set<node>,bool> insertion = bNu_u_bNNv.insert(nn_main);
                auto insertion = bNu_u_bNNv.insert(nn_main);


                //if (bNu.find(nn_main) != bNu.end()) {
                if (insertion.second == false){
                    bNu_n_bNNv.insert(nn_main);
                }
            }
        }

        //    // build egonets

        //    // egonet sets

        us_isect(bNNuNv_n_bNuNNv, bNNu_u_bNv, bNu_u_bNNv); // intersection builder
        us_isect(bNu_n_bNv, bNu, bNv); // intersection builder

        bNNuNv_u_bNuNNv.insert(bNNu_u_bNv.begin(),bNNu_u_bNv.end());
        bNNuNv_u_bNuNNv.insert(bNu_u_bNNv.begin(),bNu_u_bNNv.end());
        bNu_u_bNNu.insert(bNu.begin(), bNu.end());
        bNu_u_bNNu.insert(bNNu.begin(), bNNu.end());
        bNv_u_bNNv.insert(bNv.begin(), bNv.end());
        bNv_u_bNNv.insert(bNNv.begin(), bNNv.end());
        bNu_u_bNv.insert(bNu.begin(), bNu.end());
        bNu_u_bNv.insert(bNv.begin(), bNv.end());

        // egonets from node sets 
        Graph egonet_NNuNv_u_NuNNv = GraphTools::subgraphFromNodes(history.main_graph, bNNuNv_u_bNuNNv);
        Graph egonet_NNuNv_n_NuNNv = GraphTools::subgraphFromNodes(egonet_NNuNv_u_NuNNv, bNNuNv_u_bNuNNv);
        Graph egonet_Nu_u_NNu = GraphTools::subgraphFromNodes(egonet_NNuNv_u_NuNNv, bNu_u_bNNu);
        Graph egonet_Nv_u_NNv = GraphTools::subgraphFromNodes(egonet_NNuNv_u_NuNNv, bNv_u_bNNv);
        Graph egonet_Nu_u_Nv = GraphTools::subgraphFromNodes(egonet_NNuNv_u_NuNNv, bNu_u_bNv);
        Graph egonet_Nu = GraphTools::subgraphFromNodes(egonet_Nu_u_Nv, bNu);
        Graph egonet_Nv = GraphTools::subgraphFromNodes(egonet_Nu_u_Nv, bNv);

        //    /******************/
        //    /*local clustering*/
        //    /******************/

        //LocalClusteringCoefficient clustering(history.main_graph, true);
        std::pair<double, double> clustering = localClustering(history.main_graph, u_main, v_main);
        doubleResults["clustering u"] = clustering.first;

        doubleResults["clustering v"] = clustering.second;
        //clustering.run();

        //doubleResults["clustering u"] = clustering.score(u_main);
        //doubleResults["clustering v"] = clustering.score(v_main);




        /*********************************/
        /*Jaccard index and neighborhoods*/
        /*********************************/
        //JaccardIndex jaccard(history.main_graph);
        if (bNu_u_bNv.size() > 0) {
           doubleResults["Jaccard index u v"] = std::abs(static_cast<int>(bNu_n_bNv.size() / bNu_u_bNv.size()));
           
        }
        // TODO replace Jaccard by actual computation
        // TODO Jaccard (& Jaccard bipartite on actual neighborhoods or bounded neighborhoods ?
        if (bNu_u_bNNv.size() > 0) {
            doubleResults["Jaccard index bipartite u"] = std::abs(static_cast<int>(bNu_n_bNNv.size() / bNu_u_bNNv.size()));
        } else {
            doubleResults["Jaccard index bipartite u"] = 0;
        }

        if (bNNu_u_bNv.size() > 0) {
            doubleResults["Jaccard index bipartite v"] = std::abs(static_cast<int>(bNNu_n_bNv.size() / bNNu_u_bNv.size()));
        } else {
            doubleResults["Jaccard index bipartite v"] = 0; 
        }
        if (bNu.size() > 0) {
            doubleResults["ratio neighbors of u in N(v)"] = bNu_n_bNv.size() / bNu.size();
        } else {
            doubleResults["ratio neighbors of u in N(v)"] = 0;
        }
        if (bNv.size() > 0) {
            doubleResults["ratio neighbors of v in N(u)"] = bNu_n_bNv.size() / bNv.size();
        } else {
            doubleResults["ratio neighbors of v in N(u)"] = 0; 
        }
        doubleResults["number common neighbors expected from degree"] = bNu_n_bNv.size() / std::sqrt(degree(u) * degree(v));

        integerResults["neighborhood overlap"] = bNu_u_bNv.size();
        //integerResults["top nehborhood size"] =  
        //integerResults["top nehborhood size"]

        integerResults["egonet Nu number of links"] = egonet_Nu.numberOfEdges();
        integerResults["egonet Nu number of nodes"] = egonet_Nu.numberOfNodes();
        integerResults["egonet Nv number of links"]  = egonet_Nv.numberOfEdges();
        integerResults["egonet Nv number of nodes"]  = egonet_Nv.numberOfNodes();

        integerResults["egonet Nu u Nv number of links"] = egonet_Nu_u_Nv.numberOfEdges();
        integerResults["egonet Nu u Nv number of nodes"]  = egonet_Nu_u_Nv.numberOfNodes();

        integerResults["egonet Nv u NNv number of links"] = egonet_Nv_u_NNv.numberOfEdges();
        integerResults["egonet Nv u NNv number of nodes"]  = egonet_Nv_u_NNv.numberOfNodes();
        integerResults["egonet Nu u NNu number of links"] = egonet_Nu_u_NNu.numberOfEdges();
        integerResults["egonet Nu u NNu number of nodes"]  = egonet_Nu_u_NNu.numberOfNodes();
        integerResults["egonet (NNu u Nv) u (Nu u NNv) number of links"] = egonet_NNuNv_u_NuNNv.numberOfEdges();
        integerResults["egonet (NNu u Nv) u (Nu u NNv) number of nodes"] = egonet_NNuNv_u_NuNNv.numberOfNodes();
        integerResults["egonet (NNu u Nv) n (Nu u NNv) number of links"] = egonet_NNuNv_n_NuNNv.numberOfEdges();
        integerResults["egonet (NNu u Nv) n (Nu u NNv) number of nodes"] = egonet_NNuNv_n_NuNNv.numberOfNodes();

        integerResults["egonet Nu u Nv maxsize"] = bNu.size() * bNv.size();
        integerResults["egonet Nv u NNv maxsize"] = bNv.size() * bNNv.size();
        integerResults["egonet Nu u NNu maxsize"] = bNu.size() * bNNu.size();
        integerResults["egonet (NNu u Nv) u (Nu u NNv) maxsize"] = bNu_u_bNNv.size() * bNNu_u_bNv.size();
        integerResults["egonet (NNu u Nv) n (Nu u NNv) maxsize"] =  bNu_n_bNNv.size() * bNNu_n_bNv.size();

        // adamic adar
        doubleResults["adamic adar"] = adamic_adar(bNu_n_bNv);
        doubleResults["adamic adar bipartite u"] = adamic_adar_bipartite_u(bNu_n_bNNv);
        doubleResults["adamic adar bipartite v"] = adamic_adar_bipartite_v(bNNu_n_bNv);

        // local metrics on projection
        if (history.use_projection) {
            LocalClusteringCoefficient top_clustering(history.top_graph, true);
            top_clustering.run();

            doubleResults["top clustering u"] = top_clustering.score(u_top);

            if (history.is_bipartite){

                LocalClusteringCoefficient bot_clustering(history.bot_graph, true);
                bot_clustering.run();
                doubleResults["bot clustering v"] = bot_clustering.score(v_bot);
            } else {
                doubleResults["top clustering v"] = top_clustering.score(v_bot);

            }
        }

    }
    if (use_global) { 
        /************/
        /*components*/
        /************/

        ConnectedComponents components_with(history.main_graph);
        components_with.run();
        std::map<index, Count> components_with_sizes = components_with.getComponentSizes();
        Count max_size_component_with = 0;
        for (std::map<index, Count>::iterator it = components_with_sizes.begin(); it != components_with_sizes.end(); ++it) {
            Count component_size = it->second;
            if (component_size > max_size_component_with) {
                max_size_component_with = component_size;
            }
        }
        std::vector<std::vector<node>> all_components = components_with.getComponents();
        std::unordered_set<node> component_u_nodes_with(all_components[components_with.componentOfNode(u_main)].begin(), 
                                                        all_components[components_with.componentOfNode(u_main)].end());

        Graph link_component = GraphTools::subgraphFromNodes(history.main_graph, component_u_nodes_with);

        integerResults["link component number of nodes"] = link_component.numberOfNodes();
        integerResults["link component number of links"] = link_component.numberOfEdges();
        integerResults["number of components G"] = components_with.numberOfComponents();
        integerResults["size largest component G"] = max_size_component_with;

        /*************/
        /*BFS metrics*/
        /*************/
        // parallelize BFS ? 
        //std::vector<BFS> BFS_with;
        //BFS_with.push_back(BFS(history.main_graph, u_main, true));
        //BFS_with.push_back(BFS(history.main_graph, v_main, true));
        //#pragma omp parallel for num_threads(2)
        //for (omp_index u = 0; u < static_cast<omp_index>(2); ++u){
        //    BFS_with[u].run();
        //}
        //BFS BFSu_with = BFS_with.back();
        //BFS_with.pop_back();
        //BFS BFSv_with = BFS_with.back();
        //BFS_with.pop_back();

        BFS BFSu_with = BFS(history.main_graph, u_main, true);
        BFS BFSv_with = BFS(history.main_graph, v_main, true);
        BFSu_with.run();
        BFSv_with.run();
        
        /********************/
        /*Core Decomposition*/
        /********************/

        CoreDecomposition core_with = CoreDecomposition(history.main_graph, false);
        core_with.run();

        //if (use_nonLinear) {
        /**********/
        /*PageRank*/
        /**********/
        PageRank pagerank_with(history.main_graph);
        pagerank_with.run();
        //}

        /************************/
        /************************/
        /* Compute metrics on G-*/
        /************************/
        /************************/
        // Remove (u,v)
        history.main_graph.removeEdge(u_main, v_main);

        //if (use_global) {
        //std::vector<BFS> BFS_without;
        //BFS_without.push_back(BFS(history.main_graph, u_main, true));
        //BFS_without.push_back(BFS(history.main_graph, v_main, true));

        //#pragma omp parallel for num_threads(1)
        //for (omp_index u = 0; u < static_cast<omp_index>(2); ++u){
        //    BFS_without[u].run();
        //}

        //BFS BFSu_with = all_BFS.back(); 
        //all_BFS.pop_back();
        //BFS BFSv_with = all_BFS.back();
        //all_BFS.pop_back();

        BFS BFSu_without = BFS(history.main_graph, u_main, true);
        BFS BFSv_without = BFS(history.main_graph, v_main, true);
        BFSu_without.run();
        BFSv_without.run();

        /********************/
        /*Core Decomposition*/
        /********************/
        CoreDecomposition core_without = CoreDecomposition(history.main_graph, false);
        core_without.run();

        /**********/
        /*PageRank*/
        /**********/
        PageRank pagerank_without(history.main_graph);
        pagerank_without.run();


        /**********************/
        /*Connected Components*/
        /**********************/
        // TODO run or just deduce from BFS ? ..
        ConnectedComponents components_without(history.main_graph);
        components_without.run();

        //# compute distance changes from u and v
        Count N = history.main_graph.numberOfNodes(); 

        Count dist_change_u = 0;
        Count dist_change_v = 0;
        Count max_dist_change_u = 0;
        Count max_dist_change_v = 0;
        Count N_dist_change_u = 0;
        Count N_dist_change_v = 0;
        Count sum_diff_dist_u_v_with = 0;
        Count sum_diff_dist_u_v_without = 0;

        // keep list of distances to accessible nodes
        //std::vector<Count> dist_to_u_with(N, 0);
        //std::vector<Count> dist_to_v_with(N, 0);
        //std::vector<Count> dist_to_u_without(N, 0);
        //std::vector<Count> dist_to_v_without(N, 0);
        Count* dist_to_u_with = new Count[N]{0};
        Count* dist_to_v_with = new Count[N]{0};
        Count* dist_to_u_without = new Count[N]{0};
        Count* dist_to_v_without = new Count[N]{0};


        double pagerank_variation = 0;
        double pagerank_max_with = 0;
        double pagerank_max_without = 0;
        double core_variation = 0;
        Count idx = 0;
        for (const auto n : history.main_graph.nodeRange()){

            // Loop over all nodes to compute variations in several metrics
            if (BFSu_without.numberOfPaths(n) > 0){
                Count dist_u_with = BFSu_with.distance(n);
                Count dist_u_without = BFSu_without.distance(n);

                //dist_to_u_with.push_back(dist_u_with);
                //dist_to_u_without.push_back(dist_u_without);
                dist_to_u_with[idx] = dist_u_with;
                dist_to_u_without[idx] = dist_u_without;
                Count dist_change = std::abs(static_cast<int>(dist_u_with - dist_u_without));
                dist_change_u += dist_change; 
                max_dist_change_u = (max_dist_change_u < dist_change) ? dist_change : max_dist_change_u;
                N_dist_change_u = (dist_change > 0) ? N_dist_change_u + 1 :N_dist_change_u;

                // if there's a path to u in G, there's one to v ..
                sum_diff_dist_u_v_with += std::abs(static_cast<int>(dist_u_with - BFSv_with.distance(n)));

            }

            if (BFSv_without.numberOfPaths(n) > 0){
                Count dist_v_with = BFSv_with.distance(n);
                Count dist_v_without = BFSv_without.distance(n);

                //dist_to_v_with.push_back(dist_v_with);
                //dist_to_v_without.push_back(dist_v_without);
                dist_to_v_with[idx] = dist_v_with;
                dist_to_v_without[idx] = dist_v_without;
                Count dist_change = std::abs(static_cast<int>(dist_v_with - dist_v_without));
                dist_change_v += dist_change; 
                max_dist_change_v = (max_dist_change_v < dist_change) ? dist_change : max_dist_change_v;
                N_dist_change_v = (dist_change > 0) ? N_dist_change_v + 1 :N_dist_change_v;
            }
            ++idx;
            if (BFSu_without.numberOfPaths(n) > 0 && BFSv_without.numberOfPaths(n) > 0){

                sum_diff_dist_u_v_without += std::abs(static_cast<int>(BFSu_without.distance(n) - BFSv_without.distance(n)));
            }
         //}
          //if (use_nonLinear) {
              // pagerank variations
            pagerank_variation += std::abs(static_cast<int>(pagerank_with.score(n) - pagerank_without.score(n)));
            pagerank_max_with = (pagerank_max_with > pagerank_with.score(n))? pagerank_max_with : pagerank_with.score(n);
            pagerank_max_without = (pagerank_max_without > pagerank_without.score(n))? pagerank_max_without : pagerank_without.score(n);

          //}

            // Core number variations
            core_variation += std::abs(static_cast<int>(core_with.score(n) - core_without.score(n)));
        }

        // Component sizes 
        std::map<index, Count> components_without_sizes = components_without.getComponentSizes();
        Count max_size_component_without = 0;
        for (std::map<index, Count>::iterator it = components_without_sizes.begin(); it != components_without_sizes.end(); ++it) {
            Count component_size = it->second;
            if (component_size > max_size_component_without) {
                max_size_component_without = component_size;
            }

        }
        // Components of u and v TODO can be deduced from previous BFS for nodes, ... for links .. ? 
        std::vector<std::vector<node>> all_components_without = components_without.getComponents();
        std::unordered_set<node> component_u_nodes_without(all_components_without[components_without.componentOfNode(u_main)].begin(), all_components_without[components_without.componentOfNode(u_main)].end());
            std::unordered_set<node> component_v_nodes_without(all_components_without[components_without.componentOfNode(v_main)].begin(), all_components_without[components_without.componentOfNode(v_main)].end());

        //  Graph u_component_without = GraphTools::subgraphFromNodes(history.main_graph, component_u_nodes_without);
        //  Graph v_component_without = GraphTools::subgraphFromNodes(history.main_graph, component_v_nodes_without);
        //*std::max_element(dist_to_u_with.begin(), dist_to_u_with.end());
        //std::max_element(dist_to_v_with.begin(), dist_to_v_with.end());
        //dist_to_u_without.begin(), dist_to_u_without.end());
        //dist_to_v_without.begin(), dist_to_v_without.end());

        integerResults["eccentricity u in G"] = *std::max_element(dist_to_u_with, dist_to_u_with + N);
        integerResults["eccentricity v in G"] = *std::max_element(dist_to_v_with, dist_to_v_with +N);
        integerResults["eccentricity u in G-"] = *std::max_element(dist_to_u_without, dist_to_u_without+N);

        integerResults["eccentricity v in G-"] = *std::max_element(dist_to_v_without, dist_to_v_without+N);
        integerResults["total distance change u"] = dist_change_u;
        integerResults["total distance change v"] = dist_change_v;
        integerResults["max distance change u"] = max_dist_change_u;
        integerResults["max distance change v"] = max_dist_change_v;
        integerResults["number of distance change u"] = N_dist_change_u;
        integerResults["number of distance change v"] = N_dist_change_v;
        integerResults["sum difference of distances to u and v in G-"] = sum_diff_dist_u_v_without;
        integerResults["sum difference of distances to u and v in G"] = sum_diff_dist_u_v_with;

        Count sum_distance_u_with = std::accumulate(dist_to_u_with, dist_to_u_with+N, 0);
        Count sum_distance_v_with = std::accumulate(dist_to_v_with, dist_to_v_with+N, 0);
        Count sum_distance_u_without = std::accumulate(dist_to_u_without, dist_to_u_without+N, 0);
        Count sum_distance_v_without = std::accumulate(dist_to_v_without, dist_to_v_without+N, 0);
        
        integerResults["closeness u G"] = (sum_distance_u_with > 0) ? 1 / sum_distance_u_with : 0;        
        integerResults["closeness v G"] = (sum_distance_v_with > 0) ? 1 / sum_distance_v_with : 0;        
        integerResults["closeness u G-"] = (sum_distance_u_without > 0) ? 1 / sum_distance_u_without : 0;        
        integerResults["closeness v G-"] = (sum_distance_v_without > 0) ? 1 / sum_distance_v_without : 0;        
        integerResults["distance u v in G-"] = BFSu_without.distance(v_main);

        std::vector<double> core_with_scores = core_with.scores();
        std::vector<double> core_without_scores = core_without.scores();

        doubleResults["degeneracy G"] = *std::max_element(core_with_scores.begin(), core_with_scores.end());
        doubleResults["degeneracy G-"] = *std::max_element(core_without_scores.begin(), core_without_scores.end());
        doubleResults["sum core scores G"] = std::accumulate(core_with_scores.begin(), core_with_scores.end(), 0);
        doubleResults["sum core scores G-"] = std::accumulate(core_without_scores.begin(), core_without_scores.end(), 0);
        doubleResults["core number u in G"] = core_with.score(u_main);
        doubleResults["core number v in G"] = core_with.score(v_main);
        doubleResults["core number u in G-"] = core_without.score(u_main);
        doubleResults["core number v in G-"] = core_without.score(v_main);
        doubleResults["sum core variations"] = core_variation;

        integerResults["number of components G-"] = components_without.numberOfComponents();
        integerResults["size largest component G-"] = max_size_component_without;

        integerResults["link component number of nodes"] = link_component.numberOfNodes();
        integerResults["link component number of links"] = link_component.numberOfEdges();

        doubleResults["pagerank of u in G"] = pagerank_with.score(u_main);
        doubleResults["pagerank of v in G"] = pagerank_with.score(v_main);
        doubleResults["pagerank of u in G-"] = pagerank_without.score(u_main);
        doubleResults["pagerank of v in G-"] = pagerank_without.score(v_main);
        doubleResults["pagerank max G"] = pagerank_max_with;
        doubleResults["pagerank max G-"] = pagerank_max_without;
        doubleResults["pagerank variation"] = pagerank_variation;


        //}
        //if (use_nonLinear) { 
        /**********/
        /*PageRank*/
        /**********/

        //}
        /*************/
        /*PLM on proj*/
        /*************/
        //std::vector<PLM> plms_top;
        //std::vector<PLM> plms_bot;
        //std::vector<Partition> partitions_top(N_plm, Partition(0));
        //std::vector<Count> N_changes_top(N_plm, 0);
        //std::vector<Partition> partitions_bot(N_plm, Partition(0));
        //std::vector<Count> N_changes_bot(N_plm, 0);

        //std::vector<Count> N_subsets_bot(N_plm, 0);
        //std::vector<Count> max_subset_size_bot(N_plm, 0);
        //std::vector<bool> same_subset_uv_bot(N_plm, false);
        //std::vector<Count> u_subset_size_bot(N_plm, 0);
        //std::vector<Count> v_subset_size_bot(N_plm, 0);
        //std::vector<Count> N_subsets_top(N_plm, 0);
        //std::vector<Count> max_subset_size_top(N_plm, 0);
        //std::vector<bool> same_subset_uv_top(N_plm, false);
        //std::vector<Count> u_subset_size_top(N_plm, 0);
        //std::vector<Count> v_subset_size_top(N_plm, 0);

        //// parallelize PLM
        //#pragma omp parallel for num_threads(1)
        //for (omp_index omp_idx = 0; omp_idx < static_cast<omp_index>(N_plm); ++omp_idx){
        //    //std::pair<Partition, count> estimation = estimate(partitions[u]);
        //    //PLM plm(projection, false, 1.0, "balanced",32,true,true,partitions[u]);
        //    // top graph
        //    Count z_top = history.top_graph.upperNodeIdBound();
        //    Partition zeta_top(z_top);
        //    zeta_top.allToSingletons();
        //    PLM plm_top(history.top_graph, false, 1.0, "none randomized",32,false,false,zeta_top);
        //    plm_top.run();
        //    partitions_top[omp_idx] = plm_top.getPartition();
        //    N_changes_top[omp_idx] = plm_top.getNumberChanges();


        //    Partition partition_top = plm_top.getPartition();
        //    partitions_top[omp_idx] = partition_top;
        //    N_subsets_top[omp_idx] = partition_top.numberOfSubsets();
        //    std::vector<Count> subset_sizes_top = partition_top.subsetSizes(); 
        //    max_subset_size_top[omp_idx] = *std::max_element(subset_sizes_top.begin(), subset_sizes_top.end());
        //    same_subset_uv_top[omp_idx] = partition_top.inSameSubset(u_main, v_main);
        //    index u_subset_top = partition_top.subsetOf(u_main);
        //    index v_subset_top = partition_top.subsetOf(v_main);
        //    u_subset_size_top[omp_idx] = partition_top.subsetSizeMap()[u_subset_top];
        //    v_subset_size_top[omp_idx] = partition_top.subsetSizeMap()[v_subset_top];
        //    N_changes_top[omp_idx] = plm_top.getNumberChanges();

        //    // bot graph
        //    Count z_bot = history.bot_graph.upperNodeIdBound();
        //    Partition zeta_bot(z_bot);
        //    zeta_bot.allToSingletons();
        //    PLM plm_bot(history.bot_graph, false, 1.0, "none randomized",32,false,false,zeta_bot);
        //    plm_bot.run();
        //    partitions_bot[omp_idx] = plm_bot.getPartition();
        //    N_changes_bot[omp_idx] = plm_bot.getNumberChanges();


        //    Partition partition_bot = plm_bot.getPartition();
        //    partitions_bot[omp_idx] = partition_bot;
        //    N_subsets_bot[omp_idx] = partition_bot.numberOfSubsets();
        //    std::vector<Count> subset_sizes_bot = partition_bot.subsetSizes(); 
        //    max_subset_size_bot[omp_idx] = *std::max_element(subset_sizes_bot.begin(), subset_sizes_bot.end());
        //    same_subset_uv_bot[omp_idx] = partition_bot.inSameSubset(u_main, v_main);
        //    index u_subset_bot = partition_bot.subsetOf(u_main);
        //    index v_subset_bot = partition_bot.subsetOf(v_main);
        //    u_subset_size_bot[omp_idx] = partition_bot.subsetSizeMap()[u_subset_bot];
        //    v_subset_size_bot[omp_idx] = partition_bot.subsetSizeMap()[v_subset_bot];
        //    N_changes_bot[omp_idx] = plm_bot.getNumberChanges();

        //}
        ////// sort partitions according to the number of subsets and output features
        ////// using this order
        //std::vector<index> sorting_indexes_top(N_plm, 0);
        //std::iota(sorting_indexes_top.begin(), sorting_indexes_top.end(), 0);
        //sort(sorting_indexes_top.begin(),
        //     sorting_indexes_top.end(),
        //     [&](int i,int j){return N_subsets_top[i]<N_subsets_top[j];});
        //std::vector<index> sorting_indexes_bot(N_plm, 0);
        //std::iota(sorting_indexes_bot.begin(), sorting_indexes_top.end(), 0);
        //sort(sorting_indexes_bot.begin(),
        //     sorting_indexes_bot.end(),
        //     [&](int i,int j){return N_subsets_bot[i]<N_subsets_bot[j];});

        //for (size_t idx=0; idx < sorting_indexes_top.size(); ++idx) {
        //    integerResults["number of subsets in top graph" + std::to_string(idx)] = N_subsets_top[sorting_indexes_top[idx]];
        //    integerResults["max subset size in top graph" + std::to_string(idx)] = max_subset_size_top[sorting_indexes_top[idx]];
        //    integerResults["u subset size in top graph" + std::to_string(idx)] = u_subset_size_top[sorting_indexes_top[idx]];
        //    integerResults["v  subset size in top graph" + std::to_string(idx)] = v_subset_size_top[sorting_indexes_top[idx]];
        //    integerResults["u v in same community in top graph"+ std::to_string(idx)] = same_subset_uv_top[sorting_indexes_top[idx]];

        //    integerResults["number of subsets in bot graph" + std::to_string(idx)] = N_subsets_bot[sorting_indexes_bot[idx]];
        //    integerResults["max subset size in bot graph" + std::to_string(idx)] = max_subset_size_bot[sorting_indexes_bot[idx]];
        //    integerResults["u subset size in bot graph" + std::to_string(idx)] = u_subset_size_bot[sorting_indexes_bot[idx]];
        //    integerResults["v  subset size in bot graph" + std::to_string(idx)] = v_subset_size_bot[sorting_indexes_bot[idx]];
        //    integerResults["u v in same community in bot graph"+ std::to_string(idx)] = same_subset_uv_bot[sorting_indexes_bot[idx]];

        //}
 
        /***********/
        /*PLM on G-*/
        /***********/
        // init N_PLM louvains 
        auto t0 = high_resolution_clock::now();

        ////PLM* plms_without = new PLM[N_plm]{0};
        std::vector<PLM> plms_without;
        std::vector<Partition> partitions_without(N_plm, Partition(0));
        std::vector<Count> N_changes_without(N_plm, 0);
        std::vector<Partition> partitions_with(N_plm, Partition(0));
        std::vector<Count> N_changes_with(N_plm, 0);

        std::vector<PLM> plms_with;


        std::vector<Count> N_subsets_with(N_plm, 0);
        std::vector<Count> max_subset_size_with(N_plm, 0);
        std::vector<bool> same_subset_uv_with(N_plm, false);
        std::vector<Count> u_subset_size_with(N_plm, 0);
        std::vector<Count> v_subset_size_with(N_plm, 0);
        std::vector<Count> N_subsets_without(N_plm, 0);
        std::vector<Count> max_subset_size_without(N_plm, 0);
        std::vector<bool> same_subset_uv_without(N_plm, false);
        std::vector<Count> u_subset_size_without(N_plm, 0);
        std::vector<Count> v_subset_size_without(N_plm, 0);

        // parallelize PLM
        #pragma omp parallel for num_threads(1)
        for (omp_index omp_idx = 0; omp_idx < static_cast<omp_index>(N_plm); ++omp_idx){
            //std::pair<Partition, count> estimation = estimate(partitions[u]);
            //PLM plm(projection, false, 1.0, "balanced",32,true,true,partitions[u]);
            Count z = history.main_graph.upperNodeIdBound();
            Partition zeta(z);
            zeta.allToSingletons();
            PLM plm_without(history.main_graph, false, 1.0, "none randomized",32,false,false,zeta);
            plm_without.run();
            partitions_without[omp_idx] = plm_without.getPartition();
            N_changes_without[omp_idx] = plm_without.getNumberChanges();


            Partition partition_without = plm_without.getPartition();
            partitions_without[omp_idx] = partition_without;
            N_subsets_without[omp_idx] = partition_without.numberOfSubsets();
            std::vector<Count> subset_sizes = partition_without.subsetSizes(); 
            max_subset_size_without[omp_idx] = *std::max_element(subset_sizes.begin(), subset_sizes.end());
            same_subset_uv_without[omp_idx] = partition_without.inSameSubset(u_main, v_main);
            index u_subset = partition_without.subsetOf(u_main);
            index v_subset = partition_without.subsetOf(v_main);
            u_subset_size_without[omp_idx] = partition_without.subsetSizeMap()[u_subset];
            v_subset_size_without[omp_idx] = partition_without.subsetSizeMap()[v_subset];
            N_changes_without[omp_idx] = plm_without.getNumberChanges();

        }

        //// sort partitions according to the number of subsets and output features
        //// using this order
        std::vector<index> sorting_indexes_without(N_plm, 0);
        std::iota(sorting_indexes_without.begin(), sorting_indexes_without.end(), 0);
        sort(sorting_indexes_without.begin(),
             sorting_indexes_without.end(),
             [&](int i,int j){return N_subsets_without[i]<N_subsets_without[j];});

        for (size_t idx=0; idx < sorting_indexes_without.size(); ++idx) {
            integerResults["number of subsets in G- " + std::to_string(idx)] = N_subsets_without[sorting_indexes_without[idx]];
            integerResults["max subset size in G- " + std::to_string(idx)] = max_subset_size_without[sorting_indexes_without[idx]];
            integerResults["u subset size in G- " + std::to_string(idx)] = u_subset_size_without[sorting_indexes_without[idx]];
            integerResults["v  subset size in G- " + std::to_string(idx)] = v_subset_size_without[sorting_indexes_without[idx]];
            integerResults["u v in same community in G- "+ std::to_string(idx)] = same_subset_uv_without[sorting_indexes_without[idx]];
        }
        //}


        // Restore (u,v)
        history.main_graph.addEdge(u_main, v_main);
        history.main_graph.setWeight(u_main, v_main, history.counter[e]);

        //if (use_nonLinear) {
        /**********/
        /*PLM on G*/
        /**********/
        //TODO TODO PLM ON PROJECTION !
        #pragma omp parallel for num_threads(N_plm)
        for (omp_index omp_idx = 0; omp_idx < static_cast<omp_index>(N_plm); ++omp_idx){
            PLM plm_with(history.main_graph, false, 1.0, "none randomized",32,false,true,partitions_without[omp_idx]);
            //std::cout
            plm_with.run();
            Partition partition_with = plm_with.getPartition();
            partitions_with[omp_idx] = partition_with;
            N_subsets_with[omp_idx] = partition_with.numberOfSubsets();
            std::vector<Count> subset_sizes = partition_with.subsetSizes(); 
            max_subset_size_with[omp_idx] = *std::max_element(subset_sizes.begin(), subset_sizes.end());
            same_subset_uv_with[omp_idx] = partition_with.inSameSubset(u_main, v_main);
            index u_subset = partition_with.subsetOf(u_main);
            index v_subset = partition_with.subsetOf(v_main);
            u_subset_size_with[omp_idx] = partition_with.subsetSizeMap()[u_subset];
            v_subset_size_with[omp_idx] = partition_with.subsetSizeMap()[v_subset];
            N_changes_with[omp_idx] = plm_with.getNumberChanges();
        }

        std::sort(N_changes_with.begin(), N_changes_with.end());
        // sort partitions according to the number of subsets and output features
        // using this order
        std::vector<index> sorting_indexes_with(N_plm, 0);
        std::iota(sorting_indexes_with.begin(), sorting_indexes_with.end(), 0);
        sort(sorting_indexes_with.begin(),
             sorting_indexes_with.end(),
             [&](int i,int j){return N_subsets_without[i]<N_subsets_without[j];});

        // get number of nodes that changed partition
        std::vector<Count> N_nodes_same_subset(N_plm, 0);
        for (size_t idx=0; idx < partitions_with.size(); ++idx) {
            Count N_nodes;
            for (const auto n : history.top_graph.nodeRange()){
                index n_subset_with = partitions_with[idx].subsetOf(n);
                index n_subset_without = partitions_without[idx].subsetOf(n);
                N_nodes = (n_subset_without == n_subset_with)? N_nodes + 1 : N_nodes;
            }
            N_nodes_same_subset[idx] = N_nodes;
        }
        for (size_t idx=0; idx < sorting_indexes_with.size(); ++idx) {
            integerResults["number of subsets in G " + std::to_string(idx)] = N_subsets_without[sorting_indexes_with[idx]];
            integerResults["max subset size in G " + std::to_string(idx)] = max_subset_size_without[sorting_indexes_with[idx]];
            integerResults["u subset size in G " + std::to_string(idx)] = u_subset_size_without[sorting_indexes_with[idx]];
            integerResults["v  subset size in G " + std::to_string(idx)] = v_subset_size_without[sorting_indexes_with[idx]];
            integerResults["u v in same community in G "+ std::to_string(idx)] = same_subset_uv_without[sorting_indexes_with[idx]];
            integerResults["number of nodes changing partition " + std::to_string(idx)] = N_nodes_same_subset[sorting_indexes_with[idx]]; 
        }

    }
    //std::ofstream myfile;
    Bound window = history.getWindow();
    //std::string filename = "cppOutput_" + std::to_string(window) + ".csv";
    //myfile.open (filename, std::ios::app);
    std::string output = "";

    if  (!header_done) {

        for(std::map<std::string,Count>::iterator it = integerResults.begin(); it != integerResults.end(); ++it){
            std::string k =  it->first;
            Count v = it->second;
            //std::cout << k << ":" <<v <<",";
            output += k + ",";
        }
        for(std::map<std::string,double>::iterator it = doubleResults.begin(); it != doubleResults.end(); ++it){
            std::string k =  it->first;
            double v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            //myfile <<k <<",";
            output += k + ",";
        }
    }

    ////myfile <<"\n";
    output += "\n";
    header_done = true;
    //}
    for(std::map<std::string,Count>::iterator it = integerResults.begin(); it != integerResults.end(); ++it){
            std::string k =  it->first;
            Count v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            //myfile <<v <<",";
            output += std::to_string(v) + ",";
        }
    for(std::map<std::string,double>::iterator it = doubleResults.begin(); it != doubleResults.end(); ++it){
            std::string k =  it->first;
            double v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            //myfile <<v <<",";
            output += std::to_string(v) + ",";
        }

    ////myfile << "\n";
    //output += "\n";
    return output;

}
}
