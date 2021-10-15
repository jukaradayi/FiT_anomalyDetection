// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "G_graph.hpp"
#include "history_graph.hpp"
#include "metrics.hpp"
#include "PLM.hpp"
#include <omp.h>
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

Metrics::Metrics(HistoryGraph& history, bool use_basic, bool use_local, bool use_global, bool use_nonLinear): history(history), use_basic(use_basic), use_local(use_local), use_global(use_global), use_nonLinear(use_nonLinear), header_done(false) {}

Metrics::~Metrics() {}

void Metrics::run(node u, node v) {

    // initialisations
    Edge e(u, v);
    node u_main = history.node2main[u];
    node u_top = history.node2top[u];
    node v_main = history.node2main[v];
    node v_bot = history.node2bot[v];



        

    if (use_basic) {
        /***************/
        /*Basic metrics*/
        /***************/
        integerResults["number of nodes"] = history.main_graph.numberOfNodes(); 
        integerResults["number of links"] = history.main_graph.numberOfEdges();
        integerResults["count nodes degree 1"] = history.main_degree_distribution[0];
        integerResults["count nodes degree 2"] = history.main_degree_distribution[1];
        integerResults["degree u"] = degree(u);
        integerResults["degree v"] = degree(v);
        integerResults["weighted degree u"] = weighted_degree(u);
        integerResults["weighted degree v"] = weighted_degree(v);
        integerResults["max degree"] = max_degree();
        integerResults["max weighted degree"] = max_weighted_degree();
        integerResults["degree absolute difference"] = degree_absolute_difference(u, v);
        integerResults["weighted degree absolute difference"] = weighted_degree_absolute_difference(u, v);
        integerResults["weight"] = history.counter[e];

        /********************/
        /*Projection metrics*/
        /********************/
        integerResults["number of top nodes"] = history.top_graph.numberOfNodes();
        integerResults["number of bot nodes"] = history.bot_graph.numberOfNodes();
        integerResults["number of top links"] = history.top_graph.numberOfEdges();
        integerResults["number of bot links"] = history.bot_graph.numberOfEdges();
        integerResults["total weight"] = total_weight();
        integerResults["top degree"] = top_degree(u);
        integerResults["bot degree"] = bot_degree(v);
        integerResults["top weighted degree"] = top_weighted_degree(u);
        integerResults["bot weighted degree"] = bot_weighted_degree(v);
        integerResults["top max weighted degree"] = top_max_weighted_degree();
        integerResults["bot max weighted degree"] = bot_max_weighted_degree();
    };

    if (use_local) {
        //// sets
        // get distance 1 and 2 neighbor sets, basic sets
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

        // fill sets
        for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, u_main, history.main_bound).begin();
            bN != BoundedNeighborRange(history.main_graph, u_main, history.main_bound).end(); ++bN) {
            node n_main = *bN;
            bNu.insert(n_main);
            bNu_u_bNNv.insert(n_main);

            for (BoundedNeighborIterator bNN = BoundedNeighborRange(history.main_graph, n_main, history.main_bound).begin();
                bNN != BoundedNeighborRange(history.main_graph, n_main, history.main_bound).end(); ++bNN) {
                node nn_main = *bNN;
                bNNu.insert(nn_main);
                bNNu_u_bNv.insert(nn_main);
            }
        }

        for (BoundedNeighborIterator bN = BoundedNeighborRange(history.main_graph, v_main, history.main_bound).begin();
            bN != BoundedNeighborRange(history.main_graph, v_main, history.main_bound).end(); ++bN) {
            node n_main = *bN;
            bNv.insert(n_main);
            //std::pair<std::unordered_set<node>,bool> insertion = bNNu_u_bNv.insert(n_main);
            auto insertion = bNNu_u_bNv.insert(n_main);


            //if (bNNu.find(n_main) != bNNu.end()) {
            if (insertion.second == false){
                bNNu_n_bNv.insert(n_main);
            }
            for (BoundedNeighborIterator bNN = BoundedNeighborRange(history.main_graph, n_main, history.main_bound).begin();
                bNN != BoundedNeighborRange(history.main_graph, n_main, history.main_bound).end(); ++bNN) {
                node nn_main = *bNN;
                bNNv.insert(nn_main);
                //std::pair<std::unordered_set<node>,bool> insertion = bNu_u_bNNv.insert(nn_main);
                auto insertion = bNu_u_bNNv.insert(nn_main);


                //if (bNu.find(nn_main) != bNu.end()) {
                if (insertion.second == false){
                    bNu_n_bNNv.insert(nn_main);
                }
            }
        }

        // build egonets

        // egonet sets

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

        /******************/
        /*local clustering*/
        /******************/
        LocalClusteringCoefficient clustering(history.main_graph, true);
        clustering.run();

        doubleResults["clustering u"] = clustering.score(u_main);
        doubleResults["clustering v"] = clustering.score(v_main);

        LocalClusteringCoefficient top_clustering(history.top_graph, true);
        top_clustering.run();
        LocalClusteringCoefficient bot_clustering(history.bot_graph, true);
        bot_clustering.run();

        doubleResults["top clustering u"] = top_clustering.score(u_top);
        doubleResults["bot clustering v"] = bot_clustering.score(v_bot);

        /*********************************/
        /*Jaccard index and neighborhoods*/
        /*********************************/
        JaccardIndex jaccard(history.main_graph);
        doubleResults["Jaccard index u v"] = jaccard.run(u_main, v_main);
        doubleResults["Jaccard index bipartite u"] = std::abs(bNu_u_bNNv.size() / bNu_n_bNNv.size());
        doubleResults["Jaccard index bipartite v"] = std::abs(bNNu_u_bNv.size() / bNNu_n_bNv.size());

        doubleResults["ratio neighbors of u in N(v)"] = bNu_n_bNv.size() / bNu.size();
        doubleResults["ratio neighbors of v in N(u)"] = bNu_n_bNv.size() / bNv.size();
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
    }
    //if (use_global) { 
        /************/
        /*components*/
        /************/
    std::vector<CoreDecomposition> all_coreDecomposition;

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
    std::vector<BFS> all_BFS;


        all_BFS.push_back(BFS(history.main_graph, u_main, true));
        all_BFS.push_back(BFS(history.main_graph, v_main, true));

        /********************/
        /*Core Decomposition*/
        /********************/
        //std::vector<CoreDecomposition> all_coreDecomposition;
        all_coreDecomposition.push_back(CoreDecomposition(history.main_graph, false));

    //}
    //if (use_nonLinear) {
        /**********/
        /*PageRank*/
        /**********/
        //std::vector<PageRank> all_pagerank;
        std::vector<PageRank> all_pagerank;

        all_pagerank.push_back(PageRank(history.main_graph));
    //}
    /************************/
    /************************/
    /* Compute metrics on G-*/
    /************************/
    /************************/
    // Remove (u,v)
    history.main_graph.removeEdge(u_main, v_main);

    //if (use_global) {
        /*****/
        /*BFS*/
        /*****/
        all_BFS.push_back(BFS(history.main_graph, u_main, true));
        all_BFS.push_back(BFS(history.main_graph, v_main, true));

        #pragma omp parallel for
        for (omp_index u = 0; u < static_cast<omp_index>(4); ++u){
            all_BFS[u].run();
        }

        BFS BFSu_with = all_BFS.back(); 
        all_BFS.pop_back();
        BFS BFSv_with = all_BFS.back();
        all_BFS.pop_back();
        BFS BFSu_without = all_BFS.back();
        all_BFS.pop_back();
        BFS BFSv_without = all_BFS.back();
        all_BFS.pop_back();

        /********************/
        /*Core Decomposition*/
        /********************/
        all_coreDecomposition.push_back(CoreDecomposition(history.main_graph, false));

        #pragma omp parallel for
        for (omp_index u = 0; u < static_cast<omp_index>(2); ++u){
            all_coreDecomposition[u].run();
        }

        CoreDecomposition core_with = all_coreDecomposition.back();
        all_coreDecomposition.pop_back();

        CoreDecomposition core_without = all_coreDecomposition.back();
        all_coreDecomposition.pop_back();

        /**********/
        /*PageRank*/
        /**********/
        all_pagerank.push_back(PageRank(history.main_graph));
        //PageRank pagerank_without(history.main_graph);
        //pagerank_without.run();
        std::cout << "PR \n";
        #pragma omp parallel for
        for (omp_index u = 0; u < static_cast<omp_index>(2); ++u){
            all_pagerank[u].run();
        }
        PageRank pagerank_with = all_pagerank.back();
        all_pagerank.pop_back();
        PageRank pagerank_without = all_pagerank.back();
        all_pagerank.pop_back();


        /**********************/
        /*Connected Components*/
        /**********************/
        // TODO run or just deduce from BFS ? ..
        ConnectedComponents components_without(history.main_graph);
        components_without.run();

        //# compute distance changes from u and v
        Count dist_change_u = 0;
        Count dist_change_v = 0;
        Count max_dist_change_u = 0;
        Count max_dist_change_v = 0;
        Count N_dist_change_u = 0;
        Count N_dist_change_v = 0;
        Count sum_diff_dist_u_v_with = 0;
        Count sum_diff_dist_u_v_without = 0;

        // keep list of distances to accessible nodes
        std::vector<Count> dist_to_u_with;
        std::vector<Count> dist_to_v_with;
        std::vector<Count> dist_to_u_without;
        std::vector<Count> dist_to_v_without;

        double pagerank_variation = 0;
        double pagerank_max_with = 0;
        double pagerank_max_without = 0;
        double core_variation = 0;

        for (const auto n : history.main_graph.nodeRange()){
            // Loop over all nodes to compute variations in several metrics
            //
            if (BFSu_without.numberOfPaths(n) > 0){
                Count dist_u_with = BFSu_with.distance(n);
                Count dist_u_without = BFSu_without.distance(n);

                dist_to_u_with.push_back(dist_u_with);
                dist_to_u_without.push_back(dist_u_without);

                Count dist_change = std::abs(dist_u_with - dist_u_without);
                dist_change_u += dist_change; 
                max_dist_change_u = (max_dist_change_u < dist_change) ? dist_change : max_dist_change_u;
                N_dist_change_u = (dist_change > 0) ? N_dist_change_u + 1 :N_dist_change_u;

                // if there's a path to u in G, there's one to v ..
                sum_diff_dist_u_v_with += std::abs(dist_u_with - BFSv_with.distance(n));

            }

            if (BFSv_without.numberOfPaths(n) > 0){
                Count dist_v_with = BFSv_with.distance(n);
                Count dist_v_without = BFSv_without.distance(n);

                dist_to_v_with.push_back(dist_v_with);
                dist_to_v_without.push_back(dist_v_without);

                Count dist_change = std::abs(dist_v_with - dist_v_without);
                dist_change_v += dist_change; 
                max_dist_change_v = (max_dist_change_v < dist_change) ? dist_change : max_dist_change_v;
                N_dist_change_v = (dist_change > 0) ? N_dist_change_v + 1 :N_dist_change_v;
            }
            if (BFSu_without.numberOfPaths(n) > 0 && BFSv_without.numberOfPaths(n) > 0){
                sum_diff_dist_u_v_without += std::abs(BFSu_without.distance(n) - BFSv_without.distance(n));
            }
            if (use_nonLinear) {
                // pagerank variations
                pagerank_variation += std::abs(pagerank_with.score(n) - pagerank_without.score(n));
                pagerank_max_with = (pagerank_max_with > pagerank_with.score(n))? pagerank_max_with : pagerank_with.score(n);
                pagerank_max_without = (pagerank_max_without > pagerank_without.score(n))? pagerank_max_without : pagerank_without.score(n);

            }

            // Core number variations
            core_variation += std::abs(core_with.score(n) - core_without.score(n));
            
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

        Graph u_component_without = GraphTools::subgraphFromNodes(history.main_graph, component_u_nodes_without);
        Graph v_component_without = GraphTools::subgraphFromNodes(history.main_graph, component_v_nodes_without);

        integerResults["eccentricity u in G"] = *std::max_element(dist_to_u_with.begin(), dist_to_u_with.end());
        integerResults["eccentricity v in G"] = *std::max_element(dist_to_v_with.begin(), dist_to_v_with.end());
        integerResults["eccentricity u in G-"] = *std::max_element(dist_to_u_without.begin(), dist_to_u_without.end());
        integerResults["eccentricity v in G-"] = *std::max_element(dist_to_v_without.begin(), dist_to_v_without.end());
        integerResults["total distance change u"] = dist_change_u;
        integerResults["total distance change v"] = dist_change_v;
        integerResults["max distance change u"] = max_dist_change_u;
        integerResults["max distance change v"] = max_dist_change_v;
        integerResults["number of distance change u"] = N_dist_change_u;
        integerResults["number of distance change v"] = N_dist_change_v;
        integerResults["sum difference of distances to u and v in G-"] = sum_diff_dist_u_v_without;
        integerResults["sum difference of distances to u and v in G"] = sum_diff_dist_u_v_with;

        Count sum_distance_u_with = std::accumulate(dist_to_u_with.begin(), dist_to_u_with.end(), 0);
        Count sum_distance_v_with = std::accumulate(dist_to_v_with.begin(), dist_to_v_with.end(), 0);
        Count sum_distance_u_without = std::accumulate(dist_to_u_without.begin(), dist_to_u_without.end(), 0);
        Count sum_distance_v_without = std::accumulate(dist_to_v_without.begin(), dist_to_v_without.end(), 0);
        
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

        integerResults["number of components G"] = components_without.numberOfComponents();
        integerResults["size largest component G"] = max_size_component_without;

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


        /***********/
        /*PLM on G-*/
        /***********/
        // init N_PLM louvains 
        int N_plm = 5;
        std::vector<PLM> plms_without;
        std::vector<Partition> partitions_without;
        std::vector<Count> N_changes_without;
        std::vector<Partition> partitions_with;
        std::vector<Count> N_changes_with;

        std::vector<PLM> plms_with;


        std::vector<Count> N_subsets_with;
        std::vector<Count> max_subset_size_with;
        std::vector<bool> same_subset_uv_with;
        std::vector<Count> u_subset_size_with;
        std::vector<Count> v_subset_size_with;
        std::vector<Count> N_subsets_without;
        std::vector<Count> max_subset_size_without;
        std::vector<bool> same_subset_uv_without;
        std::vector<Count> u_subset_size_without;
        std::vector<Count> v_subset_size_without;
        std::cout << "plm on G-\n";

        #pragma omp parallel for
        for (omp_index omp_idx = 0; omp_idx < static_cast<omp_index>(N_plm); ++omp_idx){
            //std::pair<Partition, count> estimation = estimate(partitions[u]);
            //PLM plm(projection, false, 1.0, "balanced",32,true,true,partitions[u]);
            Count z = history.main_graph.upperNodeIdBound();
            Partition zeta(z);
            zeta.allToSingletons();
            //std::cout << "init PLMs"<<u<<"\n";
            PLM plm_without(history.main_graph, false, 1.0, "none randomized",32,false,false,zeta);
            //std::cout << "run PLMs"<<u<<"\n";

            plm_without.run();
            //std::cout << "ran PLM" <<u "\n";
            partitions_without.push_back(plm_without.getPartition());
            N_changes_without.push_back(plm_without.getNumberChanges());


            Partition partition_without = plm_without.getPartition();
            partitions_without.push_back(partition_without);
            N_subsets_without.push_back(partition_without.numberOfSubsets());
            std::vector<Count> subset_sizes = partition_without.subsetSizes(); 
            max_subset_size_without.push_back(*std::max_element(subset_sizes.begin(), subset_sizes.end()));
            same_subset_uv_without.push_back(partition_without.inSameSubset(u_main, v_main));
            index u_subset = partition_without.subsetOf(u_main);
            index v_subset = partition_without.subsetOf(v_main);
            u_subset_size_without.push_back(partition_without.subsetSizeMap()[u_subset]);
            v_subset_size_without.push_back(partition_without.subsetSizeMap()[v_subset]);
            N_changes_without.push_back(plm_without.getNumberChanges());

        }
        std::cout << "plm on G- finished\n";

        // sort partitions according to the number of subsets and output features
        // using this order
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
    history.main_graph.setWeight(u_main, v_main, history. counter[e]);
    //if (use_nonLinear) {
        /**********/
        /*PLM on G*/
        /**********/

        std::cout << "plm on G\n";
        #pragma omp parallel for
        for (omp_index omp_idx = 0; omp_idx < static_cast<omp_index>(N_plm); ++omp_idx){
            PLM plm_with(history.main_graph, false, 1.0, "none randomized",32,false,true,partitions_without[omp_idx]);
            //std::cout
            plm_with.run();
            Partition partition_with = plm_with.getPartition();
            partitions_with.push_back(partition_with);
            N_subsets_with.push_back(partition_with.numberOfSubsets());
            std::vector<Count> subset_sizes = partition_with.subsetSizes(); 
            max_subset_size_with.push_back(*std::max_element(subset_sizes.begin(), subset_sizes.end()));
            same_subset_uv_with.push_back(partition_with.inSameSubset(u_main, v_main));
            index u_subset = partition_with.subsetOf(u_main);
            index v_subset = partition_with.subsetOf(v_main);
            u_subset_size_with.push_back(partition_with.subsetSizeMap()[u_subset]);
            v_subset_size_with.push_back(partition_with.subsetSizeMap()[v_subset]);
            N_changes_with.push_back(plm_with.getNumberChanges());
        }
        std::cout << "plm on G finished\n";

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
    //}




    //}
    std::ofstream myfile;
    Bound window = history.getWindow();
std::string filename = "cppOutput_" + std::to_string(window) + ".csv";
    myfile.open (filename, std::ios::app);

    if  (!header_done) {

        for(std::map<std::string,Count>::iterator it = integerResults.begin(); it != integerResults.end(); ++it){
            std::string k =  it->first;
            Count v = it->second;
            //std::cout << k << ":" <<v <<",";
            myfile << k <<",";
        }
        for(std::map<std::string,double>::iterator it = doubleResults.begin(); it != doubleResults.end(); ++it){
            std::string k =  it->first;
            double v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            myfile <<k <<",";
        }

    myfile <<"\n";
    header_done = true;
    }
    for(std::map<std::string,Count>::iterator it = integerResults.begin(); it != integerResults.end(); ++it){
            std::string k =  it->first;
            Count v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            myfile <<v <<",";
        }
    for(std::map<std::string,double>::iterator it = doubleResults.begin(); it != doubleResults.end(); ++it){
            std::string k =  it->first;
            double v = it->second;
            //std::ofstream myfile;
            //myfile.open ("cppOutput.csv", std::ios::app);

            myfile <<v <<",";
        }

    myfile << "\n";

}

//Count Metrics::number_of_nodes(){
//    std::cout << history.main_graph.numberOfNodes();
//    return history.main_graph.numberOfNodes();
//}

//Count Metrics::top_number_of_nodes(){
//    return history.top_graph.numberOfNodes();
//}
//
//Count Metrics::bot_number_of_nodes(){
//    return history.bot_graph.numberOfNodes();
//}
//
//Count Metrics::top_number_of_links(){
//    return top_graph.numberOfEdges();
//}
//
//Count Metrics::bot_number_of_links(){
//    return history.bot_graph.NumberOfEdges();
//}
//
//Count Metrics::number_node_degree_1(){
//    return history.degree_distribution[0];
//}
//
//Count Metrics::number_node_degree_2(){
//    return history.degree_distribution[1];
//}
//
//Count Metrics::degree_u() {
//    return history.main_graph.degree(u);
//}
//
//Count Metrics::degree_v() {
//    return history.main_graph.degree(v);
//}
//
//Count Metrics::top_degree_u() {
//    return history.top_graph.degree(u);
//}
//
//Count Metrics::bot_degree_v(){
//    return history.bot_graph.degree(v);
//}
//
//Count Metrics::weighted_degree_u(){
//    return history.main_graph.weightedDegree(u);
//}
//
//Count Metrics::weighted_degree_v(){
//    return history.main_graph.weightedDegree(v);
//}
//
//Count Metrics::top_weighted_degree_u(){
//    return history.top_graph.weightedDegree(u);
//}
//
//Count Metrics::bot_weighted_degree_v(){
//    return history.bot_graph.weightedDegree(v);
//}
//
//double Metrics::average_degree(){
//    return 1;
//}
//
//double Metrics::top_average_degree(){
//    return 1;
//}
//
//double Metrics::bot_average_degree(){
//    return 1;
//}
//
//Count Metrics::max_degree(){
//    return std::max_element(history.degree_distribution.begin(), history.degree_distribution.end());
//}
//
//Count Metrics::top_max_degree(){
//    return 1;
//}
//Count Metrics::bot_max_degree(){
//    return 1;
//}
//Count Metrics::max_weighted_degree(){
//    return 1
//}
//Count Metrics::top_max_weighted_degree(){
//    return 1;
//}
//Count Metrics::bot_max_weighted_degree(){
//    return 1;
//}
//
//Count Metrics::degree_absolute_difference(){
//    return history.main_graph.degree(u) - history.main_graph.degree(v);
//}
//Count Metrics::weighted_degree_absolute_difference(){
//    return history.main_graph.weightedDegree(u) - history.main_graph.weightedDegree(v);
//}



}
