/*
 * Compute features on static graphs.
 * Given a history graph H (resp. G) at a given interaction (t,u,v),
 * compute a set of metrics and write the results in the output csv.
 * The metrics are in three sets, each enabled by a boolean (use_basic, use_local and use_nonLinear). 
 * use_basic enables basic graph features, that we compute in O(1) time, like the degree of u and v in main_graph, the weight of (u,v) in main_graph etc...
 * use_local enables graph features computed in O(k) or O(k²) time, where k is the degree bound, which include the adamic adar coefficient, the Jaccard coefficient of u and v, the densities of various local neighborhoods of u and v (like N(u) \union N(v), N(N(u)) \union N(N(v)) etc...).
 */

#ifndef STREAMGRAPHS_METRICS_HPP_
#define STREAMGRAPHS_METRICS_HPP_

#include "globals.hpp"
#include "history_graph.hpp"
#include "G_graph.hpp"
#include "H_graph.hpp"


#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>

#include <limits>
#include <utility>
#include <cstdint>
#include <map>
#include <string>
#include <cmath>


using namespace NetworKit;
//using namespace StreamGraphs;
namespace StreamGraphs {

class Metrics final {

    // constructor & destructor

public:

    /**
     * Constructor of metrics class, used to compute features on history graphs.
     *
     * @briefCompute metrics on main_graph, top_graph and bot_graph.
     *
     * @param history (HistoryGraph) : history graph object on which we compute metrics
     * @param use_basic (bool) : enable to compute metrics in O(1) (number of nodes/links, total weight, degrees and weighted degrees of u and v, weight of (u,v) etc...)
     * @param use_local (bool) : enable to compute metrics in O(deg(u)+deg(v)) and O(deg(u)² + deg(v)²) (local clustering, Jaccard index of u and v, Adamic Adar coefficient, and tailored neighborhoods etc..)
     * @param use_nonLinear (bool) : enable to compute intensive metrics (PageRank, Louvain method, Distance based metrics, coreness etc...)
     */
    Metrics(HistoryGraph& history, const bool use_basic, const bool use_local, const bool use_global, const bool use_nonLinear); 

    void BiFS(const Graph &G, node source, bool storePaths); 

    void init();

    ~Metrics();

    /**
     * Run metric computations.
     *
     * @param u (node) : first node of the current link
     * @param v (node) : second node of the current link
     *
     */
    std::string run(const node u, const node v, const int interaction_id); // run function stored in function-vectors

    /**
     * return the degree of a node in the main_graph 
     */
    inline Count degree(node u){node u_main = history.node2main[u]; return history.main_graph.degree(u_main);};

    /**
     * return the degree of a node in the top_graph 
     */
    inline Count top_degree(node u){node u_top = history.node2top[u];  return history.top_graph.degree(u_top);};

    /**
     * return the degree of a node in the bot_graph 
     */
    inline Count bot_degree(node v){node v_bot = history.node2bot[v]; return history.bot_graph.degree(v_bot);};

    /**
     *return the weighted degree of a node in the main graph
     */
    inline Count weighted_degree(node u){node u_main = history.node2main[u]; return history.main_weightedDegree_sequence[u_main];};

    /**
     *return the weighted degree of a node in the top graph
     */
    inline Count top_weighted_degree(node u){node u_top = history.node2top[u]; return history.top_weightedDegree_sequence[u_top];};

    /**
     *return the weighted degree of a node in the bot graph
     */
    inline Count bot_weighted_degree(node v){node v_bot = history.node2bot[v]; return history.bot_weightedDegree_sequence[v_bot];};


    //inline Count weighted_degree(node u){node u_main = history.node2main[u]; return history.main_graph.weightedDegree(u_main);};
    //inline Count top_weighted_degree(node u){node u_top = history.node2top[u]; return history.top_graph.weightedDegree(u_top);};

    //inline Count bot_weighted_degree(node v){node v_bot = history.node2bot[v]; return history.bot_graph.weightedDegree(v_bot);};


    /**
     *return the maximum degree in the graph
     */
    inline Count max_degree(){
        Count max_degree = 0;
        for (size_t i=0; i < history.main_degree_distribution.size(); ++i) {
            if (history.main_degree_distribution[i] >0) {
                max_degree = i;
            }
        }
        return max_degree;
    };

    inline Count top_max_degree(){
        Count max_degree = 0;
        for (size_t i=0; i < history.top_degree_distribution.size(); ++i) {
            if (history.top_degree_distribution[i] >0) {
                max_degree = i;
            }
        }
        return max_degree;
    };

    inline Count bot_max_degree(){
        Count max_degree = 0;
        for (size_t i=0; i < history.bot_degree_distribution.size(); ++i) {
            if (history.bot_degree_distribution[i] >0) {
                max_degree = i;
            }
        }
        return max_degree;        
    };

    /**
     *return the maximum weighted degree in the graph
     */
    inline Count max_weighted_degree(){return *std::max_element(history.main_weightedDegree_sequence.begin(),history.main_weightedDegree_sequence.end());}
    
    inline Count top_max_weighted_degree(){return *std::max_element(history.top_weightedDegree_sequence.begin(),history.top_weightedDegree_sequence.end());};
    
    inline Count bot_max_weighted_degree(){return *std::max_element(history.bot_weightedDegree_sequence.begin(),history.bot_weightedDegree_sequence.end());};

    /**
     * return the absolute difference of the degrees of u and v
     */
    inline Count degree_absolute_difference(node u, node v){
        node u_main = history.node2main[u];
        node v_main = history.node2main[v];
        return std::abs(static_cast<int>(history.main_graph.degree(u_main) - history.main_graph.degree(v_main)));
    };

    /**
     * return the absolute difference of the weighted degrees of u and v
     */
    inline Count weighted_degree_absolute_difference(node u, node v){
        node u_main = history.node2main[u];
        node v_main = history.node2main[v];
        return std::abs(static_cast<int>(history.main_graph.weightedDegree(u_main) - history.main_graph.weightedDegree(v_main)));
    };

    /**
     * Return the sum of the weights of all the links in main_graph
     */
    Count total_weight(){
        //Count total_weight = 0;
        //for(std::map<Edge,Count>::iterator it = history.counter.begin(); it != history.counter.end(); ++it){
        //    //Key k =  iter->first;
        //    Count v = it->second;
        //    total_weight += v;
        //}
        //return total_weight;
        return history.getTotalWeight();
    };

    /**
     * Compute the local clustering of nodes u and v of current link.
     */
    std::pair<double, double> localClustering(const NetworKit::Graph& G, const node u, const node v);

    /**
     * Count the number of links between the nodes in a set of nodes.
     *
     * @param G (NetworKit::Graph) : graph on which we count the links
     * @param nodes (unordered_set) : set of nodes between which we count the links
     * @param bound (int) : bound on the degree of a node, past that bound, we don't count the links of this node.
     */ 
    int countLinks(const NetworKit::Graph& G, const std::unordered_set<node>& nodes, int bound);

    //double weight_concentration();
    //double difference_to_average_weight();
    //double weight_relative_u();
    //double weight_relative_v();
    //double fraction_total_weight();

    //double density();

    //Count redundancy_top();
    //Count redundancy_bot();

    //// Local features in O(k) cost

    //Count neighborhood_overlap();
    //Count neighborhood_size_top();
    //Count neighborhood_size_bot();

    //Count neighborhood_distance_index();

    /**
     * Compute the adamic adar coefficient.
     *
     * set adamic_adar=0
     * Loop through the nodes n of N(u) \inter N(v), set adamic_adar = adamic_adar + 1/log( degree(n) )
     * @param bNu_n_bNv (unordered_set:) intersection of the neighborhoods of u and v
     */
    inline double adamic_adar(std::unordered_set<node> bNu_n_bNv) {
        double adamic_adar = 0;
        for (std::unordered_set<node>::const_iterator it = bNu_n_bNv.begin(); it != bNu_n_bNv.end(); it++){
            // skip nodes of degree 1 -> can happen because u and v are in the sets..
            if (history.main_graph.degree(*it) <= 1) {
                continue;
            }
            adamic_adar += 1 / std::log(static_cast<double>(history.main_graph.degree(*it))); 
        }
        return adamic_adar; 

    };

    /**
     * Compute a bipartite version of the adamic adar coefficient.
     * Same as the "normal" version but on neighborhood N(u) \inter N(N(v))
     * @param bNu_n_bNNv (unordered_set:) intersection of the neighborhoods of u and N(v)
     */
    double adamic_adar_bipartite_u(std::unordered_set<node> bNu_n_bNNv){
        double adamic_adar = 0;
         for (std::unordered_set<node>::const_iterator it = bNu_n_bNNv.begin(); it != bNu_n_bNNv.end(); it++){
            if (history.main_graph.degree(*it) <= 1) {
                continue;
            }
            adamic_adar += 1 / std::log(static_cast<double>(history.main_graph.degree(*it))); 
        }
        return adamic_adar;
       
    };

    /**
     * Compute a bipartite version of the adamic adar coefficient.
     * Same as the "normal" version but on neighborhood N(N(u)) \inter N(v)
     * @param bNNu_n_bNv (unordered_set:) intersection of the neighborhoods of N(u) and v
     */
    double adamic_adar_bipartite_v(std::unordered_set<node> bNNu_n_bNv){
        double adamic_adar = 0;
        for (std::unordered_set<node>::const_iterator it = bNNu_n_bNv.begin(); it != bNNu_n_bNv.end(); it++){
            if (history.main_graph.degree(*it) <= 1) {
                continue;
            }
            adamic_adar += 1 / std::log(static_cast<double>(history.main_graph.degree(*it))); 
        }
        return adamic_adar; 
    };

    //// Louvain
    //Count top_number_communityChange_PLM();
    //Count top_number_communities_PLM_with();
    //Count top_number_communities_PLM_without();
    //Count top_max_community_size_PLM();
    //Count top_u_community_size_PLM_with();
    //Count top_u_community_size_PLM_without();
    //Count top_new_node_uCommunity_PLM();
    //Count top_leaving_node_uCommunity_PLM();
    //Count bot_number_communityChange_PLM();
    //Count bot_number_communities_PLM_with();
    //Count bot_number_communities_PLM_without();
    //Count bot_max_community_size_PLM();
    //Count bot_v_community_size_PLM_with();
    //Count bot_v_community_size_PLM_without();
    //Count bot_new_node_vCommunity_PLM();
    //Count bot_leaving_node_vCommunity_PLM();

    //Count common_neighbors_index();
    //Count common_neighbors_bipartite_u();
    //Count common_neighbors_bipartite_v();
    //Count same_community_index();
    //Count katz_index();

    // attributes
    std::map<std::string, Count> integerResults; ///< map that stores features that are integers
    std::map<std::string, double> doubleResults; ///< map that stores features that are floats

    HistoryGraph& history; ///< history graph
    bool header_done; ///< boolean, used to print the header
 
private:

    // window size
    Bound window; ///< window size of the history graph
    bool use_basic; ///< compute O(1) metrics
    bool use_local; ///< compute O(k + k²) metrics
    bool use_global;
    bool use_nonLinear; ///< compute non linear metrics

    
};
}
#endif
