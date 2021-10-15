/*
 *  metrics.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
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
    Metrics(HistoryGraph& history, bool use_basic, bool use_local, bool use_global, bool use_nonLinear); 

    void init();

    ~Metrics();

    void run(node u, node v); // run function stored in function-vectors

    inline Count degree(node u){node u_main = history.node2main[u]; return history.main_graph.degree(u_main);};

    //Count degree_v();
    inline Count top_degree(node u){node u_top = history.node2top[u];  return history.top_graph.degree(u_top);};
    inline Count bot_degree(node v){node v_bot = history.node2bot[v]; return history.bot_graph.degree(v_bot);};
    inline Count weighted_degree(node u){node u_main = history.node2main[u]; return history.main_graph.weightedDegree(u_main);};
    inline Count top_weighted_degree(node u){node u_top = history.node2top[u]; return history.top_graph.weightedDegree(u_top);};

    inline Count bot_weighted_degree(node v){node v_bot = history.node2bot[v]; return history.bot_graph.weightedDegree(v_bot);};

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

    inline Count max_weighted_degree(){return *std::max_element(history.main_weightedDegree_sequence.begin(),history.main_weightedDegree_sequence.end());}
    
    inline Count top_max_weighted_degree(){return *std::max_element(history.top_weightedDegree_sequence.begin(),history.top_weightedDegree_sequence.end());};
    
    inline Count bot_max_weighted_degree(){return *std::max_element(history.bot_weightedDegree_sequence.begin(),history.bot_weightedDegree_sequence.end());};

    inline Count degree_absolute_difference(node u, node v){
        node u_main = history.node2main[u];
        node v_main = history.node2main[v];
        return std::abs(history.main_graph.degree(u_main) - history.main_graph.degree(v_main));
    };

    inline Count weighted_degree_absolute_difference(node u, node v){
        node u_main = history.node2main[u];
        node v_main = history.node2main[v];
        return std::abs(history.main_graph.weightedDegree(u_main) - history.main_graph.weightedDegree(v_main));
    };

    //double average_weight();
    Count total_weight(){
        Count total_weight = 0;
        for(std::map<Edge,Count>::iterator it = history.counter.begin(); it != history.counter.end(); ++it){
            //Key k =  iter->first;
            Count v = it->second;
            total_weight += v;
        }
        return total_weight;
    };
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
    std::map<std::string, Count> integerResults;
    std::map<std::string, double> doubleResults;

    HistoryGraph& history;
    bool header_done; // TODO remove at end
 
private:

    // window size
    Bound window;
    bool use_basic;
    bool use_local;
    bool use_global;
    bool use_nonLinear;

    
};
}
#endif
