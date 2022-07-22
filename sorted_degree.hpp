//
// Created by ZHEN HOU on 30/05/2022.
//

#ifndef FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP
#define FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP


#include "globals.hpp"
#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
#include <vector>
#include <unordered_map>
using namespace NetworKit;
namespace StreamGraphs {
    class ConstantTimeMax{
    public:
        //std::vector<uint64_t> deg_seq; // queue?
        std::deque<uint64_t> deg_seq;

        std::unordered_map<node,uint64_t> node2pos;
        //std::vector<uint64_t> pos2node; //queue?
        std::deque<uint64_t> pos2node;

        std::deque<uint64_t> deg2pos;
        std::unordered_set<uint64_t> deg_set;
        std::vector<uint64_t> deg_distrib;
        std::unordered_set<uint64_t> Nodes;

        ConstantTimeMax(std::unordered_set<uint64_t> Nodes);
        ConstantTimeMax();
        ~ConstantTimeMax();
        void addNode(node u);
        void increaseNodeDeg(node u);
        void decreaseNodeDeg(node u);
        std::string affichage();
    };

}
#endif //FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP
