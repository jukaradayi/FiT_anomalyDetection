//
// Created by Apple on 30/05/2022.
//

#ifndef FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP
#define FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP


#include "globals.hpp"
#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
#include <vector>
using namespace NetworKit;
namespace StreamGraphs {
    class ConstantTimeMax{
    public:
        std::vector<uint64_t> deg_seq;
        std::vector<uint64_t> node2pos;
        std::vector<uint64_t> pos2node;
        std::vector<uint64_t> deg2pos;
        std::unordered_set<uint64_t> deg_set;
        std::vector<uint64_t> deg_distrib;
        std::unordered_set<uint64_t> Nodes;

        ConstantTimeMax();
        ~ConstantTimeMax();
        void addNode(node u);
        void increaseNodeDeg(node u);
        void decreaseNodeDeg(node u);
    };

}
#endif //FIT_ANOMALYDETECTION_SORTED_DEGREE_HPP
