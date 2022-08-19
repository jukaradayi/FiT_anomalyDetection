/*
 *  history_graph.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
 */

#ifndef STREAMGRAPHS_GLOBALS_HPP_
#define STREAMGRAPHS_GLOBALS_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
//#include "hash_graph.hpp"

#include <limits>
#include <utility>
#include <cstdint>

using namespace NetworKit;


namespace StreamGraphs {

/******************/
/* struct & types */
/******************/
using Time = uint64_t; // for timestamps
using Bound = uint64_t; // for limits
using Count = uint64_t; // for degree, number of things etc...
using edgeweight = uint64_t;
using edgeid = uint64_t;

constexpr int none = std::numeric_limits<uint64_t>::max(); // like in networkit, for nodes not presents

// (t,u,v) triplet, u v in same order as csv
struct Interaction {
    Time t;
    node u, v;
    
    Interaction() : t(none), u(none), v(none) {}

    Interaction(Time _t, node _u, node _v) {
        t = _t;
        //u = std::min(_u, _v);
        //v = std::max(_u, _v);
        u = _u;
        v = _v;
    }
};

// edge (u,v) with u <= v (order useful for counter map)
struct Edge {
    node u, v;

    Edge() : u(none), v(none) {}

    Edge(node _u, node _v) {
        u = std::min(_u, _v);
        v = std::max(_u, _v);
    }
    bool operator <( const Edge &rhs ) const
    {
       if (u < rhs.u) {
            return true;
       } else if (u == rhs.u && v < rhs.v) {
           return true;
       } else {
           return false;
       }
    }
    bool operator ==( const Edge &rhs ) const
    {
        if (u == rhs.u && v == rhs.v ) {
            return true;
        } else {
            return false;
        }
    }
};
inline void us_isect(std::unordered_set<node> &out,
        const std::unordered_set<node> &in1,
        const std::unordered_set<node> &in2)
{   // get intersection of two unordered_set 
    out.clear();
    if (in2.size() < in1.size()) {
        us_isect(out, in2, in1);
        return;
    }
    for (std::unordered_set<node>::const_iterator it = in1.begin(); it != in1.end(); it++)
    {
        if (in2.find(*it) != in2.end())
            out.insert(*it);
    }
};

}

namespace std{

    template <>
    struct hash<std::pair<uint64_t, uint64_t>>
    {
        std::size_t operator()(const std::pair<node, node>& k) const
        {
            using std::size_t;
            using std::hash;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return hash<uint64_t>()(k.first) ^ hash<uint64_t>()(k.second) ;
        }
    };
    template <>
    struct hash<StreamGraphs::Edge>
    {
        std::size_t operator()(const StreamGraphs::Edge& k) const
        {
            using std::size_t;
            using std::hash;

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return hash<uint64_t>()(k.u) ^ hash<uint64_t>()(k.v) ;
        }
    };

}

#endif
