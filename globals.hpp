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

/*****************************/
/* Bounded neighbor iterator */
/*****************************/
class BoundedNeighborIterator { // deprecated, use forBoundedNodes 
    const NetworKit::Graph *G;
    const Bound bound = 0;


    //std::vector<node>::const_iterator nIter;
    NetworKit::Graph::NeighborIterator nIter;
    NetworKit::Graph::NeighborIterator end;
    node u;

public:
    // The value type of the neighbors (i.e. nodes). Returned by
    // operator*().
    using value_type = node;

    // Reference to the value_type, required by STL.
    using reference = value_type &;

    // Pointer to the value_type, required by STL.
    using pointer = value_type *;

    // STL iterator category.
    using iterator_category = std::forward_iterator_tag;

    // Signed integer type of the result of subtracting two pointers,
    // required by STL.
    using difference_type = ptrdiff_t;

    // Own type.
    using self = BoundedNeighborIterator;

    //BoundedNeighborIterator(std::vector<node>::const_iterator nodesIter, Bound bound) : nIter(nodesIter), bound(bound) {}
    BoundedNeighborIterator(const NetworKit::Graph &G, NetworKit::Graph::NeighborIterator neighborsIter, Bound bound, NetworKit::Graph::NeighborIterator end, node u) : G(&G), nIter(neighborsIter), bound(bound), end(end) , u(u){}

    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    BoundedNeighborIterator() {}

    BoundedNeighborIterator operator++() {
        //const auto tmp = *this;
        if (nIter == end) {
            //std::cout << "returning end \n";
            return *this;
        }
        
        do {
            ++nIter;
        } while ( nIter != end && G->degree(*nIter) > bound);

        return *this;
    }

    BoundedNeighborIterator operator++(int) {

        const auto tmp = *this;
        if (nIter == end) {
            return *this;
        }
        
        do {
            ++nIter;
        } while ( nIter != end && G->degree(*nIter) > bound);
        return tmp;


    }

    //BoundedNeighborIterator operator--() {
    //    const auto tmp = *this;
    //    while (G->hasNode(*--nIter) && G->degree(*--nIter) > bound) {}
    //    //--nIter;
    //    return tmp;
    //}

    //BoundedNeighborIterator operator--(int) {
    //    while (G->hasNode(*--nIter) && G->degree(*--nIter) > bound) {}
    //    //--nIter;
    //    return *this;
    //}

    bool operator==(const BoundedNeighborIterator &rhs) const { return nIter == rhs.nIter; }

    bool operator!=(const BoundedNeighborIterator &rhs) const { return !(nIter == rhs.nIter); }

    //bool operator!=(const NetworKit::Graph::NeighborIterator &rhs) const { return !(nIter == rhs.nIter); }

    node operator*() const { return *nIter; }
};

class BoundedNeighborRange {
    const NetworKit::Graph *G;
    const Bound bound = 0;
    node u;

public:
    BoundedNeighborRange(const Graph &G, node u, Bound bound) : G(&G), u(u), bound(bound) { assert(G.hasNode(u)); };

    BoundedNeighborRange() : G(nullptr){};

    BoundedNeighborIterator begin() const {
        assert(G);
        //NetworKit::Graph::NeighborRange<false>  neighborRange = G->NeighborRange(u);
        //std::vector<node> uNeighbors(G->neighborRange(u).begin(), G->neighborRange(u).end());
        return BoundedNeighborIterator(*G, G->neighborRange(u).begin(), bound, G->neighborRange(u).end(), u);
        //return BoundedNeighborIterator(uNeighbors.begin(), bound);

    }

    BoundedNeighborIterator end() const {
    //NetworKit::Graph::NeighborIterator end() const {

        assert(G);
        //NetworKit::Graph::NeighborRange<false>  neighborRange = G->NeighborRange(u);

        return BoundedNeighborIterator(*G, G->neighborRange(u).end(), bound, G->neighborRange(u).end(), u);

    }
};





}

#endif
