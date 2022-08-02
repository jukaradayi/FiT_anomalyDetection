/*
 *  hash_graph.hpp
 *
 *  created on: 14.09.2021
 *  authors: Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
 *
 *
 */

#ifndef STREAMGRAPHS_HASH_GRAPH_HPP_
#define STREAMGRAPHS_HASH_GRAPH_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
//#include <ext/pb_ds/assoc_container.hpp>
// or sparse_hash_set, dense_hash_map, ...
//#include <sparsehash/sparse_hash_map>


#include "globals.hpp"
#include <limits>
//#include <boost>
#include <utility>
#include <cstdint>
#include <map>
#include <string>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <tsl/robin_map.h>
#include <typeindex>
//#include <string_view>

using namespace NetworKit;
namespace StreamGraphs {


class HashGraph { // TODO check who is private or public


//google::sparse_hash_set<int, int> table;
public:

    //google::sparse_hash_map<node, std::vector<node>,std::hash<uint64_t>> Neighbors;
    //std::unordered_map<node, std::vector<node>> Neighbors;//,std::hash<node>>  ;
    tsl::robin_map<node, std::vector<node>> Neighbors;
    // hash table of nodes
    //google::sparse_hash_map<std::pair<node,node>, int,std::hash<std::pair<uint64_t, uint64_t>>> Edges;

    //std::unordered_map<std::pair<node,node>, int,std::hash<std::pair<node, node>>> Edges;
    //std::unordered_map<Edge, int,std::hash<Edge>> Edges;
    tsl::robin_map<Edge, int,std::hash<Edge>> Edges;

    //google::sparse_hash_set<Edge, std::pair<int,int>> Edges;
    //google::sparse_hash_map<Edge, int,std::hash<StreamGraphs::Edge>> Edges_weights;
    //std::unordered_map<Edge, int,std::hash<Edge>> Edges_weights;
    tsl::robin_map<Edge, int,std::hash<Edge>> Edges_weights;

    // number of nodes
    Count N_nodes;
    //current upper bound of node ids, z will be the id of the next node
    node z;
    //NetworKit::Graph g;
    std::vector<bool> exists;

    class NeighborIterator {

        std::vector<node>::const_iterator nIter;

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
            using self = NeighborIterator;

            NeighborIterator(std::vector<node>::const_iterator nodesIter) : nIter(nodesIter) {}

            /**
             * @brief WARNING: This contructor is required for Python and should not be used as the
             * iterator is not initialized.
             */
            NeighborIterator() {}

            NeighborIterator &operator++() {
                ++nIter;
                return *this;
            }

            NeighborIterator operator++(int) {
                const auto tmp = *this;
                ++nIter;
                return tmp;
            }

            NeighborIterator operator--() {
                const auto tmp = *this;
                --nIter;
                return tmp;
            }

            NeighborIterator operator--(int) {
                --nIter;
                return *this;
            }

            bool operator==(const NeighborIterator &rhs) const { return nIter == rhs.nIter; }

            bool operator!=(const NeighborIterator &rhs) const { return !(nIter == rhs.nIter); }

            node operator*() const { return *nIter; }
        };

    class NeighborRange {
            const HashGraph *G;
            node u;

        public:
            NeighborRange(const HashGraph &G, node u) : G(&G), u(u) { assert(G.hasNode(u)); };

            NeighborRange() : G(nullptr){};

            NeighborIterator begin() const {
                assert(G);
                //imEdges = Edges.get(u)
                //G->Edges.find(u)
                return  NeighborIterator(G->Neighbors.at(u).begin());
                //std::vector<node>::const_iterator nodesIter;
                //return nodesIter;
            }

            NeighborIterator end() const {
                assert(G);
                //return NeighborIterator(G->Neighbors.at(u).end());
                std::vector<node>::const_iterator nodesIter;
                return nodesIter;
            }

        };

    /*****************************/
    /* Bounded neighbor iterator */
    /*****************************/
    class BoundedNeighborIterator {
            const HashGraph *G;
            const Bound bound = 0;


            //std::vector<node>::const_iterator nIter;
            //NetworKit::Graph::NeighborIterator nIter;
            //NetworKit::Graph::NeighborIterator end;

            HashGraph::NeighborIterator nIter;
            HashGraph::NeighborIterator end;

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
            //BoundedNeighborIterator(const NetworKit::Graph &G, NetworKit::Graph::NeighborIterator neighborsIter, Bound bound, NetworKit::Graph::NeighborIterator end) : G(&G), nIter(neighborsIter), bound(bound), end(end) {}
            BoundedNeighborIterator(const HashGraph &G, HashGraph::NeighborIterator neighborsIter, Bound bound, HashGraph::NeighborIterator end) : G(&G), nIter(neighborsIter), bound(bound), end(end) {}



            /**
             * @brief WARNING: This contructor is required for Python and should not be used as the
             * iterator is not initialized.
             */
            BoundedNeighborIterator() {}

            BoundedNeighborIterator operator++() {
                //const auto tmp = *this;

                if (nIter == end) {
                    return *this;
                }

                do {
                    ++nIter;
                } while ( nIter != end && G->degree(*nIter) > bound);
                return *this;
                //// check node iterator from NK
                ////while (G->hasNode(*++nIter) && G->degree(*++nIter) > bound) {}
                //std::cout << "before while "<< *nIter << "\n";
                //++nIter;
                //while (G->degree(*nIter) > bound && std::next(nIter, 1) != end) {
                //    std::cout << "toto";
                //    ++nIter;
                //}
                //return tmp;
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

            node operator*() const { return *nIter; }
    };

    class BoundedNeighborRange {
            //const NetworKit::Graph *G;
            const HashGraph *G;
            const Bound bound = 0;
            node u;

        public:
            BoundedNeighborRange(const HashGraph &G, node u, Bound bound) : G(&G), u(u), bound(bound) { assert(G.hasNode(u)); };

            BoundedNeighborRange() : G(nullptr){};

            BoundedNeighborIterator begin() const {
                assert(G);
                //NetworKit::Graph::NeighborRange<false>  neighborRange = G->NeighborRange(u);
                //std::vector<node> uNeighbors(G->neighborRange(u).begin(), G->neighborRange(u).end());
                return BoundedNeighborIterator(*G, G->neighborRange(u).begin(), bound, G->neighborRange(u).end());
                //return BoundedNeighborIterator(uNeighbors.begin(), bound);

            }

            BoundedNeighborIterator end() const {
                assert(G);
                //NetworKit::Graph::NeighborRange<false>  neighborRange = G->NeighborRange(u);

                return BoundedNeighborIterator(*G, G->neighborRange(u).end(), bound, G->neighborRange(u).end());
            }
    };


    HashGraph(int N_nodes, bool is_weighted, bool is_directed);

    ~HashGraph();

// comparer le temps
    Count degree(const node u) const;

    Count degree(node u);

    Count weightedDegree(node u);

    void increaseWeight(node u, node v, int weight);

    node addNode();

    void addEdge(node u, node v, int weight);

    void restoreNode(node u);

    void removeNode(node u);

    void setWeight(node u, node v, int w);

    bool hasEdge(node u, node v);

    void removeEdge(node u, node v);

    index upperNodeIdBound();

    //end of comparer

    int numberOfNodes();

    Count weight(node u,node v);

    //void neighborRange();

    bool hasNode(node v) const noexcept { return (v < z) && this->exists[v]; }


    HashGraph::NeighborRange neighborRange(node u)  const;






private:
    bool is_weighted;
    bool is_directed;



};

}
#endif
