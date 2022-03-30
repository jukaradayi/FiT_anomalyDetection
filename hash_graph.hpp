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
//#include <string_view>

using namespace NetworKit;
namespace StreamGraphs {


class HashGraph { // TODO check who is private or public


//google::sparse_hash_set<int, int> table;
public:

    //google::sparse_hash_map<node, std::vector<node>,std::hash<uint64_t>> Neighbors;
    std::map<node, std::vector<node>,std::hash<uint64_t>> Neighbors;
    // hash table of nodes
    //google::sparse_hash_map<std::pair<node,node>, int,std::hash<std::pair<uint64_t, uint64_t>>> Edges;
    std::map<std::pair<node,node>, int,std::hash<std::pair<uint64_t, uint64_t>>> Edges;


    //google::sparse_hash_set<Edge, std::pair<int,int>> Edges;
    //google::sparse_hash_map<Edge, int,std::hash<StreamGraphs::Edge>> Edges_weights;
    std::map<Edge, int,std::hash<Edge>> Edges_weights;
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
            }

            NeighborIterator end() const {
                assert(G);
                return NeighborIterator(G->Neighbors.at(u).end());
            }

        };

    HashGraph(int N_nodes, bool is_weighted, bool is_directed);

    ~HashGraph();


    Count degree(node u);

    void increaseWeight(node u, node v, int weight);

    node addNode();

    void restoreNode(node u);

    void removeNode(node u);

    void setWeight(node u, node v, int w);

    bool hasEdge(node u, node v);

    void removeEdge(node u, node v);

    int numberOfNodes();

    Count weight(node u,node v);

    //void neighborRange();

    bool hasNode(node v) const noexcept { return (v < z) && this->exists[v]; }


    HashGraph::NeighborRange neighborRange(node u);






private:
    bool is_weighted;
    bool is_directed;



};

}
#endif