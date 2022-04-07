#include <networkit/graph/Graph.hpp>
#include <networkit/Globals.hpp>
//#include <ext/pb_ds/assoc_container.hpp>
//#include <sparsehash/sparse_hash_map> // or sparse_hash_set, dense_hash_map, ...
#include "hash_graph.hpp"
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
#include <limits>
#include <utility>
//#include <string_view>

using namespace NetworKit;

namespace StreamGraphs {


//google::sparse_hash_set<int, int> table;
    //google::sparse_hash_map<node, std::vector<node>,std::hash<uint64_t>> Neighbors;
    //std::map<node, std::vector<node>,std::hash<uint64_t>> Neighbors;

    HashGraph::HashGraph(int N_nodes, bool is_weighted, bool is_directed):
    N_nodes(N_nodes),
    exists(N_nodes, true),
    is_weighted(is_weighted),
    z(N_nodes),
    is_directed(is_weighted)
    {}

    HashGraph::~HashGraph(){}

    // number of nodes
    Count HashGraph::degree(node u) const{
        return HashGraph::Neighbors.at(u).size();
    }

    void HashGraph::increaseWeight(node u, node v, int weight){
        const Edge edge = Edge(u,v);
        HashGraph::Edges_weights[edge] += weight;

    }

    node HashGraph::addNode(){
        node res = z;
        N_nodes++;
        z++;
        exists.emplace_back(true);
        return res;
    }

    void HashGraph::restoreNode(node u){
        if(u < z){
            if(!exists[u]){
                exists[u] = true;
                N_nodes++;
            }
        }
        //assert(u < z);
        //assert(!exists[u]);
    }

    void HashGraph::removeNode(node u){
        if(u<z && exists[u]){
            exists[u]=false;
            N_nodes--;
        }
    }

    void HashGraph::setWeight(node u, node v, int w){
        const Edge edge = Edge(u,v);
        HashGraph::Edges_weights[edge] = w;
    }

    bool HashGraph::hasEdge(node u, node v){
        const Edge edge = Edge(u,v);
        if(HashGraph::Edges.find(edge) != HashGraph::Edges.end()){
            return true;
        }
        else return false;
    }


    void HashGraph::removeEdge(node u, node v){
        const Edge edge = Edge(u,v);
        HashGraph::Edges.erase(edge);
        HashGraph::Edges_weights.erase(edge);
    }

    int HashGraph::numberOfNodes(){
        return N_nodes;
    }

    HashGraph::NeighborRange HashGraph::neighborRange(node u)const {
        assert(exists[u]);
        return  NeighborRange(*this, u);
    }
    Count HashGraph::weight(node u,node v){
        //Edge::Edge()
        const Edge edge = Edge(u,v);
        return Edges_weights.at(edge);
    }






}