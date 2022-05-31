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

    void  HashGraph::addEdge(node u, node v, int weight){
        const Edge edge = Edge(u,v);
        const Edge edge2 = Edge(v,u);
        // verifier le node u exist ou pas
        if(HashGraph::Neighbors.find(u)==HashGraph::Neighbors.end()){
            std::vector<node> listNeighbour;
            HashGraph::Neighbors.emplace(u,listNeighbour);
        }
        // verifier le node v exist ou pas
        if(HashGraph::Neighbors.find(v)==HashGraph::Neighbors.end()){
            std::vector<node> listNeighbour;
            HashGraph::Neighbors.emplace(v,listNeighbour);
        }
        //verifier edge exist ou pas
        if(HashGraph::Edges.find(edge)==HashGraph::Edges.end()){
            int indice_u = HashGraph::Neighbors.at(v).size();
            int indice_v  = HashGraph::Neighbors.at(u).size();
            HashGraph::Neighbors.at(u).emplace_back(v);
            HashGraph::Neighbors.at(v).emplace_back(u);
            HashGraph::Edges.emplace(edge,indice_u);
            HashGraph::Edges.emplace(edge2,indice_v);
        }
        else{
            HashGraph::Edges_weights[edge] += weight;
            HashGraph::Edges_weights[edge2] += weight;
        }

        /*
        int indice_u  = 0;
        for(int i= 0; i<HashGraph::Neighbors.at(u).size();i++){
            if(HashGraph::Neighbors.at(u).at(i)==v){
                indice_u = i;
                break;
            }
        }

        if(indice_u == HashGraph::Neighbors.at(u).size()) HashGraph::Neighbors.at(u).emplace_back(v);

        int indice_v  = 0;
        for(int i= 0; i<HashGraph::Neighbors.at(v).size();i++){
            if(HashGraph::Neighbors.at(v).at(i)==u){
                indice_v = i;
                break;
            }
        }
        if(indice_v == HashGraph::Neighbors.at(v).size()) HashGraph::Neighbors.at(v).emplace_back(u);
        */
        
        //ajouter sur weightEdge
        if(HashGraph::Edges_weights.find(edge)==HashGraph::Edges_weights.end()){
            HashGraph::Edges_weights.emplace(edge,weight);
        }
        if(HashGraph::Edges_weights.find(edge2)==HashGraph::Edges_weights.end()){
            HashGraph::Edges_weights.emplace(edge2,weight);
        }

        //HashGraph::Edges_weights.at(edge) = weight;
    }

    // number of nodes
    Count HashGraph::degree(node u) const{
        if(HashGraph::Neighbors.find(u)!=HashGraph::Neighbors.end()){
            return HashGraph::Neighbors.at(u).size();
        }
        else {
            return  0;
        }
    }

    void HashGraph::increaseWeight(node u, node v, int weight){
        HashGraph::addEdge(u, v, 1);
        const Edge edge = Edge(u,v);
        //neibhours verifice
        //HashGraph::Edges_weights[edge] += weight;

    }

    node HashGraph::addNode(){
        node res = node(z);
        N_nodes++;
        z++;
        exists.emplace_back(true);
        std::vector<node> listNeighbour;
        HashGraph::Neighbors.emplace(res,listNeighbour);
        return res;
    }

    void HashGraph::restoreNode(node u){
        if(u < z){
            if(!exists[u]){
                exists[u] = true;
                N_nodes++;
                std::vector<node> listNeighbour;
                HashGraph::Neighbors.emplace(u,listNeighbour);
            }
        }
        //assert(u < z);
        //assert(!exists[u]);
    }

    void HashGraph::removeNode(node u){
        if(u<z && exists[u]){
            exists[u]=false;
            N_nodes--;
            HashGraph::Neighbors.erase(u);
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