//
// Created by Apple on 30/05/2022.
//

#include "sorted_degree.hpp"
//#include "main.cpp"
namespace StreamGraphs {
    ConstantTimeMax::ConstantTimeMax(std::unordered_set<uint64_t> Nodes) :
    Nodes(Nodes)
    {
       // this->deg_seq = new uint64_t[N] {none};
       //need initialise taille?
    }
    ConstantTimeMax::ConstantTimeMax(){
        //this->Nodes = new std::unordered_set<uint64_t> {None};

    }
    ConstantTimeMax::~ConstantTimeMax(){}
    void ConstantTimeMax::addNode(node u) {
        Nodes.emplace(u);
        //deg_seq.emplace_back(1);
        deg_seq.push_back(1);
        node2pos[u]  = deg_seq.size()-1;
        deg_set.insert(1);

        if(deg2pos.size()>0){
            deg2pos[0]+=1;
        }
        else{
            deg2pos.emplace_back(1);
            deg2pos.emplace_back(0);
        }

        if(pos2node.size()>deg_seq.size()-1){
            pos2node[deg_seq.size()-1]=u;
        }
        else{
            //pos2node.emplace_back(u);
            pos2node.push_back(u);
        }
        if(deg_distrib.size()>0){
            deg_distrib[0]+=1;
        } else{
            deg_distrib.emplace_back(1);
        }
    }
    void ConstantTimeMax::increaseNodeDeg(node u){
        int u_pos = node2pos[u];
        int u_deg = deg_seq[u_pos];
        //int u_deg = deg_seq

        int v_pos = deg2pos[u_deg];
        int v = pos2node[v_pos];

        pos2node[v_pos]=u;
        pos2node[u_pos]=v;
        node2pos[u] = v_pos;
        node2pos[v] = u_pos;

        deg_seq[v_pos]+=1;

        deg_distrib[u_deg-1]-=1;
        if(deg_distrib.size()>u_deg){
            deg_distrib[u_deg-1+1]+=1;
        }else{
            deg_distrib.emplace_back(1);
        }

        deg2pos[u_deg]+=1;
        if(u_deg == deg2pos.size()-1){
            deg2pos.emplace_back(0);
        }
        if(deg2pos[u_deg+1] == -1){
            deg2pos[u_deg+1] = 1;
        }
    }
    void ConstantTimeMax::decreaseNodeDeg(node u){
        int u_pos = node2pos[u];
        int u_deg = deg_seq[u_pos];
        //std::cout<<"===== avant la change ===="<<std::endl;
        int v_pos = deg2pos[u_deg-1]-1;
        node v = pos2node[v_pos];
        pos2node[v_pos] = u;
        pos2node[u_pos] = v;
        node2pos[u] = v_pos;
        node2pos[v] = u_pos;

        //std::cout<<"===== apres la change ===="<<std::endl;
        deg_seq[v_pos] -= 1;
        if(deg_seq[v_pos]==0){
            //std::pop_heap(deg_seq.begin(), deg_seq.end())
            deg_seq.pop_back();//？pop end?
            node2pos.erase(u);
            pos2node.pop_back();
            Nodes.erase(u);
        }
        deg_distrib[u_deg-1] -=1;
        if(u_deg-1 > 0){
            deg_distrib[u_deg-1-1]+=1;
        }
        //std::cout<<"===== u_degree ===="<<u_deg<<std::endl;
        if(deg_distrib[u_deg-1]==0 && deg2pos[u_deg]>0){
            deg2pos[u_deg-1] -=1;// 4,3,2 supprimer 3 , 2 il faut avancer
        }else if(deg_distrib[u_deg-1]==0 && deg2pos[u_deg]==0){
            //deg2pos.pop_front();
            deg2pos.pop_back();
            deg2pos[u_deg-1] = 0;
        }else{
            deg2pos[u_deg-1] -=1;
        }

    }
    std::string ConstantTimeMax::affichage(){
        std::string S ="--- degree sequence ---\n";
        for(int i=0;i<deg_seq.size();i++){
            //std::cout<<deg_seq.at(i)<<std::endl;
            //boost::lexical_cast<std::string>(val)
            S+=std::to_string(deg_seq.at(i))+",";
        }
        S+="\n";
        S+="--- node positions ---\n";
        for(auto & tmp:node2pos){
            S+=std::to_string(tmp.first)+" : "+std::to_string(tmp.second)+", ";
        }
        S+="\n";
        S+="--- position to nodes ---\n";
        for(int i=0;i<pos2node.size();i++){
            S+=std::to_string(pos2node.at(i))+",";
        }
        S+="\n";
        S+="--- degree to position ---\n";
        for(int i=0;i<deg2pos.size();i++){
            S+=std::to_string(deg2pos.at(i))+",";
        }
        S+="\n";
        S+="--- degree distribution ---\n";
        for(int i=0;i<deg_distrib.size();i++){
            S+=std::to_string(deg_distrib.at(i))+",";
        }
        S+="\n";
        return S;
    }




}