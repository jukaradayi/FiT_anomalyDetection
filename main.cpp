// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "G_graph.hpp"
#include "H_graph.hpp"
#include "metrics.hpp"

#include <networkit/graph/Graph.hpp>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class CSVRow
{
    public:
        std::string operator[](std::size_t index) const
        {
            return std::string(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
        }
        std::size_t size() const
        {
            return m_data.size() - 1;
        }
        void readNextRow(std::istream& str)
        {
            std::getline(str, m_line);

            m_data.clear();
            m_data.emplace_back(-1);
            std::string::size_type pos = 0;
            while((pos = m_line.find(' ', pos)) != std::string::npos)
            {
                m_data.emplace_back(pos);
                ++pos;
            }
            // This checks for a trailing comma with no data after it.
            pos   = m_line.size();
            m_data.emplace_back(pos);
        }
    private:
        std::string         m_line;
        std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   

class CSVIterator
{   
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) {++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

class CSVRange
{
    std::istream&   stream;
    public:
        CSVRange(std::istream& str)
            : stream(str)
        {}
        CSVIterator begin() const {return CSVIterator{stream};}
        CSVIterator end()   const {return CSVIterator{};}
};

int main() {
    using namespace StreamGraphs;
    std::ifstream file("peru_10000.csv");
    //std::ifstream file("exemple_minimal.csv");

    //CSVRow row;
    //while(file >> row)
    //{
    //    std::cout << "1th Element(" << row[0] << ")\n";
    //}
    //  
    std::unordered_set<node> node_set;
   
    // first loop, get node set  
    //for(CSVIterator node_loop(file); node_loop != CSVIterator(); ++node_loop)
    for(auto& node_loop: CSVRange(file))
    {

        //time t = *loop[0];
        Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        node_set.insert(i.u);
        node_set.insert(i.v);
        
    }
    GGraph g_graph(true, false, 10000, 10000, node_set.size(), 3600);
    //HGraph h_graph(true, false, 1000, 1000, node_set.size(), 1000);
    Metrics g_metrics(g_graph, true, true, true, true);
    //Metrics h_metrics(h_graph,true, false);
    //Motric g_metrics(g_graph, true);
    //g_metrics.add(g_metrics.number_of_nodes);
    //g_graph.nctions.push_back(g_graph.*trimQueue);
    //g_graph.nctions.push_back(toto);

    //GGraph g_graph(false, false, 1000, 1000, N, 1000);
    //std::ifstream file("peru_10000.csv");
    file.clear();
    file.seekg(0);



    // create graph
    //for(auto& main_loop: CSVRange(file))
    int idx = 0;
    double cap_shrink_to_fit = std::pow(10.0, 6);
    for(CSVIterator main_loop = CSVRange(file).begin(); main_loop != CSVRange(file).end(); ++main_loop)
    {
        ++idx;
        Interaction i(std::stoi((*main_loop)[0]), std::stoi((*main_loop)[1]), std::stoi((*main_loop)[2]));
        g_graph.updateGraph(i); 
        //for (const auto edge: NetworKit::Graph::EdgeRange(g_graph.top_graph)) {
        //    std::cout << "top graph has edge "<< g_graph.top2node[edge.u] << " " << g_graph.top2node[edge.v] << "\n";
        //}
        //std::cout << "[";
        //int idx = 0;
        //for (const auto edge: NetworKit::Graph::EdgeRange(g_graph.bot_graph)) {
        //    if (idx > 0) {
        //        std::cout <<"), ";
        //    }
        //    ++idx;
        //    std::cout << "(" << g_graph.bot2node[edge.u] << ", " << g_graph.bot2node[edge.v];
        //}
        //std::cout << ")]\n";
        //h_graph.updateGraph(i);
        if (idx > 499){
            std::cout << idx << "\n";
            g_metrics.run(i.u, i.v);
        }
        if (idx % 2 * cap_shrink_to_fit == 0) {
            g_graph.main_graph.shrinkToFit();
            g_graph.top_graph.shrinkToFit();
            g_graph.bot_graph.shrinkToFit();
        }
        //h_metrics.run(i.u, i.v);
    }




   return 0;
}

