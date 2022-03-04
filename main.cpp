// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "G_graph.hpp"
#include "H_graph.hpp"
#include "metrics.hpp"

#include <omp.h>
#include <networkit/graph/Graph.hpp>
#include <iterator>
#include <chrono>
#include <time.h>
#include <ctime>
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
class InputParser{

    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        const std::pair<std::vector<StreamGraphs::Count>, StreamGraphs::Count> getSizes() const{
            std::vector<StreamGraphs::Count> hist_sizes;
            std::vector<std::string>::const_iterator itrH =  std::find(this->tokens.begin(), this->tokens.end(), "-h");
            std::vector<std::string>::const_iterator itrG =  std::find(itrH, this->tokens.end(), "-g");
            std::vector<std::string>::const_iterator next_cmd = std::find(itrG+1, this->tokens.end(), "-k");

            StreamGraphs::Count G_index = 0;
            if (itrH != this->tokens.end()){
                for (++itrH; itrH != itrG; ++itrH) {
                    ++G_index;
                    const std::string& h_size = *itrH;
                    std::cout << h_size << "\n";

                    hist_sizes.push_back(std::stoi(h_size));
                }
            }

            if (itrG != this->tokens.end()){
                for (++itrG; itrG != next_cmd; ++itrG) {
    
                    const std::string& g_size = *itrG;
                    std::cout << g_size << "\n";

                    hist_sizes.push_back(std::stoi(g_size));
                }
            }
            std::pair<std::vector<StreamGraphs::Count>, StreamGraphs::Count> Sizes = std::make_pair(hist_sizes, G_index);
            return Sizes;
         
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

int main(int argc, char* argv[]) {
    //auto writeLog = [&](std::string output) {
    //            auto end = std::chrono::system_clock::now();
    //            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
    //            std::string current_time = std::ctime(&end_time);
    //            current_time.erase(std::remove(current_time.begin(), current_time.end(), '\n'), current_time.end());

    //            //mylog.open (logname, std::ios::app);
    //            mylog << current_time << " interaction " << line_number<<"\n" ;
    //            mylog << current_time << " average update time " << std::to_string(dur_update/100) <<"\n" ;
    //            mylog << current_time << " average metric time " << std::to_string(dur_metric/100) <<"\n" ;
    //            mylog << current_time << " overall average metric time " << std::to_string(dur_metric_whole/(line_number-499)) <<"\n" << std::flush;

    //            //mylog.close();
    //            dur_metric = 0;
    //            dur_update = 0;
    //}
    using namespace StreamGraphs;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    InputParser input(argc, argv);

    // arguments: 
    // -f filename
    // -o output folder
    // -h H window sizes
    // -g G window sizes
    // -k limit & proj limit
    // -m1 o(1)
    // -m2 local
    // -m3 global
    // -p projeté
    // -b bipartite

    // get parameters
    std::pair<std::vector<Count>, Count> Sizes = input.getSizes();
    Count G_index = Sizes.second;
    std::vector<Count> hist_sizes = Sizes.first;

    const std::string &filename = input.getCmdOption("-f");
    const std::string &output_folder = input.getCmdOption("-o");
    const Count limit = std::stoi(input.getCmdOption("-k"));
    const bool metrics1 = input.cmdOptionExists("-m1");
    const bool metrics2 = input.cmdOptionExists("-m2");
    const bool metrics3 = input.cmdOptionExists("-m3");
    //const bool metrics4 = input.cmdOptionExists("-m4");
    const bool use_proj = input.cmdOptionExists("-p");
    const bool is_bip = input.cmdOptionExists("-b");

    // check parameters
    if (filename.empty() || output_folder.empty()){
        std::cerr << "input file not provided\n";
        std::cerr << "format is: " << argv[0] << " -f [FILENAME] -o [OUTPUT FOLDER] -h [H GRAPH SIZES] -g [G GRAPH SIZES] -k [DEGREE BOUND]\n";
        return  1;
    }

    std::ifstream file(filename);
    std::unordered_set<node> node_set;
   
    // first loop, get node set // TODO REMOVE THIS, LOOP OVER THE FILE ONLY ONCE 
    for(auto& node_loop: CSVRange(file))
    {
        //time t = *loop[0];
        Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        node_set.insert(i.u);
        node_set.insert(i.v);
        
    }
    file.clear();
    file.seekg(0);
    file.close();

    // parallelize over window sizes // TODO check slower parallelization...
    omp_set_dynamic(0);
    #pragma omp parallel for num_threads(10)
    for (int graph_idx = 0; graph_idx < hist_sizes.size(); ++graph_idx) {
        std::string graph_type = (graph_idx < G_index) ? "H" : "G";

        // init NK graphs
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph;
        if (is_bip) {
            bot_graph = NetworKit::Graph(0, true, false);
        } else {
            bot_graph = top_graph;
        }

        // Init history graphs
        HistoryGraph* hist_graph = NULL; // TODO FACTORY PATTERN
        if (graph_type == "H") {
            hist_graph = new HGraph(main_graph, top_graph, bot_graph, use_proj, false, is_bip, limit, limit, node_set.size(), hist_sizes[graph_idx]);
        } else if (graph_type == "G") {
            hist_graph = new GGraph(main_graph, top_graph, bot_graph, use_proj, false, is_bip, limit, limit, node_set.size(), hist_sizes[graph_idx]);
        }

        Metrics metrics(*hist_graph, metrics1, metrics2, metrics3, false);

        std::ifstream file(filename);

        int line_number = 0;

        // metric outputs
        std::string output = "";
        std::ofstream myoutput;
        std::string outname = output_folder + "/Output_"+ graph_type +"_" + std::to_string(hist_sizes[graph_idx]) + ".csv";
        double dur_metric;
        double dur_metric_whole;
        double dur_update;

        // log outputs
        std::string log_output = "";
        std::ofstream mylog;
        std::string logname = output_folder + "/Output_"+ graph_type +"_" + std::to_string(hist_sizes[graph_idx]) + ".log";

        std::clock_t startcputime = std::clock();
        myoutput.open (outname, std::ios::app);
        mylog.open (logname, std::ios::app);
        
        for(CSVIterator main_loop = CSVRange(file).begin(); main_loop != CSVRange(file).end(); ++main_loop)
        {
            Interaction i(std::stoi((*main_loop)[0]), std::stoi((*main_loop)[1]), std::stoi((*main_loop)[2]));
            //output += std::to_string(line_number) + ","; // get interaction id 

            if (i.u == i.v)  // skip self-loop, in darpa for example
            {
                ++line_number;
                continue;
            }

            std::clock_t t1_update = std::clock(); // start clock 
            hist_graph->updateGraph(i);            // update graph
            std::clock_t t2_update = std::clock(); // end clock
            dur_update += (t2_update - t1_update) / (double)CLOCKS_PER_SEC; // time in seconds

            // start metrics at line 500 // TODO start metrics at window size
            std::clock_t t1_metric = std::clock(); // start clock
            if (line_number > 499 && (metrics1 || metrics2 || metrics3) )
            {
                output += metrics.run(i.u, i.v, line_number);
                myoutput << output << std::flush;
                output = "";
            }
            std::clock_t t2_metric = std::clock(); // end clock
            dur_metric += (t2_metric - t1_metric) / (double)CLOCKS_PER_SEC; // time in seconds
            dur_metric_whole += (t2_metric - t1_metric) / (double)CLOCKS_PER_SEC;

            //}
            // TODO tester sans shrink to fit
            //if (line_number % 100000== 0) {
            //    std::cout << "shrinking to fit " << line_number << std::endl;
            //    hist_graph->main_graph.shrinkToFit();
            //    if (use_proj){
            //        hist_graph->top_graph.shrinkToFit();
            //        hist_graph->bot_graph.shrinkToFit();
            //    }
            //}
            if (line_number % 100 == 0) { // write log every 100 interaction
                // Get current time
                auto end = std::chrono::system_clock::now();
                std::time_t end_time = std::chrono::system_clock::to_time_t(end);
                std::string current_time = std::ctime(&end_time);
                current_time.erase(std::remove(current_time.begin(), current_time.end(), '\n'), current_time.end());

                //mylog.open (logname, std::ios::app);
                mylog << current_time << " interaction " << line_number<<"\n" ;
                mylog << current_time << " average update time " << std::to_string(dur_update/100) <<"\n" ;
                mylog << current_time << " average metric time " << std::to_string(dur_metric/100) <<"\n" ;
                mylog << current_time << " overall average metric time " << std::to_string(dur_metric_whole/(line_number-499)) <<"\n" << std::flush;

                //mylog.close();
                dur_metric = 0;
                dur_update = 0;
            }   

            ++line_number;
        }
        myoutput.close();
        mylog.close();

        std::clock_t endcputime = std::clock();
        double cpu_duration = (endcputime - startcputime) / (double)CLOCKS_PER_SEC;
        mylog.open (logname, std::ios::app);
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::string current_time = std::ctime(&end_time);
        current_time.erase(std::remove(current_time.begin(), current_time.end(), '\n'), current_time.end());

        mylog << current_time << " overall runtime " << std::to_string(cpu_duration) <<"\n" ;
        mylog.close();
    }
    return 0;
}
//TODO look at https://stackoverflow.com/questions/3201538/how-to-read-a-gz-file-line-by-line-in-c/3201675 to read gzip !
