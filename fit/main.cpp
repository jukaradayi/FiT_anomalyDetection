// history graphs
// author : Nicolas Gensollen, Julien Karadayi, Matthieu Latapy
// 
// Implementation of graphs bounded on the node degree

#include "history_graph.hpp"
#include "G_graph.hpp"
#include "H_graph.hpp"
#include "metrics.hpp"
#include "csv.hpp"

//#include <omp.h>
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


int main(int argc, char* argv[]) {
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
    // -c checks

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
        Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        node_set.insert(i.u);
        node_set.insert(i.v);
        
    }
    file.clear();
    file.seekg(0);
    file.close();

    // Handle different graph sizes
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
       
        // if graph G, get first line
        Time t0 = -1;
        bool compute_metrics = false;

        /***************
         ** Main Loop **
         ***************/
        for(CSVIterator main_loop = CSVRange(file).begin(); main_loop != CSVRange(file).end(); ++main_loop)
        {
            Interaction i(std::stoi((*main_loop)[0]), std::stoi((*main_loop)[1]), std::stoi((*main_loop)[2]));
            if (line_number == 0) {
                t0 = i.t;
            }

            if (i.u == i.v)  // skip self-loop, in darpa for example
            {
                ++line_number;
                continue;
            }

            // Update graph
            std::clock_t t1_update = std::clock(); // start clock 
            hist_graph->updateGraph(i);            // update graph
            std::clock_t t2_update = std::clock(); // end clock
            dur_update += (t2_update - t1_update) / (double)CLOCKS_PER_SEC; // time in seconds

            // start metrics at window size
            std::clock_t t1_metric = std::clock(); // start clock
            if ( (graph_type == "H" && line_number > hist_sizes[graph_idx]) 
              || (graph_type == "G" && (i.t - t0 > hist_sizes[graph_idx]))) {
               compute_metrics = true; 
            }


            // Compute metrics and output to file
            if (compute_metrics && (metrics1 || metrics2 || metrics3) )
            {
                output += metrics.run(i.u, i.v, line_number);
                myoutput << output << std::flush;
                output = "";
            }
            std::clock_t t2_metric = std::clock(); // end clock
            dur_metric += (t2_metric - t1_metric) / (double)CLOCKS_PER_SEC; // time in seconds
            dur_metric_whole += (t2_metric - t1_metric) / (double)CLOCKS_PER_SEC;


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

        // Get current time
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::string current_time = std::ctime(&end_time);
        current_time.erase(std::remove(current_time.begin(), current_time.end(), '\n'), current_time.end());

        mylog << current_time << " overall runtime " << std::to_string(cpu_duration) <<"\n" ;
        mylog.close();
    }
    return 0;
}
