#define CONFIG_CATCH_MAIN
//#include "catch.hpp"
#include <catch2/catch_test_macros.hpp>
#include "history_graph.hpp"
#include "H_graph.hpp"
#include "G_graph.hpp"
#include "globals.hpp"
#include "csv.hpp"
#include "metrics.hpp"
TEST_CASE("simple graph", "[StreamGraphs::HistoryGraph]") {
    using namespace StreamGraphs;
    std::ifstream file("unit_clean"); /// using hand made dataset on which results were computed by hand to check

    SECTION("Complete run through data") {
        //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        HGraph* H5 = new HGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 5); //TODO correct size

        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));
            H5->updateGraph(i); // update graph
        }
    }

    SECTION("Check H queue size") {
        //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        HGraph* H5 = new HGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 5); //TODO correct node set size
        HGraph* H10 = new HGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 10); //TODO correct node set size

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));
            line_idx += 1;
            H5->updateGraph(i); // update graph
            H10->updateGraph(i); // update graph

            if (line_idx >= 5) {
                REQUIRE( H5->queue.size() == 5);
            }
            if (line_idx >= 10) {
                REQUIRE( H10->queue.size() == 10);
            }

        }
    }

    SECTION("Check G node degrees") {
        //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        GGraph* G31 = new GGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 31);

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));

            line_idx += 1;
            G31->updateGraph(i); // update graph
            if (line_idx == 24) {

                node u_7 = G31->node2main[7];
                node u_24 = G31->node2main[24];
                REQUIRE(G31->main_graph.degree(u_7) == 7);
                REQUIRE(G31->main_graph.degree(u_24) == 1);
            }

        }
    }

    SECTION("Check main_graph degree") {
        //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        HGraph* H10 = new HGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 10); //TODO correct node set size

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));

            line_idx += 1;
            H10->updateGraph(i); // update graph
            if (line_idx == 24) {
                node u_7 = H10->node2main[7];
                node u_24 = H10->node2main[24];

                REQUIRE(H10->main_graph.degree(u_7) == 7);
                REQUIRE(H10->main_graph.degree(u_24) == 1);
            }
        }
    }

    SECTION("Check projections" ){
        //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        HGraph* H5 = new HGraph(main_graph, top_graph, bot_graph, true, false, true, 1000, 1000, 100, 5); //TODO correct node set size

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));

            line_idx += 1;
            if (line_idx == 11) {
                node u_8 = H5->node2top[8];
                node u_10 = H5->node2top[10];
                node u_6 = H5->node2top[6];
                node u_11 = H5->node2top[11];
                REQUIRE(H5->top_graph.hasEdge(u_8, u_10));
                REQUIRE(H5->top_graph.hasEdge(u_6, u_11));
            }

            H5->updateGraph(i); // update graph
            if (line_idx == 11) {
                node u_8 = H5->node2top[8];
                node u_10 = H5->node2top[10];
                node u_6 = H5->node2top[6];
                node u_11 = H5->node2top[11];
                node u_13 = H5->node2top[13];
                REQUIRE(H5->top_graph.hasEdge(u_8, u_10));
                REQUIRE(!H5->top_graph.hasEdge(u_6, u_11));
                REQUIRE(H5->top_graph.hasEdge(u_8, u_13));
                REQUIRE(H5->top_graph.hasEdge(u_10, u_13));
            }
        }
    }

    SECTION("Check Sorted Counters"){
    
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);

        HGraph* H5 = new HGraph(main_graph, top_graph, bot_graph, true, false, true, 1000, 1000, 100, 5); //TODO correct node set size

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));
            H5->updateGraph(i); // update graph

            for (auto main_node=H5->main_graph.nodeRange().begin(); main_node!=H5->main_graph.nodeRange().end(); ++main_node) {
                node mynode = *main_node;
                REQUIRE(H5->main_graph.degree(mynode) == H5->degree_counter.get_value(mynode));
                REQUIRE(H5->main_graph.weightedDegree(mynode) == H5->weightedDegree_counter.get_value(mynode));
            }
            for (auto main_edge=H5->main_graph.edgeRange().begin(); main_edge!=H5->main_graph.edgeRange().end(); ++main_edge) {
                StreamGraphs::Edge myedge((*main_edge).u, (*main_edge).v);
                REQUIRE(H5->main_graph.weight((*main_edge).u, (*main_edge).v) == H5->weight_counter.get_value(myedge));
            }

        }       
    }

    // Check that removal of link works correctly in metrics for m3 
    SECTION("Check link removal in metrics"){
    
        NetworKit::Graph main_graph(0, true, false);
        NetworKit::Graph top_graph(0, true, false);
        NetworKit::Graph bot_graph(0, true, false);
        std::string output = "";

        HGraph* H10 = new HGraph(main_graph, top_graph, bot_graph, true, false, true, 1000, 1000, 100, 10); //TODO correct node set size

        Metrics metrics(*H10, false, false, true, false);

        int line_idx = 0;
        for(auto& main_loop: CSVRange(file))
        {
            StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));
            H10->updateGraph(i); // update graph
            for(std::map<StreamGraphs::Edge, uint64_t>::iterator counterIt = H10->counter.begin(); counterIt != H10->counter.end(); ++counterIt)
            {
                StreamGraphs::Edge e =  counterIt->first;
                uint64_t weight = counterIt->second;
                StreamGraphs::Edge main_e(H10->node2main[e.u], H10->node2main[e.v]);

                // check that edge is in graph
                REQUIRE(H10->main_graph.hasEdge(main_e.u, main_e.v));
                REQUIRE(H10->main_graph.weight(main_e.u, main_e.v) == weight);


            }

            if (line_idx >= 10){
                output += metrics.run(i.u, i.v, line_idx);
            }

            for(std::map<StreamGraphs::Edge, uint64_t>::iterator counterIt = H10->counter.begin(); counterIt != H10->counter.end(); ++counterIt)
            {
                StreamGraphs::Edge e =  counterIt->first;
                uint64_t weight = counterIt->second;

                StreamGraphs::Edge main_e(H10->node2main[e.u], H10->node2main[e.v]);

                // check that edge is in graph
                REQUIRE(H10->main_graph.hasEdge(main_e.u, main_e.v));
                REQUIRE(H10->main_graph.weight(main_e.u, main_e.v) == weight);

            }
            ++line_idx;

        }       
    }

}


