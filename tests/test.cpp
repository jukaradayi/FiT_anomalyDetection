#define CONFIG_CATCH_MAIN
//#include "catch.hpp"
#include <catch2/catch_test_macros.hpp>
#include "history_graph.hpp"
#include "H_graph.hpp"
#include "G_graph.hpp"
#include "globals.hpp"
#include "csv.hpp"
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

    SECTION("Check G queue size") {
    //StreamGraphs::Interaction i(std::stoi((node_loop)[0]), std::stoi((node_loop)[1]), std::stoi((node_loop)[2]));
    NetworKit::Graph main_graph(0, true, false);
    NetworKit::Graph top_graph(0, true, false);
    NetworKit::Graph bot_graph(0, true, false);

    GGraph* G20 = new GGraph(main_graph, top_graph, bot_graph, false, false, true, 1000, 1000, 100, 20); //TODO correct node set size

    int line_idx = 0;
    for(auto& main_loop: CSVRange(file))
    {
        StreamGraphs::Interaction i(std::stoi((main_loop)[0]), std::stoi((main_loop)[1]), std::stoi((main_loop)[2]));

        line_idx += 1;
        G20->updateGraph(i); // update graph
        if (line_idx == 24) {
//REQUIRE( hist_graph->queue.size() == 10);
            node u_7 = G20->node2main[7];
            node u_14 = G20->node2main[14];

            REQUIRE(G20->main_graph.degree(u_7) == 7);
            REQUIRE(G20->main_graph.degree(u_14) == 1);
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
//REQUIRE( hist_graph->queue.size() == 10);
            node u_7 = H10->node2main[7];
            node u_14 = H10->node2main[14];

            REQUIRE(H10->main_graph.degree(u_7) == 7);
            REQUIRE(H10->main_graph.degree(u_14) == 1);
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

}


