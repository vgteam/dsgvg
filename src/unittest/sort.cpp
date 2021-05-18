#include "catch.hpp"

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "odgi.hpp"
#include "algorithms/topological_sort.hpp"

#include <iostream>
#include <limits>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <xp.hpp>
#include <path_sgd.hpp>

namespace odgi {
namespace unittest {

using namespace std;
using namespace handlegraph;

TEST_CASE("Sorting a graph without paths", "[sort]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    // sort the graph
    auto order = algorithms::topological_order(&graph);
    graph.apply_ordering(order, true);
    SECTION("The graph is as expected when sorted") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "TTG");
        REQUIRE(graph.get_node_count() == 6);
    }
    std::reverse(order.begin(), order.end());
    graph.apply_ordering(order, true);
    SECTION("The graph is as expected when reversed") {
        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "TTG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "CAAATAAG");
        REQUIRE(graph.get_node_count() == 6);
    }
}

TEST_CASE("Sorting a graph with paths", "[sort]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    std::vector<path_handle_t> paths;
    for (uint64_t i = 0; i < 10; ++i) {
        auto path = graph.create_path_handle("x" + std::to_string(i));
        graph.append_step(path, n1);
        graph.append_step(path, n2);
        graph.append_step(path, n3);
        graph.append_step(path, n4);
        graph.append_step(path, n5);
        graph.append_step(path, n6);
        paths.push_back(path);
    }

    auto test_path =
        [&](const path_handle_t& p) {
            auto& path_meta = graph.path_metadata(p);
            uint64_t i = 0;
            graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                handle_t h = graph.get_handle_of_step(step);
                ++i;
            });
            REQUIRE(i == path_meta.length);
        };

    graph.for_each_path_handle(test_path);

    // sort the graph
    auto order = algorithms::topological_order(&graph);
    graph.apply_ordering(order, true);

    graph.for_each_path_handle(test_path);

    SECTION("The graph is as expected when sorted") {
        REQUIRE(graph.get_path_count() == 10);

        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "TTG");
        REQUIRE(graph.get_node_count() == 6);
        for (auto& path : paths) {
            REQUIRE(graph.get_handle_of_step(graph.path_begin(path)) == n1);
            std::string x;
            graph.for_each_step_in_path(
                path,
                [&](const step_handle_t& step) {
                    x.append(graph.get_sequence(graph.get_handle_of_step(step)));
                });
            REQUIRE(x == "CAAATAAGAGTCTTG");
        }
    }
    std::reverse(order.begin(), order.end());
    graph.apply_ordering(order, true);

    graph.for_each_path_handle(test_path);

    SECTION("The graph is as expected when reversed") {
        REQUIRE(graph.get_path_count() == 10);

        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "TTG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "CAAATAAG");
        REQUIRE(graph.get_node_count() == 6);
    }
}

TEST_CASE("Sorting a graph with paths 1 node long", "[sort]") {
    graph_t graph;

    std::vector<handle_t> handles;
    std::vector<path_handle_t> paths;
    for (uint64_t i = 0; i < 10; ++i) {
        auto node = graph.create_handle("C");
        auto path = graph.create_path_handle("x" + std::to_string(i));

        graph.append_step(path, node);

        handles.push_back(node);
        paths.push_back(path);
    }

    auto test_path =
            [&](const path_handle_t& p) {
                auto& path_meta = graph.path_metadata(p);
                uint64_t i = 0;
                graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                    handle_t h = graph.get_handle_of_step(step);
                    ++i;
                });
                REQUIRE(i == path_meta.length);
            };

    graph.for_each_path_handle(test_path);

    // sort the graph
    uint64_t path_sgd_iter_max = 30;
    double path_sgd_zipf_theta = 0.99;
    double path_sgd_eps = 0.01;
    double path_sgd_delta = 0;
    std::vector<handlegraph::path_handle_t> path_sgd_use_paths;
    graph.for_each_path_handle(
            [&](const handlegraph::path_handle_t &path) {
                path_sgd_use_paths.push_back(path);
            });

    // path length interrogation
    std::function<uint64_t(const std::vector<handlegraph::path_handle_t> &,
                           const xp::XP &)> get_sum_path_step_count
            = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t sum_path_step_count = 0;
                for (auto& path : path_sgd_use_paths) {
                    sum_path_step_count += path_index.get_path_step_count(path);
                }
                return sum_path_step_count;
            };
    std::function<uint64_t(const std::vector<handlegraph::path_handle_t> &,
                           const xp::XP &)> get_max_path_step_count
            = [&](const std::vector<handlegraph::path_handle_t> &path_sgd_use_paths, const xp::XP &path_index) {
                uint64_t max_path_step_count = 0;
                for (auto& path : path_sgd_use_paths) {
                    max_path_step_count = std::max(max_path_step_count, path_index.get_path_step_count(path));
                }
                return max_path_step_count;
            };

    xp::XP path_index;
    path_index.from_handle_graph(graph, 1);

    uint64_t sum_path_step_count = get_sum_path_step_count(path_sgd_use_paths, path_index);
    uint64_t path_sgd_min_term_updates = sum_path_step_count;//p_sgd_min_term_updates * sum_path_step_count;
    uint64_t max_path_step_count = get_max_path_step_count(path_sgd_use_paths, path_index);
    uint64_t path_sgd_zipf_space = std::min((uint64_t)10000, max_path_step_count);
    double path_sgd_max_eta = max_path_step_count * max_path_step_count;
    uint64_t path_sgd_zipf_space_max = 1000;
    uint64_t path_sgd_zipf_space_quantization_step = 100;
    std::string path_sgd_seed = "pangenomic!";

    uint64_t path_sgd_iter_max_learning_rate = 0; // don't use this max iter stuff

    auto order = odgi::algorithms::path_linear_sgd_order(
            graph,
            path_index,
            path_sgd_use_paths,
            path_sgd_iter_max,
            path_sgd_iter_max_learning_rate,
            path_sgd_min_term_updates,
            path_sgd_delta,
            path_sgd_eps,
            path_sgd_max_eta,
            path_sgd_zipf_theta,
            path_sgd_zipf_space,
            path_sgd_zipf_space_max,
            path_sgd_zipf_space_quantization_step,
            2,
            false,
            path_sgd_seed,
            false,
            "",
            false,
            ""
    );

    graph.apply_ordering(order, true);
    graph.for_each_path_handle(test_path);

    SECTION("The graph is as expected when sorted") {
        REQUIRE(graph.get_path_count() == 10);

        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "C");
        REQUIRE(graph.get_sequence(graph.get_handle(10)) == "C");
        REQUIRE(graph.get_node_count() == 10);
        for (uint64_t i = 0; i < paths.size(); ++i) {
            auto& node = handles[i];
            auto& path = paths[i];

            REQUIRE(graph.get_handle_of_step(graph.path_begin(path)) == node);

            std::string x;
            graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        x.append(graph.get_sequence(graph.get_handle_of_step(step)));
                    });
            REQUIRE(x == "C");
        }
    }
}

TEST_CASE("Sorting the paths in a graph", "[sort]") {
    graph_t graph;
    handle_t n1 = graph.create_handle("CAAATAAG");
    handle_t n2 = graph.create_handle("A");
    handle_t n3 = graph.create_handle("G");
    handle_t n4 = graph.create_handle("T");
    handle_t n5 = graph.create_handle("C");
    handle_t n6 = graph.create_handle("TTG");
    graph.create_edge(n1, n2);
    graph.create_edge(n2, n3);
    graph.create_edge(n3, n4);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    std::vector<path_handle_t> paths;
    for (uint64_t i = 0; i < 10; ++i) {
        auto path = graph.create_path_handle("x" + std::to_string(i));
        graph.append_step(path, n1);
        graph.append_step(path, n2);
        graph.append_step(path, n3);
        graph.append_step(path, n4);
        graph.append_step(path, n5);
        graph.append_step(path, n6);
        paths.push_back(path);
    }

    auto test_path =
            [&](const path_handle_t& p) {
                auto& path_meta = graph.path_metadata(p);
                uint64_t i = 0;
                graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
                    handle_t h = graph.get_handle_of_step(step);
                    ++i;
                });
                REQUIRE(i == path_meta.length);
            };

    graph.for_each_path_handle(test_path);

    std::reverse(paths.begin(), paths.end());
    graph.apply_path_ordering(paths);

    graph.for_each_path_handle(test_path);

    SECTION("The graph is as expected when the path order is reversed") {
        REQUIRE(graph.get_path_count() == 10);

        REQUIRE(graph.get_sequence(graph.get_handle(1)) == "CAAATAAG");
        REQUIRE(graph.get_sequence(graph.get_handle(6)) == "TTG");
        REQUIRE(graph.get_node_count() == 6);
        for (auto& path : paths) {
            REQUIRE(graph.get_handle_of_step(graph.path_begin(path)) == n1);
            std::string x;
            graph.for_each_step_in_path(
                    path,
                    [&](const step_handle_t& step) {
                        x.append(graph.get_sequence(graph.get_handle_of_step(step)));
                    });
            REQUIRE(x == "CAAATAAGAGTCTTG");
        }

        uint64_t i = 10;
        graph.for_each_path_handle([&](const path_handle_t& p) {
            i--;
            REQUIRE(graph.get_path_name(p) == "x" + std::to_string(i));
        });
        REQUIRE(i == 0);
    }
}
}
}
