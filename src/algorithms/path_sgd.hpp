#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <random>
#include <set>
#include <thread>
#include <atomic>
#include <handlegraph/path_handle_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include "xp.hpp"
#include "sgd_term.hpp"
#include "IITree.h"
#include <zipfian_int_distribution.h>
#include <iomanip>
#include <string>
#include "weakly_connected_components.hpp"
#include <sdsl/bit_vectors.hpp>
#include "dirty_zipfian_int_distribution.h"
#include "XoshiroCpp.hpp"
#include "progress.hpp"
#include "algorithms/layout.hpp"

#include <fstream>

namespace odgi {
namespace algorithms {

using namespace handlegraph;

struct handle_layout_t {
    uint64_t weak_component = 0;
    double pos = 0;
    handle_t handle = as_handle(0);
};

/// use SGD driven, by path guided, and partly zipfian distribution sampled pairwise distances to obtain a 1D linear layout of the graph that respects its topology
std::vector<double> path_linear_sgd(const graph_t &graph,
                                    const xp::XP &path_index,
                                    const std::vector<path_handle_t>& path_sgd_use_paths,
                                    const uint64_t &iter_max,
                                    const uint64_t &iter_with_max_learning_rate,
                                    const uint64_t &min_term_updates,
                                    const double &delta,
                                    const double &eps,
                                    const double &eta_max,
                                    const double &theta,
                                    const uint64_t &space,
                                    const uint64_t &space_max,
                                    const uint64_t &space_quantization_step,
                                    const uint64_t &nthreads,
                                    const bool &progress,
                                    const bool &snapshot,
                                    std::vector<std::string> &snapshots);

/// our learning schedule
std::vector<double> path_linear_sgd_schedule(const double &w_min,
                                             const double &w_max,
                                             const uint64_t &iter_max,
                                             const uint64_t &iter_with_max_learning_rate,
                                             const double &eps);

std::vector<handle_t> path_linear_sgd_order(const graph_t &graph,
                                            const xp::XP &path_index,
                                            const std::vector<path_handle_t>& path_sgd_use_paths,
                                            const uint64_t &iter_max,
                                            const uint64_t &iter_with_max_learning_rate,
                                            const uint64_t &min_term_updates,
                                            const double &delta,
                                            const double &eps,
                                            const double &eta_max,
                                            const double &theta,
                                            const uint64_t &space,
                                            const uint64_t &space_max,
                                            const uint64_t &space_quantization_step,
                                            const uint64_t &nthreads,
                                            const bool &progress,
                                            const std::string &seed,
                                            const bool &snapshot,
                                            const std::string &snapshot_prefix,
                                            const bool &write_layout,
                                            const std::string &layout_out);

}

}
