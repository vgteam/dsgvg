#include "untangle.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;

std::vector<step_handle_t> untangle_cuts(
    const PathHandleGraph& graph,
    const step_handle_t& _start,
    const step_handle_t& _end,
    const std::function<uint64_t(const step_handle_t&)>& get_step_pos,
    const std::function<bool(const handle_t&)>& is_cut) {
    auto path = graph.get_path_handle_of_step(_start);
    auto path_name = graph.get_path_name(path);
    // this assumes that the end is not inclusive
    /*
    std::cerr << "untangle_cuts(" << path_name << ", "
              << start_pos << ", "
              << end_pos << ")" << std::endl;
    */
    std::vector<step_handle_t> cut_points;
    std::deque<std::pair<step_handle_t, step_handle_t>> todo;
    todo.push_back(std::make_pair(_start, _end));
    while (!todo.empty()) {
        auto start = todo.front().first;
        auto end = todo.front().second;
        uint64_t start_pos = get_step_pos(start);
        uint64_t end_pos = get_step_pos(end);
        cut_points.push_back(start);
        todo.pop_front();
        // we go forward until we see a loop, where the other step has position < end_pos and > start_pos
        for (step_handle_t step = start; step != end; step = graph.get_next_step(step)) {
            //  we take the first and shortest loop we find
            // TODO change this, it can be computed based on the node length
            auto curr_pos = get_step_pos(step);
            handle_t handle = graph.get_handle_of_step(step);
            if (is_cut(handle)) {
                cut_points.push_back(step);
            }
            bool found_loop = false;
            step_handle_t other;
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& s) {
                    if (step != s // not the step we're on
                        && graph.get_path_handle_of_step(s) == path) {
                        auto other_pos = get_step_pos(s);
                        if (other_pos > start_pos
                            && other_pos < end_pos
                            && (!found_loop && other_pos > curr_pos
                                || (found_loop && (get_step_pos(other) > other_pos)))) {
                            found_loop = true;
                            other = s;
                        }
                    }
                });
            if (found_loop) {
                //  recurse this function into it, taking start as our current handle other side of the loop as our end
                //  to cut_points we add the start position, the result from recursion, and our end position
                todo.push_back(std::make_pair(step, other));
                //  we then step forward to the loop end and continue iterating
                step = other;
            }
        }
        // TODO this block is the same as the previous one, but in reverse
        // the differences in how positions are managed are subtle, making it hard to fold the
        // forward and reverse version together
        // now we reverse it
        step_handle_t path_begin = graph.path_begin(path);
        if (end == path_begin || !graph.has_previous_step(end)) {
            return cut_points;
        }
        //std::cerr << "reversing" << std::endl;
        for (step_handle_t step = end;
             get_step_pos(step) > start_pos;
             step = graph.get_previous_step(step)) {
            //  we take the first and shortest loop we find
            // TODO change this, it can be computed based on the node length
            auto curr_pos = get_step_pos(step);
            handle_t handle = graph.get_handle_of_step(step);
            if (is_cut(handle)) {
                cut_points.push_back(step);
            }
            //std::cerr << "rev on step " << graph.get_id(handle) << " " << curr_pos << std::endl;
            bool found_loop = false;
            step_handle_t other;
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& s) {
                    if (step != s // not the step we're on
                        && graph.get_path_handle_of_step(s) == path) {
                        auto other_pos = get_step_pos(s);
                        if (other_pos > start_pos
                            && other_pos < end_pos
                            && other_pos < curr_pos
                            && (!found_loop
                                || (found_loop && (get_step_pos(other) < other_pos)))) {
                            found_loop = true;
                            other = s;
                        }
                    }
                });
            if (found_loop) {
                //  recurse this function into it, taking start as our current handle other side of the loop as our end
                //  to cut_points we add the start position, the result from recursion, and our end position
                todo.push_back(std::make_pair(other, step));
                //  we then step forward to the loop end and continue iterating
                step = other;
            }
        }
        cut_points.push_back(end);
    }
    // and sort
    std::sort(cut_points.begin(),
              cut_points.end(),
              [&](const step_handle_t& a,
                  const step_handle_t& b) {
                  return get_step_pos(a) < get_step_pos(b);
              });
    // then take unique positions
    cut_points.erase(std::unique(cut_points.begin(),
                                 cut_points.end()),
                     cut_points.end());
    return cut_points;
}

void write_cuts(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const std::function<uint64_t(const step_handle_t&)>& get_step_pos) {
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tcut" << std::endl;
    for (auto& step : cuts) {
        std::cout << path_name << "\t" << get_step_pos(step) << std::endl;
    }
}

std::vector<step_handle_t> merge_cuts(
    const std::vector<step_handle_t>& cuts,
    const uint64_t& dist,
    const std::function<uint64_t(const step_handle_t&)>& get_step_pos) {
    std::vector<step_handle_t> merged;
    uint64_t last = 0;
    //std::cerr << "dist is " << dist << std::endl;
    for (auto& step : cuts) {
        auto pos = get_step_pos(step);
        if (pos == 0 || pos > (last + dist)) {
            merged.push_back(step);
            last = pos;
        }
    }
    return merged;
}

void self_dotplot(
    const PathHandleGraph& graph,
    const path_handle_t& path) {
    auto step_pos = make_step_index(graph, { path }, 1);
    auto path_name = graph.get_path_name(path);
    std::cout << "name\tfrom\tto" << std::endl;
    graph.for_each_step_in_path(
        path,
        [&](const step_handle_t& step) {
            // TODO this is given by the walk, no need for a hash table lookup
            const auto& curr_pos = step_pos.find(step)->second;
            handle_t handle = graph.get_handle_of_step(step);
            graph.for_each_step_on_handle(
                handle,
                [&](const step_handle_t& s) {
                    if (graph.get_path_handle_of_step(s) == path) {
                        const auto& other_pos = step_pos.find(s)->second;
                        std::cout << path_name << "\t"
                                  << curr_pos << "\t"
                                  << other_pos << std::endl;
                        /*
                        if (other_pos > curr_pos
                            && other_pos > start_pos
                            && other_pos < end_pos
                            && other_pos - start_pos > min_length) {
                        }
                        */
                    }
                });
        });
}

ska::flat_hash_map<step_handle_t, uint64_t> make_step_index(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths,
    const size_t& num_threads) {
    ska::flat_hash_map<step_handle_t, uint64_t> step_pos;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (auto& path : paths) {
        uint64_t pos = 0;
        graph.for_each_step_in_path(
            path,
            [&](const step_handle_t& step) {
#pragma omp critical (step_pos)
                step_pos[step] = pos;
                handle_t handle = graph.get_handle_of_step(step);
                pos += graph.get_length(handle);
            });
#pragma omp critical (step_pos)
        step_pos[graph.path_end(path)] = pos; // record the end position
    }
    return step_pos;
}

void show_steps(
    const PathHandleGraph& graph,
    const ska::flat_hash_map<step_handle_t, uint64_t>& steps) {
    for (auto& pos : steps) {
        auto h = graph.get_handle_of_step(pos.first);
        auto p = graph.get_path_handle_of_step(pos.first);
        auto name = graph.get_path_name(p);
        auto id = graph.get_id(h);
        auto rev = graph.get_is_reverse(h);
        std::cerr << name << " " << id << (rev?"-":"+") << " " << pos.second << std::endl;
    }
}

// compute the reference segmentations
// and map them onto the graph using a static multiset index structure based on two arrays
// we'll build up a big vector of node -> path segment pairings
// then we'll sort them and build an index
segment_map_t::segment_map_t(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& paths,
    const std::function<uint64_t(const step_handle_t&)>& get_step_pos,
    const std::function<bool(const handle_t&)>& is_cut,
    const uint64_t& merge_dist,
    const size_t& num_threads) {
    std::vector<std::pair<uint64_t, int64_t>> node_to_segment;
    std::vector<std::vector<step_handle_t>> all_cuts(paths.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (uint64_t i = 0; i < paths.size(); ++i) { //auto& path : paths) {
        auto& path = paths[i];
        all_cuts[i] =
            merge_cuts(
                untangle_cuts(graph,
                              graph.path_begin(path),
                              graph.path_back(path),
                              get_step_pos,
                              is_cut
                    ),
                merge_dist,
                get_step_pos);
    }
    // the index construction must be serial
    for (uint64_t i = 0; i < paths.size(); ++i) {
        auto& path = paths[i];
        auto& cuts = all_cuts[i];
        //std::cerr << "reference segmentation" << std::endl;
        //write_cuts(graph, path, cuts, step_pos);
        // walk the path to get the segmentation
        uint64_t curr_segment_idx = 0;
        uint64_t segment_idx = segment_cut.size();
        uint64_t* curr_length = nullptr;
        for (step_handle_t step = graph.path_begin(path);
             step != graph.path_end(path);
             step = graph.get_next_step(step)) {
            // if we are at a segment cut
            if (step == cuts[curr_segment_idx]) {
                segment_idx = segment_cut.size();
                segment_cut.push_back(step);
                segment_length.push_back(0);
                curr_length = &segment_length.back();
                ++curr_segment_idx;
            }
            handle_t h = graph.get_handle_of_step(step);
            bool is_rev = graph.get_is_reverse(h);
            node_to_segment.push_back(
                std::make_pair(graph.get_id(h),
                               (is_rev ? -segment_idx : segment_idx)));
            uint64_t node_length = graph.get_length(h);
            *curr_length += node_length;
        }
    }
    //std::cerr << "segment_cut.size() " << segment_cut.size() << std::endl;
    //std::cerr << "segment_length.size() " << segment_length.size() << std::endl;
    ips4o::parallel::sort(node_to_segment.begin(),
                          node_to_segment.end(),
                          std::less<>(),
                          num_threads);
    // make the mapping
    uint64_t prev_node = 0;
    for (auto& node_segment : node_to_segment) {
        auto& node_id = node_segment.first;
        if (node_id > prev_node) {
            while (prev_node < node_id) {
                node_idx.push_back(segments.size());
                ++prev_node;
            }
        }
        segments.push_back(node_segment.second);
        prev_node = node_id;
    }
    auto max_id = graph.get_node_count();
    while (prev_node < max_id) {
        node_idx.push_back(segments.size());
        ++prev_node;
    }
    node_idx.push_back(segments.size()); // to avoid special casing the last node
}

void segment_map_t::for_segment_on_node(
    uint64_t node_id,
    const std::function<void(const uint64_t& segment_id, const bool& is_rev)>& func) const {
    uint64_t from = node_idx[node_id-1];
    uint64_t to = node_idx[node_id];
    for (uint64_t i = from; i < to; ++i) {
        auto& j = segments[i];
        func(std::abs(j), j < 0);
    }
}

uint64_t segment_map_t::get_segment_length(const uint64_t& segment_id) const {
    return segment_length.at(segment_id);
}

struct isec_t {
    uint64_t len = 0;
    uint64_t inv = 0;
    void incr(const uint64_t& l, const bool& is_inv) {
        len += l;
        inv += (is_inv ? l : 0);
    }
};

std::vector<segment_mapping_t>
segment_map_t::get_matches(
        const PathHandleGraph& graph,
        const step_handle_t& start,
        const step_handle_t& end,
        const uint64_t& query_length) const {
    // collect the target segments that overlap our segment
    // computing the intersection size (in bp) as we go
    // our final metric is jaccard of intersection over total length for each overlapped target
    path_handle_t query_path = graph.get_path_handle_of_step(start);
    ska::flat_hash_map<uint64_t, isec_t> target_isec;
    for (step_handle_t step = start;
         step != end;
         step = graph.get_next_step(step)) {
        handle_t h = graph.get_handle_of_step(step);
        uint64_t node_id = graph.get_id(h);
        uint64_t node_length = graph.get_length(h);
        bool is_rev = graph.get_is_reverse(h);
        for_segment_on_node(
            node_id,
            [&](const uint64_t& segment_id, const bool& segment_rev) {
                // n.b. we do not skip self matches
                target_isec[segment_id].incr(node_length, is_rev != segment_rev);
            });
    }
    // compute the jaccards
    std::vector<segment_mapping_t> jaccards;
    for (auto& p : target_isec) {
        auto& segment_id = p.first;
        auto& isec = p.second.len;
        //auto& inv = p.second.inv;
        bool is_inv = (double)p.second.inv/(double)isec > 0.5;
        // intersection / union
        jaccards.push_back(
            {
                segment_id,
                is_inv,
                (double)isec
                / (double)(get_segment_length(segment_id)
                           + query_length - isec)
            });
    }
    // sort the target segments by their jaccard with the query
    std::sort(jaccards.begin(), jaccards.end(),
              [](const segment_mapping_t& a,
                 const segment_mapping_t& b) {
                  return std::tie(a.jaccard, a.segment_id, a.is_inv) >
                      std::tie(b.jaccard, b.segment_id, b.is_inv);
              });
    return jaccards;
}

const step_handle_t& segment_map_t::get_segment_cut(
    const uint64_t& idx) const {
    return segment_cut[idx];
}

double self_mean_coverage(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const step_handle_t& begin,
    const step_handle_t& end) {
    uint64_t sum = 0;
    uint64_t bp = 0;
    for (step_handle_t step = begin;
         step != end;
         step = graph.get_next_step(step)) {
        handle_t handle = graph.get_handle_of_step(step);
        uint64_t len = graph.get_length(handle);
        bp += len;
        graph.for_each_step_on_handle(
            handle,
            [&](const step_handle_t& s) {
                if (graph.get_path_handle_of_step(s) == path) {
                    sum += len;
                }
            });
    }
    return (double)sum / (double)bp;
}

void map_segments(
    const PathHandleGraph& graph,
    const path_handle_t& path,
    const std::vector<step_handle_t>& cuts,
    const segment_map_t& target_segments,
    const std::function<uint64_t(const step_handle_t&)>& get_step_pos) {
    std::string query_name = graph.get_path_name(path);
    for (uint64_t i = 0; i < cuts.size()-1; ++i) {
        auto& begin = cuts[i];
        auto& end = cuts[i+1];
        auto begin_pos = get_step_pos(begin);
        auto end_pos = get_step_pos(end);
        uint64_t length = end_pos - begin_pos;
        std::vector<segment_mapping_t> target_mapping =
            target_segments.get_matches(graph, cuts[i], cuts[i+1], length);
        if (!target_mapping.empty()) {
            // get the self coverage TODO
            double self_coverage = self_mean_coverage(graph, path, begin, end);
            auto& best = target_mapping.front();
            auto& score = best.jaccard;
            auto& idx = best.segment_id; // segment index
            auto& target_begin = target_segments.get_segment_cut(idx);
            auto target_begin_pos = get_step_pos(target_begin);
            auto target_end_pos = target_begin_pos + target_segments.get_segment_length(idx);
            path_handle_t target_path = graph.get_path_handle_of_step(target_begin);
            std::string target_name = graph.get_path_name(target_path);
#pragma omp critical (cout)
            std::cout << query_name << "\t"
                      << begin_pos << "\t"
                      << end_pos << "\t"
                      << target_name << "\t"
                      << target_begin_pos << "\t"
                      << target_end_pos << "\t"
                      << score << "\t"
                //<< "+" << "\t" // the query is always in the positive frame
                      << (best.is_inv ? "-" : "+") << "\t"
                      << self_coverage
                      << std::endl;
                //"\t" << target_mapping.size() << std::endl;
            // todo: orientation
        }
    }
}

// BEDPE (pair-BED) projection of the graph
// that describe nonlinear query : target relationships
void untangle(
    const PathHandleGraph& graph,
    const std::vector<path_handle_t>& queries,
    const std::vector<path_handle_t>& targets,
    const uint64_t& merge_dist,
    const size_t& num_threads,
    const bool progress) {

    std::cerr << "[odgi::algorithms::untangle] untangling " << queries.size() << " queries with " << targets.size() << " targets" << std::endl;
    std::vector<path_handle_t> paths;
    paths.insert(paths.end(), queries.begin(), queries.end());
    paths.insert(paths.end(), targets.begin(), targets.end());
    std::sort(paths.begin(), paths.end());
    paths.erase(std::unique(paths.begin(), paths.end()),
                paths.end());
    //std::cerr << "[odgi::algorithms::untangle] building step index" << std::endl;
    //auto step_pos = make_step_index(graph, paths, num_threads);
    step_index_t step_index(graph, paths, num_threads, progress);
    auto get_step_pos = [&](const step_handle_t& step) {
        return step_index.get_position(step);
    };
    //std::cerr << "[odgi::algorithms::untangle] step index contains " << step_pos.size() << " steps" << std::endl;
    // collect all possible cuts
    // we'll use this to drive the subsequent segmentation
    std::cerr << "[odgi::algorithms::untangle] establishing initial cuts" << std::endl;
    atomicbitvector::atomic_bv_t cut_nodes(graph.get_node_count()+1);
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (auto& path : paths) {
        std::vector<step_handle_t> cuts
            = merge_cuts(
                untangle_cuts(graph,
                              graph.path_begin(path),
                              graph.path_back(path),
                              get_step_pos,
                              [](const handle_t& h) { return false; }),
                merge_dist,
                get_step_pos);
        for (auto& step : cuts) {
            cut_nodes.set(graph.get_id(graph.get_handle_of_step(step)));
        }
        //std::cerr << "setup" << std::endl;
        //write_cuts(graph, path, cuts, step_pos);
    }

    //auto step_pos = make_step_index(graph, queries);
    // node to reference segmentation mapping
    std::cerr << "[odgi::algorithms::untangle] building target segment index" << std::endl;
    segment_map_t target_segments(graph,
                                  targets,
                                  get_step_pos,
                                  [&cut_nodes,&graph](const handle_t& h) {
                                      return cut_nodes[graph.get_id(h)] != 0;
                                  },
                                  merge_dist,
                                  num_threads);

    //show_steps(graph, step_pos);
    //std::cout << "path\tfrom\tto" << std::endl;
    //auto step_pos = make_step_index(graph, queries);
    std::cerr << "[odgi::algorithms::untangle] writing pair BED for " << queries.size() << " queries" << std::endl;
    std::cout << "#query.name\tquery.start\tquery.end\tref.name\tref.start\tref.end\tscore\tinv\tself.cov" << std::endl;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
    for (auto& query : queries) {
        std::vector<step_handle_t> cuts
            = merge_cuts(
                untangle_cuts(graph,
                              graph.path_begin(query),
                              graph.path_back(query),
                              get_step_pos,
                              [&cut_nodes,&graph](const handle_t& h) {
                                  return cut_nodes[graph.get_id(h)] != 0;
                              }),
                merge_dist,
                get_step_pos);
        map_segments(graph, query, cuts, target_segments, get_step_pos);

        //write_cuts(graph, query, cuts, step_pos);
    }
    //self_dotplot(graph, query, step_pos);
}

}
}
