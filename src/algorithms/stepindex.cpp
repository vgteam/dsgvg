#include "stepindex.hpp"
#include "progress.hpp"

namespace odgi {
namespace algorithms {

step_index_t::step_index_t(const PathHandleGraph& graph,
                           const std::vector<path_handle_t>& paths,
                           const uint64_t& nthreads,
                           const bool progress) {
    // iterate through the paths, recording steps in the structure we'll use to build the mphf
    std::vector<step_handle_t> steps;
    for (auto& path : paths) {
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) { steps.push_back(step); });
        steps.push_back(graph.path_end(path));
    }
    // sort the steps
    ips4o::parallel::sort(steps.begin(), steps.end(), std::less<>(), nthreads);
    // build the hash function (quietly)
    step_mphf = new boophf_step_t(steps.size(), steps, nthreads, 2.0, false, false);
    // use the hash function to record the step positions
    pos.resize(steps.size());
	std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
	if (progress) {
		progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
				paths.size(), "[odgi::algorithms::stepindex] Building Progress:");
	}
#pragma omp parallel for schedule(dynamic,1)
    for (auto& path : paths) {
        uint64_t offset = 0;
        graph.for_each_step_in_path(
            path, [&](const step_handle_t& step) {
                pos[step_mphf->lookup(step)] = offset;
                offset += graph.get_length(graph.get_handle_of_step(step));
            });
        pos[step_mphf->lookup(graph.path_end(path))] = offset;
        if (progress) {
        	progress_meter->increment(1);
        }
    }
	if (progress) {
		progress_meter->finish();
	}
}

const uint64_t& step_index_t::get_position(const step_handle_t& step) {
    return pos[step_mphf->lookup(step)];
}

step_index_t::~step_index_t(void) {
    delete step_mphf;
}

}
}
