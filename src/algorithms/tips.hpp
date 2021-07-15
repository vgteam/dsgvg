#pragma once

#include "utils.hpp"
#include "algorithms/stepindex.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"
#include "odgi.hpp"
#include <omp.h>

/**
 * \file tips.hpp
 *
 * Defines algorithms for identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s).
 */

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

		/// Iterate over the given paths. We walk from a chosen end (path_start or path_back) until we hit a node of the given query (reference) path(s).
		/// We record the hit as a BED output.
		/// #chrom #start #end #median_range #path_name #path_pos
		/// #chrom: The query path name.
		/// #start: The 0-based start position of the query we hit in the node.
		/// #end: The 1-based end position of the query we hit in the node.
		/// #median_range: The 0-based median of the whole query path range of the node we hit. It is possible that a node contains several steps, so we want to mirror that here.
		/// #path: The name of the path we walked.
		/// #path_pos: The 0-based position of the path we walked when we hit the node of the query path.
		void walk_tips(const graph_t& graph,
				 const std::vector<path_handle_t>& paths,
				 const path_handle_t& target_path_t,
				 const std::vector<bool>& target_handles,
				 algorithms::step_index_t& step_index,
				 const uint64_t& num_threads,
				 const std::function<step_handle_t(const path_handle_t&)>& get_path_end,
				 const std::function<step_handle_t(const step_handle_t&)>& get_step,
				 const std::function<bool(const step_handle_t&)>& has_step,
				 algorithms::tips_bed_writer& bed_writer_thread,
				 const bool& progress,
				 const bool& walk_from_front,
				 ska::flat_hash_set<std::string>& not_visited_set);
	}
}
