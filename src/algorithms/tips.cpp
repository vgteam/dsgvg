#include "tips.hpp"

namespace odgi {
	namespace algorithms {

		using namespace handlegraph;

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
					   ska::flat_hash_set<std::string>& not_visited_set) {

			const std::string target_path = graph.get_path_name(target_path_t);

			std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
			if (progress) {
				std::string end;
				if (walk_from_front) {
					end = "front";
				} else {
					end = "back";
				}
				std::string progress_message = "[odgi::tips::walk_tips] BED Progress for query path '"
											   + graph.get_path_name(target_path_t) + "' walking from the " + end + " of all paths:";
				progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
						paths.size(), progress_message);
			}

#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
			for (auto path : paths) {
				// prevent self tips
				if (path == target_path_t) {
					continue;
				}
				std::string query_path_name = graph.get_path_name(path);
				bool tip_reached_query = false;
				// collecting tips from the front
				step_handle_t cur_step = get_path_end(path);
				handle_t cur_h = graph.get_handle_of_step(cur_step);
				while (!tip_reached_query) {
					// did we already hit the given reference path?
					if (target_handles[number_bool_packing::unpack_number(cur_h)]) {
						std::vector<uint64_t> target_path_start_positions;
						std::vector<uint64_t> target_path_end_positions;
						graph.for_each_step_on_handle(
								cur_h,
								[&](const step_handle_t& s) {
									if (graph.get_path_handle_of_step(s) == target_path_t) {
//#pragma omp critical (cout)
//									std::cerr << "[odgi::tips] error: We reached a tip from the front!" << std::endl;

										// we collect all positions
										// later we will get the smallest, the highest and the median from these
										uint64_t pos = step_index.get_position(s);
										target_path_start_positions.push_back(pos);
										target_path_end_positions.push_back(pos + graph.get_length(cur_h) - 1); // 0-based
									}
								});
//#pragma omp critical (cout)
//					std::cerr << "target_path_start_positions.size(): " << target_path_start_positions.size() << std::endl;
						std::sort(target_path_start_positions.begin(),
								  target_path_start_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						std::sort(target_path_end_positions.begin(),
								  target_path_end_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						std::vector<uint64_t> target_path_positions;
						target_path_positions.insert(target_path_positions.end(), target_path_start_positions.begin(), target_path_start_positions.end());
						target_path_positions.insert(target_path_positions.end(), target_path_end_positions.begin(), target_path_end_positions.end());
						std::sort(target_path_positions.begin(),
								  target_path_positions.end(),
								  [&](const uint64_t & pos_a,
									  const uint64_t & pos_b) {
									  return pos_a < pos_b;
								  });
						double target_pos_median = utils::median_of_sorted_vec(target_path_positions);
						uint64_t target_min_pos = target_path_start_positions[0]; // 0-based starting position in BED
						uint64_t target_max_pos = target_path_end_positions[target_path_end_positions.size() - 1] + 1; // 1-based ending position in BED
//#pragma omp critical (cout)
						//std::cout << target_path << "\t" << target_min_pos << "\t" << target_max_pos << "\t"
						//		  << target_pos_median << "\t" << graph.get_path_name(path) << "\t" << step_index.get_position(cur_step) << std::endl;
						/// add BED record to queue of BED_writer_thread
						bed_writer_thread.append(target_path, target_min_pos, target_max_pos,
												 target_pos_median, query_path_name, step_index.get_position(cur_step), walk_from_front);
						tip_reached_query = true;
					}
					if (has_step(cur_step)) {
						cur_step = get_step(cur_step);
						cur_h = graph.get_handle_of_step(cur_step);
					} else {
						// did we iterate over all steps and we did not hit the query path?
						tip_reached_query = true;
#pragma omp critical (not_visited_set)
						not_visited_set.insert(query_path_name);
						// do we still want this?
//						if (progress) {
//#pragma omp critical (cout)
//							std::cerr << "[odgi::tips::walk_tips] warning: For target path '" << query_path_name << "' there was no hit!" << std::endl;
//						}
					}
				}
				if (progress) {
					progress_meter->increment(1);
				}
			}
			if (progress) {
				progress_meter->finish();
			}
		}

	}
}
