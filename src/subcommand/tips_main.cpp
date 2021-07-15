#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include <omp.h>
#include "utils.hpp"
#include "algorithms/stepindex.hpp"
#include "algorithms/tips.hpp"
#include "algorithms/tips_bed_writer_thread.hpp"
#include "vector"
#include "hash_map.hpp"

namespace odgi {

	using namespace odgi::subcommand;

	int main_tips(int argc, char **argv) {

		// trick argumentparser to do the right thing with the subcommand
		for (uint64_t i = 1; i < argc - 1; ++i) {
			argv[i] = argv[i + 1];
		}
		std::string prog_name = "odgi tips";
		argv[0] = (char *) prog_name.c_str();
		--argc;

		args::ArgumentParser parser(
				"Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s). Prints BED records to stdout.");
		args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
		args::ValueFlag<std::string> og_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "input"});
		args::Group tips_opts(parser, "[ Tips Options ]");
		args::ValueFlag<std::string> _query_path(tips_opts, "NAME", "Use this query path.",
												 {'q', "query-path"});
		args::ValueFlag<std::string> _target_path(tips_opts, "NAME", "Use this target (reference) path.",
												  {'r', "target-path"});
		args::ValueFlag<std::string> _query_paths(tips_opts, "FILE", "Use query paths listed (one per line) in FILE.",
												  {'Q', "query-paths"});
		args::ValueFlag<std::string> _target_paths(tips_opts, "FILE", "Use target (reference) paths listed (one per line) in FILE.",
												   {'R', "target-paths"});
		args::ValueFlag<std::string> _not_visited_tsv(tips_opts, "FILE", "Write target path(s) that do not visit the query path(s) to this FILE.",
												   {'v', "not-visited-tsv"});
		args::Group threading(parser, "[ Threading ]");
		args::ValueFlag<uint64_t> nthreads(threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
		args::Group processing_info_opts(parser, "[ Processing Information ]");
		args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
		args::Group program_information(parser, "[ Program Information ]");
		args::HelpFlag help(program_information, "help", "Print a help message for odgi tips.", {'h', "help"});

		try {
			parser.ParseCLI(argc, argv);
		} catch (args::Help) {
			std::cout << parser;
			return 0;
		} catch (args::ParseError e) {
			std::cerr << e.what() << std::endl;
			std::cerr << parser;
			return 1;
		}
		if (argc == 1) {
			std::cout << parser;
			return 1;
		}

		if (!og_file) {
			std::cerr << "[odgi::tips] error: please specify a graph to walk the tips in via -i=[FILE], --idx=[FILE]."
					  << std::endl;
			return 1;
		}

		const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

		odgi::graph_t graph;
		assert(argc > 0);
		std::string infile = args::get(og_file);
		if (!infile.empty()) {
			if (infile == "-") {
				graph.deserialize(std::cin);
			} else {
				utils::handle_gfa_odgi_input(infile, "tips", args::get(progress), num_threads, graph);
			}
		}

		omp_set_num_threads(num_threads);

		// path loading
		auto load_paths = [&](const std::string& path_names_file) {
			std::ifstream path_names_in(path_names_file);
			uint64_t num_of_paths_in_file = 0;
			std::vector<bool> path_already_seen;
			path_already_seen.resize(graph.get_path_count(), false);
			std::string line;
			std::vector<path_handle_t> paths;
			while (std::getline(path_names_in, line)) {
				if (!line.empty()) {
					if (graph.has_path(line)) {
						const path_handle_t path = graph.get_path_handle(line);
						const uint64_t path_rank = as_integer(path) - 1;
						if (!path_already_seen[path_rank]) {
							path_already_seen[path_rank] = true;
							paths.push_back(path);
						} else {
							std::cerr << "[odgi::tips] error: in the path list there are duplicated path names."
									  << std::endl;
							exit(1);
						}
					}
					++num_of_paths_in_file;
				}
			}
			path_names_in.close();
			if (progress) {
				std::cerr << "[odgi::tips] found " << paths.size() << "/" << num_of_paths_in_file
						  << " paths to consider." << std::endl;
			}
			if (paths.empty()) {
				std::cerr << "[odgi::tips] error: no path to consider." << std::endl;
				exit(1);
			}
			return paths;
		};

		std::vector<path_handle_t> target_paths;
		std::vector<path_handle_t> query_paths;
		if (_target_path) {
			auto& path_name = args::get(_target_path);
			if (graph.has_path(path_name)) {
				target_paths.push_back(graph.get_path_handle(path_name));
			} else {
				std::cerr << "[odgi::tips] error: no target path '"
						  << path_name << "' found in graph." << std::endl;
				exit(1);
			}
		} else if (_target_paths) {
			target_paths = load_paths(args::get(_target_paths));
		} else {
			target_paths.reserve(graph.get_path_count());
			graph.for_each_path_handle([&](const path_handle_t path) {
				target_paths.push_back(path);
			});
		}
		if (_query_path) {
			auto& path_name = args::get(_query_path);
			if (graph.has_path(path_name)) {
				query_paths.push_back(graph.get_path_handle(path_name));
			} else {
				std::cerr << "[odgi::tips] error: no query path '"
						  << path_name << "' found in graph." << std::endl;
				exit(1);
			}
		} else if (_query_paths) {
			query_paths = load_paths(args::get(_query_paths));
		} else {
			query_paths.reserve(graph.get_path_count());
			graph.for_each_path_handle([&](const path_handle_t path) {
				query_paths.push_back(path);
			});
		}

		std::vector<path_handle_t> paths;
		paths.reserve(graph.get_path_count());
		graph.for_each_path_handle([&](const path_handle_t &path) {
			paths.push_back(path);
		});
		algorithms::step_index_t step_index(graph, paths, num_threads, progress);
		// open the bed writer thread
		algorithms::tips_bed_writer bed_writer_thread;
		bed_writer_thread.open_writer();
		ofstream not_visited_out;
		if (_not_visited_tsv) {
			not_visited_out = ofstream(args::get(_not_visited_tsv));
		}

		auto get_path_begin = [&](const path_handle_t& path) {
			return graph.path_begin(path);
		};
		auto get_next_step = [&](const step_handle_t& step) {
			return graph.get_next_step(step);
		};
		auto get_path_back = [&](const path_handle_t& path) {
			return graph.path_back(path);
		};
		auto get_prev_step = [&](const step_handle_t& step) {
			return graph.get_previous_step(step);
		};
		auto has_next_step = [&](const step_handle_t& step) {
			return graph.has_next_step(step);
		};
		auto has_previous_step = [&](const step_handle_t& step) {
			return graph.has_previous_step(step);
		};

		// TODO loop over query paths here
		for (auto target_path_t : target_paths) {
			// make bit vector across nodes to tell us if we have a hit
			// this is a speed up compared to iterating through all steps of a potential node for each walked step
			std::vector<bool> target_handles;
			target_handles.resize(graph.get_node_count(), false);
			graph.for_each_step_in_path(target_path_t, [&](const step_handle_t &step) {
				handle_t h = graph.get_handle_of_step(step);
				target_handles[number_bool_packing::unpack_number(h)] = true;
			});
			ska::flat_hash_set<std::string> not_visited_set;
			/// walk from the front
			algorithms::walk_tips(graph, query_paths, target_path_t, target_handles, step_index, num_threads, get_path_begin,
								  get_next_step, has_next_step, bed_writer_thread, progress, true, not_visited_set);
			std::vector<path_handle_t> visitable_query_paths;
			for (auto query_path : query_paths) {
				if (!not_visited_set.count(graph.get_path_name(query_path))) {
					visitable_query_paths.push_back(query_path);
				}
			}
			/// walk from the back
			algorithms::walk_tips(graph, visitable_query_paths, target_path_t, target_handles, step_index, num_threads, get_path_back,
								  get_prev_step, has_previous_step, bed_writer_thread, progress, false, not_visited_set);
			/// let's write our paths we did not visit
			std::string query_path = graph.get_path_name(target_path_t);
			for (auto not_visited_path : not_visited_set) {
				not_visited_out << query_path << "\t" << not_visited_path << std::endl;
			}
		}
		bed_writer_thread.close_writer();
		if (_not_visited_tsv) {
			not_visited_out.close();
		}

		exit(0);
	}

	static Subcommand odgi_tips("tips",
									"Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s). Prints BED records to stdout.",
									PIPELINE, 3, main_tips);

}
