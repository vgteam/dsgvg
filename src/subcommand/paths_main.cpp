#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "split.hpp"
#include "position.hpp"
#include <omp.h>
#include "utils.hpp"

namespace odgi {

using namespace odgi::subcommand;

int main_paths(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi paths";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("Interrogate the embedded paths of a graph. Does not print anything to stdout by default!");
    args::Group mandatory_opts(parser, "[ MANDATORY ARGUMENTS ]");
    args::ValueFlag<std::string> dg_in_file(mandatory_opts, "FILE", "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!", {'i', "idx"});
    args::Group path_investigation_opts(parser, "[ Path Investigation Options ]");
    args::ValueFlag<std::string> overlaps_file(path_investigation_opts, "FILE", "Read in the path grouping *FILE* to generate the overlap statistics"
                                                               " from. The file must be tab-delimited. The first column lists a"
                                                               " grouping and the second the path itself. Each line has one path entry."
                                                               " For each group the pairwise overlap statistics for each pairing will"
                                                               " be calculated and printed to stdout.", {'O', "overlaps"});
    args::Flag list_names(path_investigation_opts, "list-names", "Print the paths in the graph to stdout. Each path is printed in its"
                                                " own line.", {'L', "list-paths"});
	args::Flag list_path_start_end(path_investigation_opts, "list-path-start-end", "If -L,--list-paths was specified, this additionally prints the start and end positions of each path in additional, tab-delimited coloumns."
						   , {'l', "list-path-start-end"});
    args::ValueFlag<std::string> path_delim(path_investigation_opts, "CHAR", "The part of each path name before this delimiter CHAR is a group"
                                                                             " identifier. This parameter should only be set in combination with"
                                                                             " **-H, --haplotypes**. Prints an additional, first column"
                                                                             " **group.name** to stdout.", {'D', "delim"});
    args::Flag haplo_matrix(path_investigation_opts, "haplo", "Print to stdout the paths in an approximate binary haplotype matrix"
                                                              " based on the graph’s sort order. The output is tab-delimited:"
                                                              " *path.name*, *path.length*, *node.count*, *node.1*,"
                                                              " *node.2*, *node.n*. Each path entry is printed in its own line.", {'H', "haplotypes"});
    args::Flag distance_matrix(path_investigation_opts, "distance", "Provides a sparse distance matrix for paths. If **-D, --delim** is"
                                                                    " set, it will be path groups distances. Each line prints in a tab-delimited format to stdout:"
                                                                    " *path.a*, *path.b*, *path.a.length*, *path.b.length*, *intersection*, *jaccard*, *euclidean*." , {'d', "distance"});
    args::Flag write_fasta(path_investigation_opts, "fasta", "Print paths in FASTA format to stdout. One line for the FASTA header, another line for the whole sequence.", {'f', "fasta"});
    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> threads(threading_opts, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
	args::Group processing_info_opts(parser, "[ Processing Information ]");
	args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.", {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi paths.", {'h', "help"});

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
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!dg_in_file) {
        std::cerr << "[odgi::paths] error: please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

	if (list_path_start_end && !list_names) {
		std::cerr << "[odgi::paths] error: please specify also -L,--list-path with the -l,--list-path-start-end option!" << std::endl;
		return 1;
	}

	const uint64_t num_threads = args::get(threads) ? args::get(threads) : 1;
    omp_set_num_threads(num_threads);

	graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
			utils::handle_gfa_odgi_input(infile, "paths", args::get(progress), num_threads, graph);
        }
    }

    if (list_path_start_end && list_names) {
    	std::vector<path_handle_t> paths;
		graph.for_each_path_handle([&](const path_handle_t& p) {
			paths.push_back(p);
		});
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
		for (auto path : paths) {
			uint64_t path_len = 0;
			graph.for_each_step_in_path(path, [&](const step_handle_t& s) {
				handle_t h = graph.get_handle_of_step(s);
				path_len += graph.get_length(h);
			});
#pragma omp critical (cout)
			std::cout << graph.get_path_name(path) << "\t" << 1 << "\t" << path_len << std::endl;
		}
	} else if (args::get(list_names)) {
        graph.for_each_path_handle([&](const path_handle_t& p) {
                std::cout << graph.get_path_name(p) << std::endl;
            });
    }

    if (args::get(write_fasta)) {
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::cout << ">" << graph.get_path_name(p) << std::endl;
                graph.for_each_step_in_path(
                    p, [&](const step_handle_t& s) {
                           std::cout << graph.get_sequence(graph.get_handle_of_step(s));
                       });
                std::cout << std::endl;
            });
    }


    if (args::get(haplo_matrix)) {
        char delim = '\0';
        if (!args::get(path_delim).empty()) {
            delim = args::get(path_delim).at(0);
        }
        { // write the header
            stringstream header;
            if (delim) {
                header << "group.name" << "\t";
            }
            header << "path.name" << "\t"
                   << "path.length" << "\t"
                   << "node.count";
            graph.for_each_handle(
                [&](const handle_t& handle) {
                    header << "\t" << "node." << graph.get_id(handle);
                });
            std::cout << header.str() << std::endl;
        }
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string full_path_name = graph.get_path_name(p);
                std::string group_name = (delim ? split(full_path_name, delim)[0] : "");
                std::string path_name = (delim ? full_path_name.substr(full_path_name.find(delim)+1) : full_path_name);
                uint64_t path_length = 0;
                uint64_t path_step_count = 0;
                std::vector<bool> row(graph.get_node_count());
                graph.for_each_step_in_path(
                    p,
                    [&](const step_handle_t& s) {
                        const handle_t& h = graph.get_handle_of_step(s);
                        path_length += graph.get_length(h);
                        ++path_step_count;
                        row[graph.get_id(h)-1] = 1;
                    });
                if (delim) {
                    std::cout << group_name << "\t";
                }
                std::cout << path_name << "\t"
                          << path_length << "\t"
                          << path_step_count;
                for (uint64_t i = 0; i < row.size(); ++i) {
                    std::cout << "\t" << row[i];
                }
                std::cout << std::endl;
            });
    }

    if (args::get(distance_matrix)) {
        // TODO this breaks on a simple graph
        // WHY?!
        bool using_delim = !args::get(path_delim).empty();
        char delim = '\0';
        ska::flat_hash_map<std::string, uint64_t> path_group_ids;
        ska::flat_hash_map<path_handle_t, uint64_t> path_handle_group_ids;
        std::vector<std::string> path_groups;
        if (using_delim) {
            delim = args::get(path_delim).at(0);
            uint64_t i = 0;
            graph.for_each_path_handle(
                [&](const path_handle_t& p) {
                    std::string group_name = split(graph.get_path_name(p), delim)[0];
                    auto f = path_group_ids.find(group_name);
                    if (f == path_group_ids.end()) {
                        path_group_ids[group_name] = i++;
                        path_groups.push_back(group_name);
                    }
                    path_handle_group_ids[p] = path_group_ids[group_name];
                });
        }

        auto get_path_name
            = (using_delim ?
               (std::function<std::string(const uint64_t&)>)
               [&](const uint64_t& id) { return path_groups[id]; }
               :
               (std::function<std::string(const uint64_t&)>)
               [&](const uint64_t& id) { return graph.get_path_name(as_path_handle(id)); });

        auto get_path_id
            = (using_delim ?
               (std::function<uint64_t(const path_handle_t&)>)
               [&](const path_handle_t& p) {
                   return path_handle_group_ids[p];
               }
               :
               (std::function<uint64_t(const path_handle_t&)>)
               [&](const path_handle_t& p) {
                   return as_integer(p);
               });

        std::vector<uint64_t> bp_count;
        if (using_delim) {
            bp_count.resize(path_groups.size());
        } else {
            bp_count.resize(graph.get_path_count());
        }

        uint64_t path_max = 0;
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                path_max = std::max(path_max, as_integer(p));
            });

#pragma omp parallel for
        for (uint64_t i = 0; i < path_max; ++i) {
            path_handle_t p = as_path_handle(i + 1);
            uint64_t path_length = 0;
            graph.for_each_step_in_path(
                p,
                [&](const step_handle_t& s) {
                    path_length += graph.get_length(graph.get_handle_of_step(s));
                });
#pragma omp critical (bp_count)
            bp_count[get_path_id(p)] += path_length;
        }

        ska::flat_hash_map<std::pair<uint64_t, uint64_t>, uint64_t> path_intersection_length;
        graph.for_each_handle(
            [&](const handle_t& h) {
                uint64_t paths_here = 0;
                ska::flat_hash_map<uint64_t, uint64_t> local_path_lengths;
                size_t l = graph.get_length(h);
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& s) {
                        local_path_lengths[get_path_id(graph.get_path_handle_of_step(s))] += l;
                    });
#pragma omp critical (path_intersection_length)
                for (auto& p : local_path_lengths) {
                    for (auto& q : local_path_lengths) {
                        path_intersection_length[std::make_pair(p.first, q.first)] += std::min(p.second, q.second);
                    }
                }
            }, true);

        if (using_delim) {
            std::cout << "group.a" << "\t"
                      << "group.b" << "\t";
        } else {
            std::cout << "path.a" << "\t"
                      << "path.b" << "\t";
        }
        std::cout << "path.a.length" << "\t"
                  << "path.b.length" << "\t"
                  << "intersection" << "\t"
                  << "jaccard" << "\t"
                  << "euclidean"
                  << std::endl;
        for (auto& p : path_intersection_length) {
            auto& id_a = p.first.first;
            auto& id_b = p.first.second;
            auto& intersection = p.second;
            std::cout << get_path_name(id_a) << "\t"
                      << get_path_name(id_b) << "\t"
                      << bp_count[id_a] << "\t"
                      << bp_count[id_b] << "\t"
                      << intersection << "\t"
                      << (double)intersection / (double)(bp_count[id_a] + bp_count[id_b] - intersection) << "\t"
                      << std::sqrt((double)((bp_count[id_a] + bp_count[id_b] - intersection) - intersection)) << std::endl;
        }
    }

    if (!args::get(overlaps_file).empty()) {
        std::string line;
        ska::flat_hash_map<std::string, std::vector<std::string> > path_sets;
        ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;
        auto& x = args::get(overlaps_file);
        std::ifstream overlaps_in(x);
        while (std::getline(overlaps_in, line)) {
            // This file should contain tab-delimited lists of path names, one per line
            // the first field is an identifier
            std::vector<string> fields = split(line, '\t');
            assert(fields.size()==2);
            auto& name = fields[0];
            auto& pathset = path_sets[name];
            pathset.push_back(fields[1]);
            path_decomposition[fields[1]] = {};
        }

        std::vector<std::string> path_names;
        for (auto& p : path_decomposition) {
            path_decomposition[p.first] = {};
            path_names.push_back(p.first);
        }

        // the header
        std::cout << "group.name" << "\t"
                  << "query" << "\t"
                  << "target" << "\t"
                  << "overlap" << "\t"
                  << "overlap.frac" << std::endl;

#pragma omp parallel for
        for (uint64_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& decomposition = path_decomposition[path_name];
            // walk the path, adding each position to the decomposition
            path_handle_t path = graph.get_path_handle(path_name);
            uint64_t pos = 0;
            graph.for_each_step_in_path(path, [&](const step_handle_t& occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    nid_t id = graph.get_id(h);
                    uint64_t len = graph.get_length(h);
                    for (uint64_t i = 0; i < len; ++i) {
                        decomposition.push_back(make_pos_t(id, i, graph.get_is_reverse(h)));
                    }
                });
        }

        for (auto& s : path_sets) {
            auto& group_name = s.first;
            auto& paths = s.second;
            for (uint64_t i = 0; i < paths.size(); ++i) {
                for (uint64_t j = i+1; j < paths.size(); ++j) {
                    auto& v1 = path_decomposition[paths[i]];
                    auto& v2 = path_decomposition[paths[j]];
                    std::vector<pos_t> v3;
                    std::sort(v1.begin(), v1.end());
                    std::sort(v2.begin(), v2.end());
                    std::set_intersection(v1.begin(),v1.end(),
                                          v2.begin(),v2.end(),
                                          back_inserter(v3));
                    //ska::flat_hash_map<std::string, std::vector<pos_t> > path_decomposition;                    
                    std::cout << group_name << "\t" << paths[i] << "\t" << paths[j]
                              << "\t" << v3.size()
                              << "\t" << (float)v3.size()/((float)(v1.size()+v2.size())/2) << std::endl;
                }
            }
        }
    }

    return 0;
}

static Subcommand odgi_paths("paths", "Interrogate the embedded paths of a graph. Does not print anything to stdout by default!",
                              PIPELINE, 3, main_paths);


}
