#include "algorithms/untangle.hpp"
#include "args.hxx"
#include "odgi.hpp"
#include "subcommand.hpp"
#include "utils.hpp"
#include <omp.h>

namespace odgi {

using namespace odgi::subcommand;

int main_untangle(int argc, char **argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc - 1; ++i) {
        argv[i] = argv[i + 1];
    }
    const std::string prog_name = "odgi untangle";
    argv[0] = (char *)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("Project paths into reference-relative BEDPE, to decompose paralogy relationships.");
    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::ValueFlag<std::string> og_in_file(
        mandatory_opts, "FILE",
        "Load the succinct variation graph in ODGI format from this *FILE*. The file name usually "
        "ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format "
        "requires additional time!",
        {'i', "idx"});
    args::Group untangling_opts(parser, "[ Untangling Options ]");
    args::ValueFlag<std::string> _query_path(untangling_opts, "NAME", "Use this query path.",
                       {'q', "query-path"});
    args::ValueFlag<std::string> _target_path(untangling_opts, "NAME", "Use this target (reference) path.",
                       {'r', "target-path"});
    args::ValueFlag<std::string> _query_paths(untangling_opts, "FILE", "Use query paths listed (one per line) in FILE.",
                       {'Q', "query-paths"});
    args::ValueFlag<std::string> _target_paths(untangling_opts, "FILE", "Use target (reference) paths listed (one per line) in FILE.",
                       {'R', "target-paths"});
    args::ValueFlag<uint64_t> merge_dist(untangling_opts, "N", "Merge segments shorter than this length into previous segments.",
                                         {'m', "merge-dist"});
    args::Group debugging_opts(parser, "[ Debugging Options ]");
    args::Flag make_self_dotplot(debugging_opts, "DOTPLOT", "Render a table showing the positional dotplot of the query against itself.",
                                 {'s', "self-dotplot"});
    args::Group threading(parser, "[ Threading ]");
    args::ValueFlag<uint64_t> nthreads(
        threading, "N", "Number of threads to use for parallel operations.", {'t', "threads"});
    args::Group processing_info_opts(parser, "[ Processing Information ]");
    args::Flag progress(processing_info_opts, "progress", "Write the current progress to stderr.",
                        {'P', "progress"});
    args::Group program_info_opts(parser, "[ Program Information ]");
    args::HelpFlag help(program_info_opts, "help", "Print a help message for odgi untangle.",
                        {'h', "help"});

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

    if (!og_in_file) {
        std::cerr << "[odgi::untangle] error: please specify an input file from where to load the "
                     "graph via -i=[FILE], --idx=[FILE]."
                  << std::endl;
        return 1;
    }

    const uint64_t num_threads = args::get(nthreads) ? args::get(nthreads) : 1;

    graph_t graph;
    assert(argc > 0);
    {
        const std::string infile = args::get(og_in_file);
        if (!infile.empty()) {
            if (infile == "-") {
                graph.deserialize(std::cin);
            } else {
                utils::handle_gfa_odgi_input(infile, "untangle", args::get(progress), num_threads,
                                             graph);
            }
        }
    }    

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
                        std::cerr << "[odgi::untangle] error: in the path list there are duplicated path names."
                                  << std::endl;
                        exit(1);
                    }
                }
                ++num_of_paths_in_file;
            }
        }
        path_names_in.close();
        std::cerr << "[odgi::untangle] found " << paths.size() << "/" << num_of_paths_in_file
                  << " paths to consider." << std::endl;
        if (paths.empty()) {
            std::cerr << "[odgi::untangle] error: no path to consider." << std::endl;
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
            std::cerr << "[odgi::untangle] error: no path "
                      << path_name << " found in graph." << std::endl;
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
            std::cerr << "[odgi::untangle] error: no path "
                      << path_name << " found in graph." << std::endl;
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

    if (make_self_dotplot) {
        for (auto& query : query_paths) {
            algorithms::self_dotplot(graph, query);
        }
    } else {
        algorithms::untangle(graph,
                             query_paths,
                             target_paths,
                             args::get(merge_dist),
                             num_threads,
                             progress);
    }

    return 0;
}

static Subcommand odgi_untangle("untangle", "Project paths into reference-relative BEDPE, to decompose paralogy relationships.", PIPELINE, 3, main_untangle);

} // namespace odgi
