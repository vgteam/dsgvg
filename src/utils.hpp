#include "odgi.hpp"
#include "gfa_to_handle.hpp"

#include <filesystem>

using namespace handlegraph;

#ifndef ODGI_UTILS_H
#define ODGI_UTILS_H

namespace utils {
    bool is_number(const std::string &s);
    void graph_deep_copy(const odgi::graph_t &source,
                         odgi::graph_t* target);
	bool ends_with(const std::string &fullString, const std::string &ending);
	int handle_gfa_odgi_input(const std::string infile, const std::string subcommmand_name, const bool progress,
							  const uint64_t num_threads, odgi::graph_t &graph);
	double median_of_sorted_vec(const std::vector<uint64_t>& vec);
}


#endif //ODGI_UTILS_H
