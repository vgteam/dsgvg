#include <iostream>
#include <sdsl/enc_vector.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include "atomic_image.hpp"
#include "weakly_connected_components.hpp"

namespace odgi {

namespace algorithms {

namespace layout {

using namespace handlegraph;

void to_tsv(std::ostream &out,
            const std::vector<double> &X,
            const std::vector<double> &Y,
            const std::vector<std::vector<handlegraph::handle_t>> weak_components);

double coord_dist(const xy_d_t, const xy_d_t);

union conv_t { uint64_t i; double d; };

class Layout {
    sdsl::enc_vector<> xy;
    double min_value = std::numeric_limits<double>::max();
public:
    Layout() { }
    Layout(const std::vector<double> &X, const std::vector<double> &Y);
    void serialize(std::ostream& out);
    void load(std::istream& in);
    void to_tsv(std::ostream &out);
    xy_d_t coords(const handle_t& handle);
    size_t size();
    double get_x(uint64_t i) const;
    double get_y(uint64_t i) const;
    std::vector<double> get_X();
    std::vector<double> get_Y();
};

}

}

}
