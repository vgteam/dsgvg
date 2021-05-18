#include "layout.hpp"
#include "draw.hpp"

namespace odgi {
namespace algorithms {
namespace layout {

using namespace handlegraph;

void to_tsv(std::ostream &out, const std::vector<double> &X, const std::vector<double> &Y, const std::vector<std::vector<handlegraph::handle_t>> weak_components) {
    /*std::vector<std::pair<double, uint64_t>> weak_component_order;
    for (uint64_t i = 0; i < weak_components.size(); ++i) {
        auto& weak_component = weak_components[i];
        uint64_t id_sum = 0;
        for (auto node_id : weak_component) {
            id_sum += node_id;
        }
        double avg_id = id_sum / (double) weak_component.size();
        weak_component_order.push_back(std::make_pair(avg_id, i));
    }
    std::sort(weak_component_order.begin(), weak_component_order.end());*/

    out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    out << "idx" << "\t" << "X" << "\t" << "Y" << "\t" << "component" << std::endl;

    uint64_t idx = 0;
    for (uint64_t num_component = 0; num_component < weak_components.size(); ++num_component) {
        for (auto& handle :  weak_components[num_component]) {
            uint64_t pos = 2 * number_bool_packing::unpack_number(handle);

            out << idx++ << "\t" << X[pos] << "\t" << Y[pos] << "\t" << num_component <<std::endl;
            out << idx++ << "\t" << X[pos + 1] << "\t" << Y[pos + 1] << "\t" << num_component<< std::endl;
        }
    }

}

double coord_dist(xy_d_t point1, xy_d_t point2) {
    double dx = point1.x - point2.x;
    double dy = point1.y - point2.y;
    return sqrt(dx * dx + dy * dy);
}

Layout::Layout(const std::vector<double> &X, const std::vector<double> &Y) {
    std::vector<uint64_t> vals;
    //assert(X.size)( == Y.size());
    for (auto& v : X) { min_value = std::min(v, min_value); }
    for (auto& v : Y) { min_value = std::min(v, min_value); }
    conv_t x, y;
    for (uint64_t i = 0; i < X.size(); ++i) {
        x.d = X[i] - min_value;
        y.d = Y[i] - min_value;
        vals.push_back(x.i);
        vals.push_back(y.i);
    }
    sdsl::util::assign(xy, sdsl::enc_vector<>(vals));
}

void Layout::serialize(std::ostream& out) {
    sdsl::write_member(min_value, out);
    xy.serialize(out);
}

void Layout::load(std::istream& in) {
    sdsl::read_member(min_value, in);
    xy.load(in);
}

void Layout::to_tsv(std::ostream &out) {
    out << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    out << "idx" << "\t" << "X" << "\t" << "Y" << std::endl;
    for (uint64_t i = 0; i < size(); ++i) {
        out << i << "\t" << get_x(i) << "\t" << get_y(i) << std::endl;
    }
}

xy_d_t Layout::coords(const handle_t& handle) {
    uint64_t idx = 2 * number_bool_packing::unpack_number(handle)
        + number_bool_packing::unpack_bit(handle);
    return { get_x(idx), get_y(idx) };
}

size_t Layout::size() {
    return xy.size()/2;
}

double Layout::get_x(uint64_t i) const {
    conv_t x;
    x.i = xy[2*i];
    return x.d + min_value;
}

double Layout::get_y(uint64_t i) const {
    conv_t y;
    y.i = xy[2*i+1];
    return y.d + min_value;
}

std::vector<double> Layout::get_X() {
    std::vector<double> X(size());
    for (uint64_t i = 0; i < size(); ++i) {
        X[i] = get_x(i);
    }
    return X;
}

std::vector<double> Layout::get_Y() {
    std::vector<double> Y(size());
    for (uint64_t i = 0; i < size(); ++i) {
        Y[i] = get_y(i);
    }
    return Y;
}

}
}
}
