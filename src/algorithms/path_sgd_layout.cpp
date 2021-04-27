#include "path_sgd_layout.hpp"
#include "algorithms/layout.hpp"

namespace odgi {
    namespace algorithms {

        void path_linear_sgd_layout(const PathHandleGraph &graph,
                                    const xp::XP &path_index,
                                    const std::vector<path_handle_t> &path_sgd_use_paths,
                                    const uint64_t &iter_max,
                                    const uint64_t &iter_with_max_learning_rate,
                                    const uint64_t &min_term_updates,
                                    const double &delta,
                                    const double &eps,
                                    const double &eta_max,
                                    const double &theta,
                                    const uint64_t &space,
                                    const uint64_t &space_max,
                                    const uint64_t &space_quantization_step,
                                    const uint64_t &nthreads,
                                    const bool &progress,
                                    const bool &snapshot,
                                    const std::string &snapshot_prefix,
                                    std::vector<std::atomic<double>> &X,
                                    std::vector<std::atomic<double>> &Y) {
#ifdef debug_path_sgd
            std::cerr << "iter_max: " << iter_max << std::endl;
            std::cerr << "min_term_updates: " << min_term_updates << std::endl;
            std::cerr << "delta: " << delta << std::endl;
            std::cerr << "eps: " << eps << std::endl;
            std::cerr << "theta: " << theta << std::endl;
            std::cerr << "space: " << space << std::endl;
#endif

            uint64_t total_term_updates = iter_max * min_term_updates;
            std::unique_ptr<progress_meter::ProgressMeter> progress_meter;
            if (progress) {
                progress_meter = std::make_unique<progress_meter::ProgressMeter>(
                    total_term_updates, "[odgi::path_linear_sgd_layout] 2D path-guided SGD:");
            }
            using namespace std::chrono_literals; // for timing stuff
            uint64_t num_nodes = graph.get_node_count();
            // is a snapshot in progress?
            atomic<bool> snapshot_in_progress;
            snapshot_in_progress.store(false);
            // here we record which snapshots were already processed
            std::vector<atomic<bool>> snapshot_progress(iter_max);
            // we will produce one less snapshot compared to iterations
            snapshot_progress[0].store(true);
            // seed them with the graph order
            uint64_t len = 0;
            bool at_least_one_path_with_more_than_one_step = false;

            for (auto &path : path_sgd_use_paths) {
                if (path_index.get_path_step_count(path) > 1){
                    at_least_one_path_with_more_than_one_step = true;
                    break;
                }
            }

            if (at_least_one_path_with_more_than_one_step){
                double w_min = (double) 1.0 / (double) (eta_max);

#ifdef debug_path_sgd
                std::cerr << "w_min " << w_min << std::endl;
#endif
                double w_max = 1.0;
                // get our schedule
                std::vector<double> etas = path_linear_sgd_layout_schedule(w_min, w_max, iter_max,
                                                                           iter_with_max_learning_rate,
                                                                           eps);

                // cache zipf zetas for our full path space (heavy, but one-off)
                std::vector<double> zetas((space <= space_max ? space : space_max + (space - space_max) / space_quantization_step + 1)+1);
                uint64_t last_quantized_i = 0;
#pragma omp parallel for schedule(static,1)
                for (uint64_t i = 1; i < space+1; ++i) {
                    uint64_t quantized_i = i;
                    uint64_t compressed_space = i;
                    if (i > space_max){
                        quantized_i = space_max + (i - space_max) / space_quantization_step + 1;
                        compressed_space = space_max + ((i - space_max) / space_quantization_step) * space_quantization_step;
                    }

                    if (quantized_i != last_quantized_i){
                        dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, compressed_space, theta);
                        zetas[quantized_i] = z_p.zeta();

                        last_quantized_i = quantized_i;
                    }
                }

                // how many term updates we make
                std::atomic<uint64_t> term_updates;
                term_updates.store(0);
                // learning rate
                std::atomic<double> eta;
                eta.store(etas.front());
                // our max delta
                std::atomic<double> Delta_max;
                Delta_max.store(0);
                // should we keep working?
                std::atomic<bool> work_todo;
                work_todo.store(true);
                // approximately what iteration we're on
                uint64_t iteration = 0;
                // launch a thread to update the learning rate, count iterations, and decide when to stop
                auto checker_lambda =
                        [&]() {
                            while (work_todo.load()) {
                                if (term_updates.load() > min_term_updates) {
                                    if (snapshot) {
                                        if (snapshot_progress[iteration].load() || iteration == iter_max) {
                                            iteration++;
                                            if (iteration == iter_max) {
                                                snapshot_in_progress.store(false);
                                            } else {
                                                snapshot_in_progress.store(true);
                                            }
                                        } else {
                                            snapshot_in_progress.store(true);
                                            continue;
                                        }
                                    } else {
                                        iteration++;
                                        snapshot_in_progress.store(false);
                                    }
                                    if (iteration > iter_max) {
                                        work_todo.store(false);
                                    } else if (Delta_max.load() <= delta) { // nb: this will also break at 0
                                        if (progress) {
                                            std::cerr << "[odgi::path_linear_sgd_layout] delta_max: " << Delta_max.load()
                                                      << " <= delta: "
                                                      << delta << ". Threshold reached, therefore ending iterations."
                                                      << std::endl;
                                        }
                                        work_todo.store(false);
                                    } else {
                                        eta.store(etas[iteration]); // update our learning rate
                                        Delta_max.store(delta); // set our delta max to the threshold
                                    }
                                    term_updates.store(0);
                                }
                                std::this_thread::sleep_for(1ms);
                            }
                        };

                auto worker_lambda =
                        [&](uint64_t tid) {
                            // everyone tries to seed with their own random data
                            const std::uint64_t seed = 9399220 + tid;
                            XoshiroCpp::Xoshiro256Plus gen(seed); // a nice, fast PRNG
                            // some references to literal bitvectors in the path index hmmm
                            const sdsl::int_vector<> &nr_iv = path_index.get_nr_iv();
                            const sdsl::int_vector<> &npi_iv = path_index.get_npi_iv();
                            // we'll sample from all path steps
                            std::uniform_int_distribution<uint64_t> dis_step = std::uniform_int_distribution<uint64_t>(0, nr_iv.size() - 1);
                            std::uniform_int_distribution<uint64_t> flip(0, 1);
                            while (work_todo.load()) {
                                if (!snapshot_in_progress.load()) {
                                    // sample the first node from all the nodes in the graph
                                    // pick a random position from all paths
                                    uint64_t step_index = dis_step(gen);
#ifdef debug_sample_from_nodes
                                    std::cerr << "step_index: " << step_index << std::endl;
#endif
                                    uint64_t path_i = npi_iv[step_index];
                                    path_handle_t path = as_path_handle(path_i);
#ifdef debug_sample_from_nodes
                                    std::cerr << "path integer: " << path_i << std::endl;
#endif
                                    step_handle_t step_a, step_b;
                                    as_integers(step_a)[0] = path_i; // path index
                                    size_t s_rank = nr_iv[step_index] - 1; // step rank in path
                                    as_integers(step_a)[1] = s_rank;
#ifdef debug_sample_from_nodes
                                    std::cerr << "step rank in path: " << nr_iv[step_index]  << std::endl;
#endif
                                    size_t path_step_count = path_index.get_path_step_count(path);
                                    if (path_step_count == 1){
                                        continue;
                                    }

                                    if (flip(gen)) {
                                        if (s_rank > 0 && flip(gen) || s_rank == path_step_count-1) {
                                            // go backward
                                            uint64_t jump_space = std::min(space, s_rank);
                                            uint64_t space = jump_space;
                                            if (jump_space > space_max){
                                                space = space_max + (jump_space - space_max) / space_quantization_step + 1;
                                            }
                                            dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, jump_space, theta, zetas[space]);
                                            dirtyzipf::dirty_zipfian_int_distribution<uint64_t> z(z_p);
                                            uint64_t z_i = z(gen);
                                            //assert(z_i <= path_space);
                                            as_integers(step_b)[0] = as_integer(path);
                                            as_integers(step_b)[1] = s_rank - z_i;
                                        } else {
                                            // go forward
                                            uint64_t jump_space = std::min(space, path_step_count - s_rank - 1);
                                            uint64_t space = jump_space;
                                            if (jump_space > space_max){
                                                space = space_max + (jump_space - space_max) / space_quantization_step + 1;
                                            }
                                            dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type z_p(1, jump_space, theta, zetas[space]);
                                            dirtyzipf::dirty_zipfian_int_distribution<uint64_t> z(z_p);
                                            uint64_t z_i = z(gen);
                                            //assert(z_i <= path_space);
                                            as_integers(step_b)[0] = as_integer(path);
                                            as_integers(step_b)[1] = s_rank + z_i;
                                        }
                                    } else {
                                        // sample randomly across the path
                                        std::uniform_int_distribution<uint64_t> rando(0, graph.get_step_count(path)-1);
                                        as_integers(step_b)[0] = as_integer(path);
                                        as_integers(step_b)[1] = rando(gen);
                                    }


                                    // and the graph handles, which we need to record the update
                                    handle_t term_i = graph.get_handle_of_step(step_a);
                                    handle_t term_j = graph.get_handle_of_step(step_b);
                                    uint64_t term_i_length = graph.get_length(term_i);
                                    uint64_t term_j_length = graph.get_length(term_j);

                                    // adjust the positions to the node starts
                                    size_t pos_in_path_a = path_index.get_position_of_step(step_a);
                                    size_t pos_in_path_b = path_index.get_position_of_step(step_b);

                                    // determine which end we're working with for each node
                                    bool term_i_is_rev = graph.get_is_reverse(term_i);
                                    bool use_other_end_a = flip(gen); // 1 == +; 0 == -
                                    if (use_other_end_a) {
                                        pos_in_path_a += term_i_length;
                                        // flip back if we were already reversed
                                        use_other_end_a = !term_i_is_rev;
                                    } else {
                                        use_other_end_a = term_i_is_rev;
                                    }
                                    bool term_j_is_rev = graph.get_is_reverse(term_j);
                                    bool use_other_end_b = flip(gen); // 1 == +; 0 == -
                                    if (use_other_end_b) {
                                        pos_in_path_b += term_j_length;
                                        // flip back if we were already reversed
                                        use_other_end_b = !term_j_is_rev;
                                    } else {
                                        use_other_end_b = term_j_is_rev;
                                    }

#ifdef debug_path_sgd
                                    std::cerr << "1. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                                    // assert(pos_in_path_a < path_index.get_path_length(path));
                                    // assert(pos_in_path_b < path_index.get_path_length(path));
#ifdef debug_path_sgd
                                    std::cerr << "2. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                                    // establish the term distance
                                    double term_dist = std::abs(
                                            static_cast<double>(pos_in_path_a) - static_cast<double>(pos_in_path_b));

                                    if (term_dist == 0) {
                                        term_dist = 1e-9;
                                    }
#ifdef eval_path_sgd
                                    std::string path_name = path_index.get_path_name(path);
                                std::cerr << path_name << "\t" << pos_in_path_a << "\t" << pos_in_path_b << "\t" << term_dist << std::endl;
#endif
                                    // assert(term_dist == zipf_int);
#ifdef debug_path_sgd
                                    std::cerr << "term_dist: " << term_dist << std::endl;
#endif
                                    double term_weight = 1.0 / (double) term_dist;

                                    double w_ij = term_weight;
#ifdef debug_path_sgd
                                    std::cerr << "w_ij = " << w_ij << std::endl;
#endif
                                    double mu = eta.load() * w_ij;
                                    if (mu > 1) {
                                        mu = 1;
                                    }
                                    // actual distance in graph
                                    double d_ij = term_dist;
                                    // identities
                                    uint64_t i = number_bool_packing::unpack_number(term_i);
                                    uint64_t j = number_bool_packing::unpack_number(term_j);
#ifdef debug_path_sgd
                                    #pragma omp critical (cerr)
                                std::cerr << "nodes are " << graph.get_id(term_i) << " and " << graph.get_id(term_j) << std::endl;
#endif
                                    // distance == magnitude in our 2D situation
                                    uint64_t offset_i = 0;
                                    uint64_t offset_j = 0;
                                    if (use_other_end_a) {
                                        offset_i += 1;
                                    }
                                    if (use_other_end_b) {
                                        offset_j += 1;
                                    }
                                    double dx = X[2 * i + offset_i].load() - X[2 * j + offset_j].load();
                                    double dy = Y[2 * i + offset_i].load() - Y[2 * j + offset_j].load();
                                    if (dx == 0) {
                                        dx = 1e-9; // avoid nan
                                    }
#ifdef debug_path_sgd
                                    #pragma omp critical (cerr)
                                std::cerr << "distance is " << dx << " but should be " << d_ij << std::endl;
#endif
                                    //double mag = dx; //sqrt(dx*dx + dy*dy);
                                    double mag = sqrt(dx * dx + dy * dy);
#ifdef debug_path_sgd
                                    std::cerr << "mu " << mu << " mag " << mag << " d_ij " << d_ij << std::endl;
#endif
                                    // check distances for early stopping
                                    double Delta = mu * (mag - d_ij) / 2;
                                    // try until we succeed. risky.
                                    double Delta_abs = std::abs(Delta);
#ifdef debug_path_sgd
                                    #pragma omp critical (cerr)
                                std::cerr << "Delta_abs " << Delta_abs << std::endl;
#endif
                                    // todo use atomic compare and swap
                                    while (Delta_abs > Delta_max.load()) {
                                        Delta_max.store(Delta_abs);
                                    }
                                    // calculate update
                                    double r = Delta / mag;
                                    double r_x = r * dx;
                                    double r_y = r * dy;
#ifdef debug_path_sgd
                                    #pragma omp critical (cerr)
                                std::cerr << "r_x is " << r_x << std::endl;
#endif
                                    // update our positions (atomically)
#ifdef debug_path_sgd
                                    std::cerr << "before X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                                    X[2 * i + offset_i].store(X[2 * i + offset_i].load() - r_x);
                                    Y[2 * i + offset_i].store(Y[2 * i + offset_i].load() - r_y);
                                    X[2 * j + offset_j].store(X[2 * j + offset_j].load() + r_x);
                                    Y[2 * j + offset_j].store(Y[2 * j + offset_j].load() + r_y);
#ifdef debug_path_sgd
                                    std::cerr << "after X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                                    term_updates++; // atomic
                                    if (progress) {
                                        progress_meter->increment(1);
                                    }
                                }
                            }
                        };

                auto snapshot_lambda =
                        [&]() {
                            uint64_t iter = 0;
                            while (snapshot && work_todo.load()) {
                                if ((iter < iteration) && iteration != iter_max) {
                                    std::cerr << "[odgi::path_linear_sgd_layout] snapshot thread: Taking snapshot!" << std::endl;
                                    // drop out of atomic stuff... maybe not the best way to do this
                                    std::vector<double> X_iter(X.size());
                                    uint64_t i = 0;
                                    for (auto &x : X) {
                                        X_iter[i++] = x.load();
                                    }
                                    std::vector<double> Y_iter(Y.size());
                                    i = 0;
                                    for (auto &y : Y) {
                                        Y_iter[i++] = y.load();
                                    }
                                    algorithms::layout::Layout layout(X_iter, Y_iter);
                                    std::string local_snapshot_prefix = snapshot_prefix + std::to_string(iter + 1);
                                    ofstream snapshot_out(local_snapshot_prefix);
                                    // write out
                                    layout.serialize(snapshot_out);
                                    iter = iteration;
                                    snapshot_in_progress.store(false);
                                    snapshot_progress[iter].store(true);
                                }
                                std::this_thread::sleep_for(1ms);
                            }

                        };

                std::thread checker(checker_lambda);
                std::thread snapshot_thread(snapshot_lambda);

                std::vector<std::thread> workers;
                workers.reserve(nthreads);
                for (uint64_t t = 0; t < nthreads; ++t) {
                    workers.emplace_back(worker_lambda, t);
                }

                for (uint64_t t = 0; t < nthreads; ++t) {
                    workers[t].join();
                }

                snapshot_thread.join();

                checker.join();
            }

            if (progress) {
                progress_meter->finish();
            }

            // drop out of atomic stuff... maybe not the best way to do this
            std::vector<double> X_final(X.size());
            uint64_t i = 0;
            for (auto &x : X) {
                X_final[i++] = x.load();
            }
        }

        std::vector<double> path_linear_sgd_layout_schedule(const double &w_min,
                                                            const double &w_max,
                                                            const uint64_t &iter_max,
                                                            const uint64_t &iter_with_max_learning_rate,
                                                            const double &eps) {
#ifdef debug_schedule
            std::cerr << "w_min: " << w_min << std::endl;
            std::cerr << "w_max: " << w_max << std::endl;
            std::cerr << "iter_max: " << iter_max << std::endl;
            std::cerr << "eps: " << eps << std::endl;
#endif
            double eta_max = 1.0 / w_min;
            double eta_min = eps / w_max;
            double lambda = log(eta_max / eta_min) / ((double) iter_max - 1);
#ifdef debug_schedule
            std::cerr << "eta_max: " << eta_max << std::endl;
            std::cerr << "eta_min: " << eta_min << std::endl;
            std::cerr << "lambda: " << lambda << std::endl;
#endif
            // initialize step sizes
            std::vector<double> etas;
            etas.reserve(iter_max+1);
#ifdef debug_schedule
            std::cerr << "etas: ";
#endif
            for (int64_t t = 0; t <= iter_max; t++) {
                etas.push_back(eta_max * exp(-lambda * (abs(t - (int64_t) iter_with_max_learning_rate))));
#ifdef debug_schedule
                std::cerr << etas.back() << ", ";
#endif
            }
#ifdef debug_schedule
            std::cerr << std::endl;
#endif
            return etas;
        }
/*
        void deterministic_path_linear_sgd(const PathHandleGraph &graph,
                                           const xp::XP &path_index,
                                           const std::vector<path_handle_t> &path_sgd_use_paths,
                                           const uint64_t &iter_max,
                                           const uint64_t &iter_with_max_learning_rate,
                                           const uint64_t &min_term_updates,
                                           const double &delta,
                                           const double &eps,
                                           const double &eta_max,
                                           const double &theta,
                                           const uint64_t &space,
                                           const std::string &seeding_string,
                                           const bool &progress,
                                           const bool &snapshot,
                                           std::vector<std::vector<double>> &snapshots,
                                           std::vector<std::atomic<double>> &X,
                                           std::vector<std::atomic<double>> &Y) {
            using namespace std::chrono_literals; // for timing stuff
            uint64_t num_nodes = graph.get_node_count();
            // seed them with the graph order
            uint64_t len = 0;
            graph.for_each_handle(
                    [&X, &graph, &len](const handle_t &handle) {
                        // nb: we assume that the graph provides a compact handle set
                        X[number_bool_packing::unpack_number(handle)].store(len);
                        len += graph.get_length(handle);
                    });
            // the longest path length measured in nucleotides
            size_t longest_path_in_nucleotides = 0;
            // the total path length in nucleotides
            size_t total_path_len_in_nucleotides = 0;
            // here we store all path nucleotides lengths so we know later from which path we sampled our random position from
            IITree<uint64_t, path_handle_t> path_nucleotide_tree;
            // iterate over all relevant path_handles:
            //  1. build the interval tree
            //  2. find out the longest path in nucleotides and store this number size_t
            //  3. add the current path length to the total length
            for (auto &path : path_sgd_use_paths) {
#ifdef debug_path_sgd
                std::string path_name = graph.get_path_name(path);
                std::cerr << path_name << std::endl;
                std::cerr << as_integer(path) << std::endl;
#endif
                size_t path_len = path_index.get_path_length(path);
#ifdef debug_path_sgd
                std::cerr << path_name << " has length: " << path_len << std::endl;
#endif
                path_nucleotide_tree.add(total_path_len_in_nucleotides,
                                         total_path_len_in_nucleotides + path_len,
                                         path);

                if (path_len > longest_path_in_nucleotides) {
                    longest_path_in_nucleotides = path_len;
                }
                total_path_len_in_nucleotides += path_len;
            }
            path_nucleotide_tree.index();

            double w_min = (double) 1.0 / (double) (eta_max);
#ifdef debug_path_sgd
            std::cerr << "w_min " << w_min << std::endl;
#endif
            double w_max = 1.0;
            // get our schedule
            std::vector<double> etas = path_linear_sgd_layout_schedule(w_min, w_max, iter_max,
                                                                       iter_with_max_learning_rate,
                                                                       eps);
            // initialize Zipfian distribution so we only have to calculate zeta once
            dirtyzipf::dirty_zipfian_int_distribution<uint64_t>::param_type p(1, space, theta);
            dirtyzipf::dirty_zipfian_int_distribution<uint64_t> zipfian(p);
            // how many term updates we make
            std::atomic<uint64_t> term_updates;
            term_updates.store(0);
            // learning rate
            std::atomic<double> eta;
            eta.store(etas.front());
            // our max delta
            std::atomic<double> Delta_max;
            Delta_max.store(0);
            // seed with the given string
            std::seed_seq seed(seeding_string.begin(), seeding_string.end());
            std::mt19937 gen(seed);
            std::uniform_int_distribution<uint64_t> dis(1, num_nodes);
            if (sample_from_path_steps) {
                dis = std::uniform_int_distribution<uint64_t>(0, path_index.get_np_bv().size() - 1);
            }
            if (sample_from_paths) {
                dis = std::uniform_int_distribution<uint64_t>(0, total_path_len_in_nucleotides - 1);
            }
            std::uniform_int_distribution<uint64_t> flip(0, 1);
            const sdsl::bit_vector &np_bv = path_index.get_np_bv();
            const sdsl::int_vector<> &nr_iv = path_index.get_nr_iv();
            const sdsl::int_vector<> &npi_iv = path_index.get_npi_iv();
            auto &np_bv_select = path_index.get_np_bv_select();
            uint64_t hit_num_paths = 0;
            step_handle_t s_h;
            uint64_t node_index;
            uint64_t next_node_index;
            for (uint64_t iteration = 0; iteration < iter_max; iteration++) {
                if (snapshot && iteration < iter_max - 1) {
                    // drop out of atomic stuff... maybe not the best way to do this
                    std::vector<double> X_iter(X.size());
                    uint64_t i = 0;
                    for (auto &x : X) {
                        X_iter[i++] = x.load();
                    }
                    snapshots.push_back(X_iter);
                }
                for (uint64_t term_update = 0; term_update < min_term_updates; term_update++) {
                    // pick a random position from all paths
                    uint64_t pos = dis(gen);
                    size_t pos_in_path_a;
                    size_t path_len;
                    path_handle_t path;
#ifdef debug_path_sgd
                    std::cerr << "uniform_position: " << pos << std::endl;
#endif
                    if (sample_from_paths) {
                        // use our interval tree to get the path handle and path nucleotide position of the picked position
                        //std::vector<Interval<size_t, path_handle_t> > result;
                        //result = path_nucleotide_tree.findOverlapping(pos, pos);
                        std::vector<size_t> a;
                        path_nucleotide_tree.overlap(pos, pos + 1, a);
                        if (a.empty()) {
                            std::cerr << "[odgi::path_sgd] no overlapping intervals at position " << pos
                                      << std::endl;
                            exit(1);
                        }
                        auto &p = a[0];
                        path = path_nucleotide_tree.data(p);
                        size_t path_start_pos = path_nucleotide_tree.start(p);
                        // size_t path_end_pos = result[0].stop;
                        path_len = path_index.get_path_length(path) - 1;
                        // we have a 0-based positioning in the path index
                        pos_in_path_a = pos - path_start_pos;
                    } else {
                        if (sample_from_path_steps) {
                            // did we hit a node and not a path?
                            if (np_bv[pos] == 1) {
                                continue;
                            } else {
                                uint64_t path_i = npi_iv[pos];
                                path = as_path_handle(path_i);
                                as_integers(s_h)[0] = path_i; // path index
                                as_integers(s_h)[1] = nr_iv[pos] - 1; // maybe -1?! // step rank in path
                                pos_in_path_a = path_index.get_position_of_step(s_h);
                                path_len = path_index.get_path_length(path) - 1;
#ifdef debug_sample_from_nodes
                                std::cerr << "path_len: " << path_len << std::endl;
                                std::cerr << "path id: " << (npi_iv[pos]) << std::endl;
                                std::cerr << "step rank in path: " << (nr_iv[pos] - 1) << std::endl;
                                std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
#endif
                            }
                            // default: sample the first node from all the nodes in the graph
                        } else {
                            node_index = np_bv_select(pos);
                            // did we hit the last node?
                            if (pos == num_nodes) {
                                next_node_index = np_bv.size();
                            } else {
                                next_node_index = np_bv_select(pos + 1);
                            }
                            hit_num_paths = next_node_index - node_index - 1;
                            if (hit_num_paths == 0) {
                                continue;
                            }
                            std::uniform_int_distribution<uint64_t> dis_path(1, hit_num_paths);
                            uint64_t path_pos_in_np_iv = dis_path(gen);
#ifdef debug_sample_from_nodes
                            std::cerr << "path_pos_in_np_iv first: " << path_pos_in_np_iv << std::endl;
                            std::cerr << "node_index: " << node_index << std::endl;
#endif
                            path_pos_in_np_iv = node_index + path_pos_in_np_iv;
#ifdef debug_sample_from_nodes
                            std::cerr << "path pos in np_iv: " << path_pos_in_np_iv << std::endl;
#endif
                            uint64_t path_i = npi_iv[path_pos_in_np_iv];
                            path = as_path_handle(path_i);
#ifdef debug_sample_from_nodes
                            std::cerr << "path integer: " << path_i << std::endl;
#endif
                            as_integers(s_h)[0] = path_i; // path index
                            as_integers(s_h)[1] = nr_iv[path_pos_in_np_iv] - 1; // step rank in path
#ifdef debug_sample_from_nodes
                            std::cerr << "step rank in path: " << nr_iv[path_pos_in_np_iv]  << std::endl;
#endif
                            pos_in_path_a = path_index.get_position_of_step(s_h);
#ifdef debug_sample_from_nodes
                            std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
#endif
                            path_len = path_index.get_path_length(path) - 1;
#ifdef debug_sample_from_nodes
                            std::cerr << "path_len " << path_len << std::endl;
                            std::cerr << "node count " << num_nodes << std::endl;
#endif
                        }
                    }
                    uint64_t zipf_int = zipfian(gen);
#ifdef debug_path_sgd
                    std::cerr << "random pos: " << pos << std::endl;
                    std::cerr << "path_start_pos: " << path_start_pos << std::endl;
                    std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                    std::cerr << "path_len: " << path_len << std::endl;
                    //std::cerr << "zipf: " << zipf_int << std::endl;
#endif
                    size_t pos_in_path_b = pos_in_path_a;
                    if (flip(gen)) {
                        if (zipf_int > pos_in_path_a) {
                            if (pos_in_path_a == 0) {
                                continue;
                            } else {
                                zipf_int %= pos_in_path_a;
                            }
                        }
                        pos_in_path_b -= zipf_int;
                    } else {
                        if (zipf_int > path_len - pos_in_path_a) {
                            if (path_len - pos_in_path_a == 0) {
                                continue;
                            } else {
                                zipf_int %= path_len - pos_in_path_a;
                            }
                        }
                        pos_in_path_b += zipf_int;
                    }
#ifdef debug_path_sgd
                    std::cerr << "zipf: " << zipf_int << std::endl;
                    std::cerr << "pos_in_path_a: " << pos_in_path_a << std::endl;
                    std::cerr << "pos_in_path_b: " << pos_in_path_b << std::endl;
#endif
                    // get the step handles
                    step_handle_t step_a = path_index.get_step_at_position(path, pos_in_path_a);
                    step_handle_t step_b = path_index.get_step_at_position(path, pos_in_path_b);

                    // and the graph handles, which we need to record the update
                    handle_t term_i = path_index.get_handle_of_step(step_a);
                    handle_t term_j = path_index.get_handle_of_step(step_b);

                    // adjust the positions to the node starts
                    pos_in_path_a = path_index.get_position_of_step(step_a);
                    pos_in_path_b = path_index.get_position_of_step(step_b);
#ifdef debug_path_sgd
                    std::cerr << "1. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                    // assert(pos_in_path_a < path_index.get_path_length(path));
                    // assert(pos_in_path_b < path_index.get_path_length(path));

                    // and adjust to account for our relative orientation
                    if (graph.get_is_reverse(term_i)) {
                        pos_in_path_a += graph.get_length(term_i);
                    }
                    if (graph.get_is_reverse(term_j)) {
                        pos_in_path_b += graph.get_length(term_j);
                    }
#ifdef debug_path_sgd
                    std::cerr << "2. pos in path " << pos_in_path_a << " " << pos_in_path_b << std::endl;
#endif
                    // establish the term distance
                    double term_dist = std::abs(
                            static_cast<double>(pos_in_path_a) - static_cast<double>(pos_in_path_b));

                    if (term_dist == 0) {
                        continue;
                        // term_dist = 1e-9;
                    }
#ifdef eval_path_sgd
                    std::string path_name = path_index.get_path_name(path);
                    std::cerr << path_name << "\t" << pos_in_path_a << "\t" << pos_in_path_b << "\t" << term_dist << std::endl;
#endif
                    // assert(term_dist == zipf_int);
#ifdef debug_path_sgd
                    std::cerr << "term_dist: " << term_dist << std::endl;
#endif
                    double term_weight = 1.0 / term_dist;

                    double w_ij = term_weight;
#ifdef debug_path_sgd
                    std::cerr << "w_ij = " << w_ij << std::endl;
#endif
                    double mu = eta.load() * w_ij;
                    if (mu > 1) {
                        mu = 1;
                    }
                    // actual distance in graph
                    double d_ij = term_dist;
                    // identities
                    uint64_t i = number_bool_packing::unpack_number(term_i);
                    uint64_t j = number_bool_packing::unpack_number(term_j);
#ifdef debug_path_sgd
#pragma omp critical (cerr)
                    std::cerr << "nodes are " << graph.get_id(term_i) << " and " << graph.get_id(term_j) << std::endl;
#endif
                    // distance == magnitude in our 1D situation
                    double dx = X[i].load() - X[j].load();
                    if (dx == 0) {
                        dx = 1e-9; // avoid nan
                    }
#ifdef debug_path_sgd
#pragma omp critical (cerr)
                    std::cerr << "distance is " << dx << " but should be " << d_ij << std::endl;
#endif
                    //double mag = dx; //sqrt(dx*dx + dy*dy);
                    double mag = std::abs(dx);
#ifdef debug_path_sgd
                    std::cerr << "mu " << mu << " mag " << mag << " d_ij " << d_ij << std::endl;
#endif
                    // check distances for early stopping
                    double Delta = mu * (mag - d_ij) / 2;
                    // try until we succeed. risky.
                    double Delta_abs = std::abs(Delta);
#ifdef debug_path_sgd
#pragma omp critical (cerr)
                    std::cerr << "Delta_abs " << Delta_abs << std::endl;
#endif
                    while (Delta_abs > Delta_max.load()) {
                        Delta_max.store(Delta_abs);
                    }
                    // calculate update
                    double r = Delta / mag;
                    double r_x = r * dx;
#ifdef debug_path_sgd
#pragma omp critical (cerr)
                    std::cerr << "r_x is " << r_x << std::endl;
#endif
                    // update our positions (atomically)
#ifdef debug_path_sgd
                    std::cerr << "before X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                    X[i].store(X[i].load() - r_x);
                    X[j].store(X[j].load() + r_x);
#ifdef debug_path_sgd
                    std::cerr << "after X[i] " << X[i].load() << " X[j] " << X[j].load() << std::endl;
#endif
                    term_updates.store(term_updates.load() + 1);
                }
                if (Delta_max.load() <= delta) {
                    if (progress) {
                        std::cerr << "[odgi::path_linear_sgd_layout] delta_max: " << Delta_max.load() << " <= delta: " << delta
                                  << ". Threshold reached, therefore ending iterations." << std::endl;
                    }
                    break;
                } else {
                    if (progress) {
                        double percent_progress = ((double) (iteration + 1) / (double) iter_max) * 100.0;
                        std::cerr << std::fixed << std::setprecision(2) << "[odgi::path_linear_sgd_layout]: " << percent_progress
                                  << "% progress: "
                                     "iteration: " << (iteration + 1) <<
                                  ", eta: " << eta.load() <<
                                  ", delta: " << Delta_max.load() <<
                                  ", number of updates: " << term_updates.load() << std::endl;
                    }

                    // If it is the last iteration, there is no need to update the next values, and it is avoided
                    // to request an element outside the vector
                    if (iteration + 1 < iter_max) {
                        eta.store(etas[iteration + 1]); // update our learning rate
                        Delta_max.store(delta); // set our delta max to the threshold
                    }
                }
                term_updates.store(0);
            }

            // drop out of atomic stuff... maybe not the best way to do this
            std::vector<double> X_final(X.size());
            uint64_t i = 0;
            for (auto &x : X) {
                X_final[i++] = x.load();
            }
        }*/
    }
}
