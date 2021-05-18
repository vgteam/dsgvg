#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <thread>
#include "atomic_queue.h"
#include <limits>

typedef std::numeric_limits< double > dbl;

namespace odgi {
namespace algorithms {

class bed_records_class {

private:
    // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
    struct bed_record_t {
        std::string chrom;
        uint64_t chromStart;
        uint64_t chromEnd;
        double path_layout_dist;
        uint64_t path_nuc_dist;
        double path_layout_nuc_dist_ratio;
    };

    std::thread writer_thread;
    atomic_queue::AtomicQueue2<bed_record_t*, 2 << 16> bed_record_queue;
    std::atomic<bool> work_todo;

public:

    bed_records_class(void) {
        work_todo.store(false);
    }

    ~bed_records_class(void) {
        close_writer();
    }

    void writer_func(void) {
        /*
        std::ofstream writer(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
        */
        auto* bed_record = new bed_record_t();
        while (work_todo.load() || !bed_record_queue.was_empty()) {
            if (bed_record_queue.try_pop(bed_record)) {
                do {
                    // writer.write((char*)&ival, sizeof(Interval));
                    // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
                    std::cout.precision(dbl::max_digits10);
                    std::cout << std::fixed
                              << bed_record->chrom << "\t" // chrom
                              << bed_record->chromStart << "\t" // chromStart
                              << bed_record->chromEnd << "\t" // chromEnd
                              << bed_record->path_layout_dist << "\t"
                              << bed_record->path_nuc_dist << "\t"
                              << bed_record->path_layout_nuc_dist_ratio
                              << std::endl;
                } while (bed_record_queue.try_pop(bed_record));
            } else {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
        }
        // writer.close();
    }

    // open backing file and start writer_thread
    void open_writer(void) {
        if (!work_todo.load()) {
            work_todo.store(true);
            writer_thread = std::thread(&bed_records_class::writer_func, this);
        }
    }

    void close_writer(void) {
        if (work_todo.load()) {
            work_todo.store(false);
            while (!bed_record_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            if (writer_thread.joinable()) {
                writer_thread.join();
            }
        }
    }

    /// write into our write buffer
    /// open_writer() must be called first to set up our buffer and writer
	void append(const std::string &chrom, const uint64_t &chromStart, const uint64_t &chromEnd,
          const double &path_layout_dist, const uint64_t &path_nuc_dist, const double &path_layout_nuc_dist_ratio) {
        bed_record_queue.push(new bed_record_t{chrom, chromStart, chromEnd, path_layout_dist, path_nuc_dist, path_layout_nuc_dist_ratio});
    }
};

}

}