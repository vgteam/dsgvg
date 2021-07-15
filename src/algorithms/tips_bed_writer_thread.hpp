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

		class tips_bed_writer {

		private:
			// BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
			struct tips_bed_record_t {
				std::string chrom;
				uint64_t chromStart;
				uint64_t chromEnd;
				double query_pos_median;
				std::string path;
				uint64_t path_pos;
				bool walking_dir;
			};

			std::thread writer_thread;
			atomic_queue::AtomicQueue2<tips_bed_record_t*, 2 << 16> bed_record_queue;
			std::atomic<bool> work_todo;

		public:

			tips_bed_writer(void) {
				work_todo.store(false);
			}

			~tips_bed_writer(void) {
				close_writer();
			}

			void writer_func(void) {
				auto* bed_record = new tips_bed_record_t();
				while (work_todo.load() || !bed_record_queue.was_empty()) {
					if (bed_record_queue.try_pop(bed_record)) {
						do {
							// BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
							std::cout.precision(dbl::max_digits10);
							std::cout << std::fixed
									  << bed_record->chrom << "\t" // chrom
									  << bed_record->chromStart << "\t" // chromStart
									  << bed_record->chromEnd << "\t" // chromEnd
									  << bed_record->query_pos_median << "\t"
									  << bed_record->path << "\t"
									  << bed_record->path_pos << "\t"
									  << bed_record->walking_dir
									  << std::endl;
						} while (bed_record_queue.try_pop(bed_record));
					} else {
						std::this_thread::sleep_for(std::chrono::nanoseconds(1));
					}
				}
			}

			// open backing file and start writer_thread
			void open_writer(void) {
				if (!work_todo.load()) {
					work_todo.store(true);
					writer_thread = std::thread(&tips_bed_writer::writer_func, this);
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
						const double &query_pos_median, const std::string &path, const uint64_t &path_pos, const bool walking_dir) {
				bed_record_queue.push(new tips_bed_record_t{chrom, chromStart, chromEnd, query_pos_median,
												path, path_pos, walking_dir});
			}
		};

	}

}