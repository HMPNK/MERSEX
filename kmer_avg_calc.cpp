#include <immintrin.h>
#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <sstream>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <algorithm>

// Configure batch size and number of threads
constexpr size_t BATCH_SIZE = 10000;
const unsigned int NUM_THREADS = std::thread::hardware_concurrency();


struct PartialResult {
    std::vector<double> sums;
    std::vector<uint64_t> counts;
    size_t max_field;  // max field index seen in this batch
};

// Thread-safe queue for batches
template<typename T>
class TSQueue {
    std::queue<T> q;
    std::mutex mtx;
    std::condition_variable cv;
    bool finished = false;
public:
    void push(T&& v) {
        std::unique_lock<std::mutex> lock(mtx);
        q.push(std::move(v));
        cv.notify_one();
    }
    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mtx);
        while (q.empty() && !finished) {
            cv.wait(lock);
        }
        if (q.empty() && finished) return false;
        out = std::move(q.front());
        q.pop();
        return true;
    }
    void set_finished() {
        std::unique_lock<std::mutex> lock(mtx);
        finished = true;
        cv.notify_all();
    }
};

TSQueue<std::vector<std::string>> batch_queue;

// Utility: split line by tabs
inline void split_line(const std::string& line, std::vector<std::string>& fields) {
    fields.clear();
    size_t start = 0;
    while (true) {
        size_t pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.emplace_back(line.substr(start));
            break;
        }
        fields.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
}

// SIMD helper: parse doubles from string views (non trivial; fallback to scalar)
// For this example, we parse scalar then load into SIMD for vector ops

// Worker function
PartialResult process_batch(const std::vector<std::string>& lines) {
    // We need to know max number of fields
    size_t max_fields = 0;

    // First pass: parse each line, store numeric fields 7..NF as doubles in a 2D vector for SIMD
    // To keep memory usage manageable, we store them row-wise.

    std::vector<std::vector<double>> data;
    data.reserve(lines.size());

    std::vector<std::string> fields;
    for (const auto& line : lines) {
        split_line(line, fields);
        if (fields.size() > max_fields) max_fields = fields.size();

        std::vector<double> numeric_fields;
        numeric_fields.reserve(fields.size() - 6); // fields from 7th onward

        for (size_t i = 6; i < fields.size(); i++) {
            // Parse double
            try {
                double val = std::stod(fields[i]);
                numeric_fields.push_back(val);
            } catch (...) {
                numeric_fields.push_back(0.0); // treat parse error as 0
            }
        }
        data.emplace_back(std::move(numeric_fields));
    }

    size_t num_fields = max_fields >= 7 ? max_fields - 6 : 0;

    std::vector<double> sum(num_fields, 0.0);
    std::vector<uint64_t> count(num_fields, 0);

    // Vectorized accumulation
    // AVX2 works with 4 doubles (256 bits)
    constexpr size_t simd_width = 4;

    for (const auto& row : data) {
        size_t N = row.size();
        size_t i = 0;

        // SIMD loop for chunks of 4 fields
        for (; i + simd_width <= N; i += simd_width) {
            __m256d vals = _mm256_loadu_pd(&row[i]);
            __m256d ones = _mm256_set1_pd(1.0);
            __m256d mask = _mm256_cmp_pd(vals, ones, _CMP_GT_OQ); // vals > 1

            // Create mask to count how many pass
            int mask_int = _mm256_movemask_pd(mask);

            if (mask_int) {
                // Extract masked values
                double vals_arr[simd_width];
                _mm256_storeu_pd(vals_arr, vals);

                for (int k = 0; k < simd_width; k++) {
                    if (mask_int & (1 << k)) {
                        sum[i + k] += vals_arr[k];
                        count[i + k]++;
                    }
                }
            }
        }

        // Scalar remainder
        for (; i < N; i++) {
            if (row[i] > 1.0) {
                sum[i] += row[i];
                count[i]++;
            }
        }
    }

    return {std::move(sum), std::move(count), max_fields};
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Launch worker threads
    std::vector<std::thread> workers;
    std::mutex result_mutex;
    std::vector<double> global_sum;
    std::vector<uint64_t> global_count;
    size_t global_max_fields = 0;

    std::atomic<bool> done_reading(false);

    // Worker thread function
    auto worker = [&]() {
        std::vector<std::string> batch;
        while (batch_queue.pop(batch)) {
            auto partial = process_batch(batch);

            // Combine partial results into global (protected by mutex)
            std::lock_guard<std::mutex> lock(result_mutex);

            if (partial.max_field > global_max_fields) {
                global_sum.resize(partial.max_field - 6, 0.0);
                global_count.resize(partial.max_field - 6, 0);
                global_max_fields = partial.max_field;
            }

            // Partial sums count vectors can be smaller if max_field smaller,
            // so only add where valid
            for (size_t i = 0; i < partial.sums.size(); i++) {
                global_sum[i] += partial.sums[i];
                global_count[i] += partial.counts[i];
            }
        }
    };

    for (int i = 0; i < NUM_THREADS; i++) {
        workers.emplace_back(worker);
    }

    // Read input lines and create batches
    std::vector<std::string> batch;
    batch.reserve(BATCH_SIZE);
    std::string line;
    while (std::getline(std::cin, line)) {
        batch.push_back(std::move(line));
        if (batch.size() >= BATCH_SIZE) {
            batch_queue.push(std::move(batch));
            batch.clear();
            batch.reserve(BATCH_SIZE);
        }
    }
    if (!batch.empty()) {
        batch_queue.push(std::move(batch));
    }

    // Signal no more batches
    batch_queue.set_finished();

    // Wait for workers to finish
    for (auto& t : workers) t.join();

    // Output averages for fields from 7 to max_fields
    for (size_t i = 0; i < global_sum.size(); i++) {
        if (global_count[i] > 0) {
            std::cout << (i + 7) << '\t' << (global_sum[i] / global_count[i]) << '\n';
        }
    }

    return 0;
}
