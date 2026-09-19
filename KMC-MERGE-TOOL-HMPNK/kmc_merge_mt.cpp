// Build:
//   g++ -std=c++17 -O3 -march=native -flto -pthread kmc_merge_mt.cpp kmc_api/*.cpp -o kmc_merge_mt
//
// Usage:
//   kmc_merge_mt [options] db1 db2 ...
//
//   -o FILE   write to FILE instead of stdout
//   -z        compress with pigz (implied if FILE ends in .gz)
//   -p N      pigz threads (default: cores - 1)
//   -l N      pigz level 1-9 (default: 4)
//   -B KiB    pigz block size in KiB (default: 1024)
//   -c "CMD"  custom compressor command writing to stdout (overrides pigz options)
//   -t N      decoder threads (default: 3)
//   -m CPU    pin the merge thread to logical CPU number CPU
//
// Example (i7-5820K: merge on CPU 0, sibling CPU 6 left idle):
//   ./kmc_merge_mt -o merged.tsv.gz -m 0 -t 3 \
//     -c "taskset -c 1-5,7-11 pigz -c -p 8 -4 -b 1024 -i" db*

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <pthread.h>
#include <sched.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

// Expose CKmerAPI internals so the first 64-bit word can be used as a sort key.
// All standard headers must be included ABOVE this point.
#define protected public
#define private public
#include "kmc_file.h"
#undef protected
#undef private

static inline uint64_t key_of(const CKmerAPI& km) { return km.kmer_data[0]; }

// ============================================================
// Output sink: stdout, a file, or a compressor process
// ============================================================
class Sink {
public:
    Sink() = default;
    ~Sink() { try { close(); } catch (...) {} }
    Sink(const Sink&) = delete;
    Sink& operator=(const Sink&) = delete;

    void open(const std::string& path, bool compress, const std::vector<std::string>& cmd)
    {
        if (!compress) {
            if (path.empty()) {
                f_ = stdout;
            } else {
                f_ = std::fopen(path.c_str(), "wb");
                if (!f_) throw std::runtime_error("Cannot open output file: " + path);
            }
            return;
        }

        int fds[2];
        if (pipe(fds) != 0) throw std::runtime_error("pipe() failed");
#ifdef F_SETPIPE_SZ
        fcntl(fds[1], F_SETPIPE_SZ, 1 << 20);   // best effort
#endif

        std::vector<const char*> args;
        for (const auto& s : cmd) args.push_back(s.c_str());
        args.push_back(nullptr);

        std::fflush(stdout);
        std::fflush(stderr);

        pid_ = fork();
        if (pid_ < 0) throw std::runtime_error("fork() failed");

        if (pid_ == 0) {                        // child -> compressor
            ::close(fds[1]);
            dup2(fds[0], STDIN_FILENO);
            ::close(fds[0]);

            if (!path.empty()) {
                int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
                if (fd < 0) { std::perror("cannot open output file"); _exit(126); }
                dup2(fd, STDOUT_FILENO);
                ::close(fd);
            }
            execvp(args[0], const_cast<char* const*>(args.data()));
            std::perror("cannot execute compressor (is it installed and on PATH?)");
            _exit(127);
        }

        ::close(fds[0]);
        f_ = fdopen(fds[1], "wb");
        if (!f_) throw std::runtime_error("fdopen() failed");
    }

    void write(const char* data, size_t n)
    {
        if (n && std::fwrite(data, 1, n, f_) != n)
            throw std::runtime_error("Write to output failed (compressor exited early or disk full?)");
    }

    void close()
    {
        if (f_) {
            if (f_ == stdout) std::fflush(stdout);
            else              std::fclose(f_);
            f_ = nullptr;
        }
        if (pid_ > 0) {
            int status = 0;
            const pid_t p = pid_;
            pid_ = -1;
            waitpid(p, &status, 0);
            if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
                throw std::runtime_error("compressor process failed");
        }
    }

private:
    FILE* f_   = nullptr;
    pid_t pid_ = -1;
};

// ============================================================
// Async writer: merge thread fills buffers, writer thread drains them
// ============================================================
class AsyncWriter {
public:
    struct Buf {
        std::vector<char> data;
        size_t len = 0;
        explicit Buf(size_t n) : data(n) {}
    };

    AsyncWriter(Sink& sink, size_t buf_size, size_t nbufs = 8) : sink_(sink)
    {
        for (size_t i = 0; i < nbufs; ++i) {
            pool_.emplace_back(buf_size);
            free_.push_back(&pool_.back());
        }
        thr_ = std::thread([this] { run(); });
    }
    ~AsyncWriter() { stop(); }

    Buf* acquire()
    {
        std::unique_lock<std::mutex> lk(m_);
        cv_free_.wait(lk, [&] { return failed_ || !free_.empty(); });
        if (failed_) std::rethrow_exception(err_);
        Buf* b = free_.front();
        free_.pop_front();
        return b;
    }

    void submit(Buf* b, size_t len)
    {
        b->len = len;
        { std::lock_guard<std::mutex> lk(m_); full_.push_back(b); }
        cv_full_.notify_one();
    }

    void finish()
    {
        stop();
        if (failed_) std::rethrow_exception(err_);
    }

private:
    void stop()
    {
        { std::lock_guard<std::mutex> lk(m_); done_ = true; }
        cv_full_.notify_all();
        if (thr_.joinable()) thr_.join();
    }

    void run()
    {
        for (;;) {
            Buf* b;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_full_.wait(lk, [&] { return done_ || !full_.empty(); });
                if (full_.empty()) return;
                b = full_.front();
                full_.pop_front();
            }
            try {
                sink_.write(b->data.data(), b->len);
            } catch (...) {
                std::lock_guard<std::mutex> lk(m_);
                err_ = std::current_exception();
                failed_ = true;
                cv_free_.notify_all();
                return;
            }
            { std::lock_guard<std::mutex> lk(m_); free_.push_back(b); }
            cv_free_.notify_one();
        }
    }

    Sink& sink_;
    std::deque<Buf> pool_;
    std::deque<Buf*> free_, full_;
    std::mutex m_;
    std::condition_variable cv_free_, cv_full_;
    std::exception_ptr err_;
    bool failed_ = false, done_ = false;
    std::thread thr_;
};

// ============================================================
// Decoder pool: a few worker threads decode all databases in batches
// ============================================================
class DecoderPool {
public:
    static constexpr size_t BATCH  = 2048;   // k-mers per batch
    static constexpr size_t NBATCH = 4;      // batches in flight per database

    DecoderPool(const std::vector<std::string>& paths, int nthreads)
        : nthreads_(std::max(1, nthreads))
    {
        streams_.reserve(paths.size());
        for (const auto& path : paths) {
            std::cerr << "  " << path << std::endl;
            auto s = std::make_unique<Stream>();
            if (!s->db.OpenForListing(path))
                throw std::runtime_error("Cannot open KMC database: " + path);

            uint32_t kl, mode, cs, lut, sig, minc;
            uint64_t maxc, total;
            s->db.Info(kl, mode, cs, lut, sig, minc, maxc, total);

            if (streams_.empty()) k_ = kl;
            else if (kl != k_) throw std::runtime_error("Databases have different k");

            for (size_t j = 0; j < NBATCH; ++j) {
                s->pool.emplace_back(std::make_unique<Batch>(BATCH, kl));
                s->free.push_back(s->pool.back().get());
            }
            streams_.push_back(std::move(s));
        }
    }

    ~DecoderPool() { stop(); }

    DecoderPool(const DecoderPool&) = delete;
    DecoderPool& operator=(const DecoderPool&) = delete;

    // Start decoding. Call before advance(). Call BEFORE pinning the calling thread.
    void start()
    {
        {
            std::lock_guard<std::mutex> lk(m_);
            for (auto& s : streams_) { s->queued = true; work_.push_back(s.get()); }
        }
        for (int t = 0; t < nthreads_; ++t)
            workers_.emplace_back([this] { worker(); });
    }

    size_t   size() const { return streams_.size(); }
    uint32_t k() const { return k_; }
    uint64_t stalls() const { return stalls_; }

    // ---- consumer side (merge thread only) ----

    // Move stream i to its next k-mer. Returns false when exhausted.
    bool advance(size_t i)
    {
        Stream& s = *streams_[i];
        if (s.done) return false;
        if (s.cur && ++s.pos < s.cur->size) return true;

        if (s.cur) {
            if (s.cur->eof) { s.done = true; return false; }
            release(s, s.cur);
            s.cur = nullptr;
        }

        Batch* b = wait_full(s);
        s.cur = b;
        s.pos = 0;
        if (b->size == 0) { s.done = true; return false; }
        return true;
    }

    uint64_t         key(size_t i)   const { const Stream& s = *streams_[i]; return s.cur->keys[s.pos]; }
    uint64_t         count(size_t i) const { const Stream& s = *streams_[i]; return s.cur->counts[s.pos]; }
    const CKmerAPI&  kmer(size_t i)  const { const Stream& s = *streams_[i]; return s.cur->kmers[s.pos]; }

private:
    struct Batch {
        std::vector<CKmerAPI> kmers;
        std::vector<uint64_t> keys, counts;
        size_t size = 0;
        bool   eof  = false;

        Batch(size_t cap, uint32_t k) : keys(cap), counts(cap)
        {
            kmers.reserve(cap);                    // reserved: no reallocation, no copies
            for (size_t i = 0; i < cap; ++i) kmers.emplace_back(k);
        }
    };

    struct Stream {
        CKMCFile db;
        std::vector<std::unique_ptr<Batch>> pool;

        // guarded by m_
        std::deque<Batch*> free, full;
        bool queued = false, busy = false, eof = false;

        // consumer only
        Batch* cur = nullptr;
        size_t pos = 0;
        bool   done = false;

        ~Stream() { db.Close(); }
    };

    void release(Stream& s, Batch* b)
    {
        bool wake = false;
        {
            std::lock_guard<std::mutex> lk(m_);
            s.free.push_back(b);
            if (!s.eof && !s.queued && !s.busy) {
                s.queued = true;
                work_.push_back(&s);
                wake = true;
            }
        }
        if (wake) cv_work_.notify_one();
    }

    Batch* wait_full(Stream& s)
    {
        std::unique_lock<std::mutex> lk(m_);
        if (s.full.empty()) ++stalls_;             // merge thread had to wait for decoding
        cv_ready_.wait(lk, [&] { return !s.full.empty(); });
        Batch* b = s.full.front();
        s.full.pop_front();
        return b;
    }

    void worker()
    {
        for (;;) {
            Stream* s;
            Batch*  b;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_work_.wait(lk, [&] { return stop_ || !work_.empty(); });
                if (stop_) return;
                s = work_.front();
                work_.pop_front();
                s->queued = false;
                s->busy   = true;
                b = s->free.front();               // queued streams always have a free batch
                s->free.pop_front();
            }

            size_t n = 0;
            uint64_t c;
            while (n < BATCH && s->db.ReadNextKmer(b->kmers[n], c)) {
                b->keys[n]   = key_of(b->kmers[n]);
                b->counts[n] = c;
                ++n;
            }
            b->size = n;
            b->eof  = (n < BATCH);

            bool requeued = false;
            {
                std::lock_guard<std::mutex> lk(m_);
                s->busy = false;
                s->full.push_back(b);
                if (b->eof) {
                    s->eof = true;
                } else if (!s->free.empty()) {
                    s->queued = true;
                    work_.push_back(s);
                    requeued = true;
                }
            }
            cv_ready_.notify_one();
            if (requeued) cv_work_.notify_one();
        }
    }

    void stop()
    {
        { std::lock_guard<std::mutex> lk(m_); stop_ = true; }
        cv_work_.notify_all();
        for (auto& t : workers_) if (t.joinable()) t.join();
        workers_.clear();
    }

    int nthreads_;
    uint32_t k_ = 0;
    uint64_t stalls_ = 0;                          // merge thread only

    std::vector<std::unique_ptr<Stream>> streams_;
    std::deque<Stream*> work_;
    std::mutex m_;
    std::condition_variable cv_work_, cv_ready_;
    bool stop_ = false;
    std::vector<std::thread> workers_;
};

// ============================================================
// SIMD helpers (AVX2, with scalar fallback)
// ============================================================
#ifdef __AVX2__
static inline __m256i vmin_biased(__m256i a, __m256i b)
{
    return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(a, b));
}

// n must be a multiple of 16 and >= 16
static inline uint64_t min_key(const uint64_t* keys, size_t n)
{
    const __m256i bias = _mm256_set1_epi64x(static_cast<long long>(0x8000000000000000ULL));
    auto ld = [&](size_t i) {
        return _mm256_xor_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(keys + i)), bias);
    };
    __m256i m0 = ld(0), m1 = ld(4), m2 = ld(8), m3 = ld(12);
    for (size_t i = 16; i < n; i += 16) {
        m0 = vmin_biased(m0, ld(i));
        m1 = vmin_biased(m1, ld(i + 4));
        m2 = vmin_biased(m2, ld(i + 8));
        m3 = vmin_biased(m3, ld(i + 12));
    }
    m0 = vmin_biased(vmin_biased(m0, m1), vmin_biased(m2, m3));
    alignas(32) uint64_t t[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(t), m0);
    const uint64_t r = std::min(std::min(t[0], t[1]), std::min(t[2], t[3]));
    return r ^ 0x8000000000000000ULL;
}

// n must be a multiple of 4; output indices are ascending
static inline void collect_equal(const uint64_t* keys, size_t n, uint64_t v,
                                 std::vector<uint32_t>& out)
{
    out.clear();
    const __m256i vv = _mm256_set1_epi64x(static_cast<long long>(v));
    for (size_t i = 0; i < n; i += 4) {
        const __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(keys + i));
        int mask = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(x, vv)));
        while (mask) {
            out.push_back(static_cast<uint32_t>(i + __builtin_ctz(mask)));
            mask &= mask - 1;
        }
    }
}
#else
static inline uint64_t min_key(const uint64_t* keys, size_t n)
{ return *std::min_element(keys, keys + n); }

static inline void collect_equal(const uint64_t* keys, size_t n, uint64_t v,
                                 std::vector<uint32_t>& out)
{
    out.clear();
    for (size_t i = 0; i < n; ++i)
        if (keys[i] == v) out.push_back(static_cast<uint32_t>(i));
}
#endif

static void pin_this_thread(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
        std::cerr << "Warning: could not pin merge thread to CPU " << cpu << std::endl;
}

// ============================================================
// Merge
// ============================================================
static void merge_databases(const std::vector<std::string>& databases, Sink& sink,
                            int pin_cpu, int decoder_threads)
{
    const size_t n = databases.size();
    std::cerr << "Opening " << n << " KMC databases..." << std::endl;

    DecoderPool pool(databases, decoder_threads);
    const uint32_t k = pool.k();

    // ---- output buffers ----
    const size_t max_row  = k + 2 + n * 21;
    const size_t buf_size = std::max<size_t>(1u << 22, max_row * 2);

    {
        std::string h = "kmer";
        for (auto& d : databases) { h += '\t'; h += d; }
        h += '\n';
        sink.write(h.data(), h.size());
    }

    // Decoder and writer threads start BEFORE pinning so they may use other cores.
    pool.start();
    AsyncWriter writer(sink, buf_size);
    if (pin_cpu >= 0) pin_this_thread(pin_cpu);

    AsyncWriter::Buf* cur = writer.acquire();
    char* base    = cur->data.data();
    char* p       = base;
    char* buf_end = base + buf_size;

    auto flush = [&] {
        writer.submit(cur, static_cast<size_t>(p - base));
        cur     = writer.acquire();
        base    = cur->data.data();
        p       = base;
        buf_end = base + buf_size;
    };

    std::string zeros;
    zeros.reserve(2 * n);
    for (size_t i = 0; i < n; ++i) zeros += "\t0";

    // ---- selection state: fixed slots, padded to a multiple of 16 ----
    const size_t npad = (n + 15) & ~size_t(15);
    std::vector<uint64_t> keys(npad, UINT64_MAX);
    std::vector<uint8_t>  alive(npad, 0);
    size_t remaining = 0;

    for (uint32_t i = 0; i < n; ++i)
        if (pool.advance(i)) {
            keys[i]  = pool.key(i);
            alive[i] = 1;
            ++remaining;
        }

    // Self-check: key order must agree with CKmerAPI::operator<
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < n; ++j)
            if (alive[i] && alive[j] && keys[i] < keys[j] &&
                !(pool.kmer(i) < pool.kmer(j)))
                throw std::runtime_error(
                    "kmer_data[0] is not the most significant word in this KMC version");

    const bool multiword = k > 32;
    std::vector<uint32_t> cand, matches;
    cand.reserve(n);
    matches.reserve(n);

    uint64_t out_kmers = 0;
    uint64_t total_matches = 0;      // sum over rows of "how many databases contain this k-mer"
    const auto t0 = std::chrono::steady_clock::now();

    while (remaining) {
        // ---- SIMD: smallest first word, and every slot equal to it ----
        const uint64_t mn = min_key(keys.data(), npad);
        collect_equal(keys.data(), npad, mn, cand);

        cand.erase(std::remove_if(cand.begin(), cand.end(),
                       [&](uint32_t i) { return !alive[i]; }), cand.end());

        if (!multiword) {
            matches.swap(cand);
        } else {
            uint32_t best = cand[0];
            for (size_t j = 1; j < cand.size(); ++j)
                if (pool.kmer(cand[j]) < pool.kmer(best)) best = cand[j];
            matches.clear();
            for (uint32_t i : cand)
                if (pool.kmer(i) == pool.kmer(best)) matches.push_back(i);
        }

        total_matches += matches.size();

        // ---- format row ----
        if (buf_end - p < static_cast<ptrdiff_t>(max_row)) flush();

        pool.kmer(matches[0]).to_string(p);        // writes k chars (+NUL, overwritten below)
        p += k;

        uint32_t pos = 0;
        for (uint32_t m : matches) {
            const size_t gap = m - pos;
            if (gap) { std::memcpy(p, zeros.data(), 2 * gap); p += 2 * gap; }
            *p++ = '\t';
            const uint64_t c = pool.count(m);
            if (c < 10) *p++ = static_cast<char>('0' + c);
            else        p = std::to_chars(p, p + 20, c).ptr;
            pos = m + 1;
        }
        if (pos < n) {
            const size_t gap = n - pos;
            std::memcpy(p, zeros.data(), 2 * gap);
            p += 2 * gap;
        }
        *p++ = '\n';

        // ---- advance matched streams and refresh their keys ----
        for (uint32_t m : matches) {
            if (pool.advance(m)) {
                keys[m] = pool.key(m);
            } else {
                keys[m]  = UINT64_MAX;
                alive[m] = 0;
                --remaining;
            }
        }

        if ((++out_kmers & ((1ULL << 20) - 1)) == 0) {
            const double e = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            std::cerr << "\rMerged " << out_kmers << " unique k-mers | "
                      << static_cast<uint64_t>(out_kmers / e) << " k-mers/s"
                      << " | matches/row "
                      << static_cast<double>(total_matches) / static_cast<double>(out_kmers)
                      << " | stalls " << pool.stalls()
                      << "   " << std::flush;
        }
    }

    writer.submit(cur, static_cast<size_t>(p - base));
    writer.finish();     // drains the queue, rethrows write errors
    sink.close();        // waits for the compressor

    const double e = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cerr << "\nFinished.\nUnique k-mers written: " << out_kmers
              << "\nElapsed time: " << e << " seconds"
              << "\nAverage rate: " << (e > 0 ? out_kmers / e : 0.0)
              << " unique k-mers/s"
              << "\nAverage matches per row: "
              << (out_kmers ? static_cast<double>(total_matches) / out_kmers : 0.0)
              << "\nDecoder stalls (merge thread waited for a batch): " << pool.stalls()
              << std::endl;
}

// ============================================================
// Main
// ============================================================
static void usage(const char* prog)
{
    std::cerr << "Usage:\n  " << prog
              << " [-o FILE] [-z] [-p N] [-l N] [-B KiB] [-c \"CMD\"] [-t N] [-m CPU] db1 db2 ...\n"
              << "  -o FILE   output file (default: stdout)\n"
              << "  -z        compress with pigz (implied by .gz extension)\n"
              << "  -p N      pigz threads (default: cores - 1)\n"
              << "  -l N      pigz level 1-9 (default: 4)\n"
              << "  -B KiB    pigz block size in KiB (default: 1024)\n"
              << "  -c CMD    custom compressor command writing to stdout\n"
              << "  -t N      decoder threads (default: 3)\n"
              << "  -m CPU    pin merge thread to logical CPU\n";
}

int main(int argc, char** argv)
{
    signal(SIGPIPE, SIG_IGN);

    try {
        std::string out_path;
        bool compress = false;
        int ncpu = static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN));
        int threads = std::max(1, ncpu - 1);
        int level = 4;
        int block_kib = 1024;
        int pin_cpu = -1;
        int decoders = 3;
        std::vector<std::string> cmd;

        std::vector<std::string> dbs;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if      (a == "-o" && i + 1 < argc) out_path = argv[++i];
            else if (a == "-z")                 compress = true;
            else if (a == "-p" && i + 1 < argc) threads = std::max(1, std::atoi(argv[++i]));
            else if (a == "-l" && i + 1 < argc) level = std::clamp(std::atoi(argv[++i]), 1, 9);
            else if (a == "-B" && i + 1 < argc) block_kib = std::max(32, std::atoi(argv[++i]));
            else if (a == "-t" && i + 1 < argc) decoders = std::max(1, std::atoi(argv[++i]));
            else if (a == "-m" && i + 1 < argc) pin_cpu = std::atoi(argv[++i]);
            else if (a == "-c" && i + 1 < argc) {
                std::istringstream ss(argv[++i]);
                for (std::string w; ss >> w; ) cmd.push_back(w);
                compress = true;
            }
            else if (a == "-h" || a == "--help") { usage(argv[0]); return EXIT_SUCCESS; }
            else dbs.push_back(a);
        }

        if (dbs.empty()) { usage(argv[0]); return EXIT_FAILURE; }

        if (out_path.size() >= 3 &&
            out_path.compare(out_path.size() - 3, 3, ".gz") == 0)
            compress = true;

        if (compress && cmd.empty())
            cmd = {"pigz", "-c", "-p", std::to_string(threads), "-" + std::to_string(level),
                   "-b", std::to_string(block_kib), "-i"};

        std::cerr << "KMC streaming merge\n===================\n"
                  << "Databases: " << dbs.size() << "\n"
                  << "Decoder threads: " << decoders << "\n"
                  << "Output:    " << (out_path.empty() ? "stdout" : out_path) << "\n";
        if (compress) {
            std::cerr << "Compressor:";
            for (const auto& w : cmd) std::cerr << ' ' << w;
            std::cerr << "\n";
        }
        std::cerr << std::flush;

        Sink sink;
        sink.open(out_path, compress, cmd);
        merge_databases(dbs, sink, pin_cpu, decoders);
        return EXIT_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
