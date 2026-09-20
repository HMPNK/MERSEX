// kmc_merge: streaming merge of many KMC2/KMC3 databases into one count table.
//
// Replaces N x "kmc_tools transform <db> dump -s" + paste-style merging.
// Each database is read in globally sorted k-mer order by an internal per-database
// bin heap (the same algorithm as kmc_tools' CKMC2DbReaderSorted), prefetched by a
// small shared thread pool into bounded batches.
//
// Build:
//   g++ -std=c++17 -O3 -march=native -flto -pthread -Ikmc_api kmc_merge.cpp kmc_api/*.cpp -o kmc_merge
//   (add -DVERIFY_SCAN for a slow scalar cross-check of the SIMD minimum)
//
// Usage:
//   kmc_merge [options] db1 db2 ...
//
//   -o FILE   write to FILE instead of stdout
//   -z        compress with pigz (implied if FILE ends in .gz)
//   -p N      pigz threads (default: cores - 1)
//   -l N      pigz level 1-9 (default: 4)
//   -B KiB    pigz block size in KiB (default: 1024)
//   -c "CMD"  custom compressor command instead of pigz; must write to stdout
//   -m CPU    pin the merge thread to logical CPU number CPU
//   -t N      database reader threads (default: min(4, cores-2))
//   -b N      k-mers per prefetch batch (default: 16384)
//   -S BYTES  suffix read buffer per bin (default: 2048)
//   -L N      LUT entries buffered per bin (default: 256)
//
// Output: header line "kmer<TAB>db1<TAB>db2...", then one row per distinct k-mer in
// lexicographic order (A<C<G<T), counts 0 where a database lacks the k-mer.
//
// Only KMC2/KMC3 databases (.kmc_pre/.kmc_suf with KMCP/KMCS markers) are supported.

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
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
#include <functional>
#include <iostream>
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
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "kmc_file.h"   // public API only; used just to read the database header (Info)

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
// Small fixed thread pool
// ============================================================
class TaskPool {
public:
    explicit TaskPool(int n)
    {
        for (int i = 0; i < std::max(1, n); ++i)
            thr_.emplace_back([this] { loop(); });
    }
    ~TaskPool()
    {
        { std::lock_guard<std::mutex> lk(m_); stop_ = true; }
        cv_.notify_all();
        for (auto& t : thr_) t.join();
    }
    void post(std::function<void()> f)
    {
        { std::lock_guard<std::mutex> lk(m_); q_.push_back(std::move(f)); }
        cv_.notify_one();
    }
private:
    void loop()
    {
        for (;;) {
            std::function<void()> f;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [&] { return stop_ || !q_.empty(); });
                if (q_.empty()) return;          // drain, then exit
                f = std::move(q_.front());
                q_.pop_front();
            }
            f();
        }
    }
    std::vector<std::thread> thr_;
    std::deque<std::function<void()>> q_;
    std::mutex m_;
    std::condition_variable cv_;
    bool stop_ = false;
};

// ============================================================
// Sorted stream of one KMC2/KMC3 database (= kmc_tools "sorted" mode)
//
// File layout used (from kmc_tools' CKMC2DbReaderSorted):
//   .kmc_pre : "KMCP" | uint64 LUT[no_of_bins * lut_size + 1] (global cumulative
//              k-mer numbers) | ...
//   .kmc_suf : "KMCS" | records of (suffix_bytes + counter_size) bytes
//   Bin i owns records [LUT[i*lut_size], LUT[(i+1)*lut_size]); inside a bin, prefix p
//   owns [LUT[i*lut_size+p], LUT[i*lut_size+p+1]). Each bin is sorted, bins are not
//   ordered relative to each other, so a heap over the bins gives the global order.
//
// K-mers are produced as big-endian bytes [prefix bytes][suffix bytes], so memcmp
// order equals lexicographic A<C<G<T order.
// ============================================================
static void pread_full(int fd, void* buf, size_t n, uint64_t off)
{
    char* p = static_cast<char*>(buf);
    while (n) {
        const ssize_t r = ::pread(fd, p, n, static_cast<off_t>(off));
        if (r < 0 && errno == EINTR) continue;
        if (r <= 0) throw std::runtime_error("pread failed / unexpected EOF");
        p += r; n -= static_cast<size_t>(r); off += static_cast<uint64_t>(r);
    }
}

class SortedDb {
public:
    uint32_t k = 0, kmer_bytes = 0;

    SortedDb(const std::string& path, size_t suf_buf_bytes, size_t lut_buf_entries)
        : path_(path), lut_cfg_(std::max<size_t>(2, lut_buf_entries))
    {
        uint32_t kl = 0, mode = 0, cs = 0, lut = 0, sig = 0, minc = 0;
        uint64_t maxc = 0, total = 0;
        {   // header only; CKMCFile is closed again immediately
            CKMCFile db;
            if (!db.OpenForListing(path)) throw std::runtime_error("Cannot open KMC database: " + path);
            db.Info(kl, mode, cs, lut, sig, minc, maxc, total);
            db.Close();
        }
        k = kl;
        counter_size_ = cs;
        total_ = total;
        if (lut >= k || (k - lut) % 4) throw std::runtime_error("unsupported KMC layout: " + path);
        lut_size_     = 1ull << (2 * lut);
        prefix_bytes_ = (lut + 3) / 4;
        suffix_bytes_ = (k - lut) / 4;
        kmer_bytes    = prefix_bytes_ + suffix_bytes_;
        rec_size_     = suffix_bytes_ + counter_size_;
        cfg_recs_     = std::max<size_t>(1, suf_buf_bytes / rec_size_);

        fd_pre_ = ::open((path + ".kmc_pre").c_str(), O_RDONLY);
        fd_suf_ = ::open((path + ".kmc_suf").c_str(), O_RDONLY);
        if (fd_pre_ < 0 || fd_suf_ < 0)
            throw std::runtime_error("Cannot open .kmc_pre/.kmc_suf: " + path);
        char m[4];
        pread_full(fd_pre_, m, 4, 0);
        if (std::memcmp(m, "KMCP", 4)) throw std::runtime_error("bad .kmc_pre marker (KMC1 db?): " + path);
        pread_full(fd_suf_, m, 4, 0);
        if (std::memcmp(m, "KMCS", 4)) throw std::runtime_error("bad .kmc_suf marker: " + path);
    }
    ~SortedDb()
    {
        if (fd_pre_ >= 0) ::close(fd_pre_);
        if (fd_suf_ >= 0) ::close(fd_suf_);
    }
    SortedDb(const SortedDb&) = delete;
    SortedDb& operator=(const SortedDb&) = delete;

    // Returns kmer_bytes big-endian bytes (valid until the next call), or nullptr at end.
    const uint8_t* next(uint64_t& count)
    {
        init();
        if (pending_ >= 0) {                       // lazy advance keeps the last pointer valid
            if (advance(static_cast<uint32_t>(pending_))) {
                sift_down(0);
            } else {
                heap_[0] = heap_.back();
                heap_.pop_back();
                if (!heap_.empty()) sift_down(0);
            }
            pending_ = -1;
        }
        if (heap_.empty()) return nullptr;
        const uint32_t bi = heap_[0];
        pending_ = bi;
        count = bins_[bi].count;
        return &curs_[size_t(bi) * kmer_bytes];
    }

private:
    struct Bin {
        uint64_t rec_next = 0, rec_end = 0;     // records still to buffer
        uint64_t left = 0, left_in_prefix = 0;
        uint64_t prefix = ~0ull;                // ++ wraps to 0
        uint64_t lut0 = 0, lut_lo = 0;          // lut0 = bin * lut_size_
        uint64_t count = 0;
        uint32_t lut_n = 0, lut_cap = 0, buf_n = 0, buf_pos = 0, buf_cap = 0;
        uint8_t*  sbuf = nullptr;
        uint64_t* lbuf = nullptr;
    };

    uint64_t lut_entry_file(uint64_t idx)
    {
        uint64_t v;
        pread_full(fd_pre_, &v, 8, 4 + 8 * idx);
        return v;
    }

    void init()
    {
        if (inited_) return;
        inited_ = true;
        if (!total_) return;

        struct stat st;
        if (fstat(fd_pre_, &st) != 0) throw std::runtime_error("fstat failed");

        // Number of bins: first b whose boundary entry already equals total_kmers.
        // Trailing empty bins are harmless to drop.
        uint64_t nb = 0;
        for (uint64_t b = 1;; ++b) {
            const uint64_t idx = b * lut_size_;
            if (4 + 8 * (idx + 1) > static_cast<uint64_t>(st.st_size))
                throw std::runtime_error("cannot determine bin count: " + path_);
            const uint64_t v = lut_entry_file(idx);
            if (v == total_) { nb = b; break; }
            if (v > total_)  throw std::runtime_error("corrupt LUT: " + path_);
        }

        std::vector<uint64_t> bound(nb + 1);
        for (uint64_t b = 0; b <= nb; ++b) bound[b] = lut_entry_file(b * lut_size_);

        bins_.resize(nb);
        curs_.assign(size_t(nb) * kmer_bytes, 0);
        size_t s_tot = 0, l_tot = 0;
        for (uint64_t b = 0; b < nb; ++b) {
            Bin& x = bins_[b];
            const uint64_t cnt = bound[b + 1] - bound[b];
            x.rec_next = bound[b];
            x.rec_end  = bound[b + 1];
            x.left     = cnt;
            x.lut0     = b * lut_size_;
            if (!cnt) continue;
            x.buf_cap = static_cast<uint32_t>(std::min<uint64_t>(cfg_recs_, cnt));
            x.lut_cap = static_cast<uint32_t>(std::min<uint64_t>(lut_cfg_, lut_size_ + 1));
            s_tot += size_t(x.buf_cap) * rec_size_;
            l_tot += x.lut_cap;
        }
        sarena_.resize(s_tot + 8);
        larena_.resize(l_tot);
        size_t so = 0, lo = 0;
        for (auto& x : bins_) {
            if (!x.left) continue;
            x.sbuf = sarena_.data() + so; so += size_t(x.buf_cap) * rec_size_;
            x.lbuf = larena_.data() + lo; lo += x.lut_cap;
        }

        heap_.reserve(nb);
        for (uint32_t b = 0; b < nb; ++b)
            if (bins_[b].left && advance(b)) heap_.push_back(b);
        std::make_heap(heap_.begin(), heap_.end(),
                       [&](uint32_t a, uint32_t b) { return gt(a, b); });
    }

    uint64_t lut_at(Bin& b, uint64_t idx)       // idx in [0, lut_size_]
    {
        if (idx < b.lut_lo || idx >= b.lut_lo + b.lut_n) {
            const uint64_t n = std::min<uint64_t>(b.lut_cap, lut_size_ + 1 - idx);
            pread_full(fd_pre_, b.lbuf, n * 8, 4 + 8 * (b.lut0 + idx));
            b.lut_lo = idx;
            b.lut_n  = static_cast<uint32_t>(n);
        }
        return b.lbuf[idx - b.lut_lo];
    }

    void refill(Bin& b)
    {
        const uint64_t n = std::min<uint64_t>(b.buf_cap, b.rec_end - b.rec_next);
        pread_full(fd_suf_, b.sbuf, n * rec_size_, 4 + b.rec_next * rec_size_);
        b.rec_next += n;
        b.buf_n   = static_cast<uint32_t>(n);
        b.buf_pos = 0;
    }

    bool advance(uint32_t bi)
    {
        Bin& b = bins_[bi];
        if (!b.left) return false;
        uint8_t* cur = &curs_[size_t(bi) * kmer_bytes];
        while (!b.left_in_prefix) {
            if (++b.prefix >= lut_size_) throw std::runtime_error("corrupt LUT: " + path_);
            const uint64_t a = lut_at(b, b.prefix);
            b.left_in_prefix = lut_at(b, b.prefix + 1) - a;
            if (b.left_in_prefix)
                for (uint32_t i = 0; i < prefix_bytes_; ++i)
                    cur[i] = static_cast<uint8_t>(b.prefix >> (8 * (prefix_bytes_ - 1 - i)));
        }
        if (b.buf_pos == b.buf_n) refill(b);
        const uint8_t* r = b.sbuf + size_t(b.buf_pos) * rec_size_;
        std::memcpy(cur + prefix_bytes_, r, suffix_bytes_);       // suffix stored MSB first
        uint64_t c = 0;
        for (uint32_t i = counter_size_; i-- > 0;)
            c = (c << 8) | r[suffix_bytes_ + i];                   // little-endian counter
        b.count = counter_size_ ? c : 1;
        ++b.buf_pos; --b.left; --b.left_in_prefix;
        return true;
    }

    bool gt(uint32_t a, uint32_t b) const
    {
        return std::memcmp(&curs_[size_t(a) * kmer_bytes],
                           &curs_[size_t(b) * kmer_bytes], kmer_bytes) > 0;
    }

    void sift_down(size_t i)                    // min-heap, same layout as std::make_heap
    {
        const size_t n = heap_.size();
        const uint32_t x = heap_[i];
        for (;;) {
            size_t c = 2 * i + 1;
            if (c >= n) break;
            if (c + 1 < n && gt(heap_[c], heap_[c + 1])) ++c;
            if (!gt(x, heap_[c])) break;
            heap_[i] = heap_[c];
            i = c;
        }
        heap_[i] = x;
    }

    std::string path_;
    size_t   lut_cfg_, cfg_recs_ = 1;
    uint32_t counter_size_ = 0, prefix_bytes_ = 0, suffix_bytes_ = 0, rec_size_ = 0;
    uint64_t total_ = 0, lut_size_ = 0;
    int fd_pre_ = -1, fd_suf_ = -1;
    bool inited_ = false;
    int64_t pending_ = -1;
    std::vector<Bin>      bins_;
    std::vector<uint8_t>  curs_, sarena_;
    std::vector<uint64_t> larena_;
    std::vector<uint32_t> heap_;
};

// ============================================================
// Reader with background prefetch (bounded memory)
//   - a pool worker pulls sorted k-mers from SortedDb into batches
//   - the merge thread only reads from ready batches
//   - at most one worker touches a given database at a time
// K-mers are exposed as W_ uint64 words, word 0 most significant, right-aligned.
// ============================================================
struct Reader {
    uint64_t count = 0;
    uint32_t k = 0;
    const uint64_t* words() const { return cw_; }

    Reader(const std::string& fn, TaskPool& pool, size_t batch_kmers,
           size_t suf_bytes, size_t lut_entries)
        : sd_(fn, suf_bytes, lut_entries), pool_(pool), cap_(std::max<size_t>(batch_kmers, 1024))
    {
        k   = sd_.k;
        kb_ = sd_.kmer_bytes;
        W_  = (k + 31) / 32;
        prev_.resize(kb_);
        store_.resize(kBatches);
        for (auto& b : store_) {
            b.words.resize(cap_ * W_);
            b.counts.resize(cap_);
            free_.push_back(&b);
        }
        std::lock_guard<std::mutex> lk(m_);
        schedule_locked();                       // start prefetching immediately
    }

    ~Reader()
    {
        std::unique_lock<std::mutex> lk(m_);
        closing_ = true;
        cv_.wait(lk, [&] { return !scheduled_; });   // no worker touches us after this
    }
    Reader(const Reader&) = delete;
    Reader& operator=(const Reader&) = delete;

    bool less(const Reader& o) const
    {
        for (size_t i = 0; i < W_; ++i)
            if (cw_[i] != o.cw_[i]) return cw_[i] < o.cw_[i];
        return false;
    }
    bool equal(const Reader& o) const { return std::memcmp(cw_, o.cw_, W_ * 8) == 0; }

    bool next()
    {
        if (cur_ && pos_ < cur_->n) { load(); return true; }   // hot path, no locks
        return next_slow();
    }

private:
    static constexpr size_t kBatches = 3;      // 1 being consumed + 2 in flight

    struct Batch {
        std::vector<uint64_t> words;
        std::vector<uint64_t> counts;
        size_t n = 0;
    };

    inline void load()
    {
        cw_   = &cur_->words[pos_ * W_];
        count = cur_->counts[pos_];
        ++pos_;
    }

    bool next_slow()
    {
        std::unique_lock<std::mutex> lk(m_);
        if (cur_) {
            cur_->n = 0;
            free_.push_back(cur_);
            cur_ = nullptr;
            schedule_locked();
        }
        cv_.wait(lk, [&] { return !ready_.empty() || eof_ || err_; });
        if (ready_.empty()) {
            if (err_) std::rethrow_exception(err_);
            return false;                                  // clean EOF
        }
        cur_ = ready_.front();
        ready_.pop_front();
        pos_ = 0;
        lk.unlock();
        load();
        return true;
    }

    void schedule_locked()
    {
        if (!scheduled_ && !eof_ && !closing_ && !err_ && !free_.empty()) {
            scheduled_ = true;
            pool_.post([this] { fill(); });
        }
    }

    // big-endian bytes -> right-aligned words (word 0 most significant)
    void pack(const uint8_t* kb, uint64_t* dst) const
    {
        std::memset(dst, 0, W_ * 8);
        const size_t off = 8 * W_ - kb_;
        for (size_t q = 0; q < kb_; ++q) {
            const size_t p = off + q;
            dst[p >> 3] |= uint64_t(kb[q]) << (8 * (7 - (p & 7)));
        }
    }

    // runs on a pool thread; only one instance per Reader at any time
    void fill()
    {
        try {
            for (;;) {
                Batch* b;
                {
                    std::lock_guard<std::mutex> lk(m_);
                    if (closing_ || eof_ || free_.empty()) {
                        scheduled_ = false;
                        cv_.notify_all();          // notify under lock: ~Reader may follow
                        return;
                    }
                    b = free_.back();
                    free_.pop_back();
                }

                size_t n = 0;
                bool eof = false;
                uint64_t c = 0;
                while (n < cap_) {
                    const uint8_t* kb = sd_.next(c);
                    if (!kb) { eof = true; break; }
                    if (have_prev_ && std::memcmp(prev_.data(), kb, kb_) >= 0)
                        throw std::runtime_error(
                            "internal error: k-mers not strictly increasing in a database "
                            "(bin order assumption violated)");
                    std::memcpy(prev_.data(), kb, kb_);
                    have_prev_ = true;
                    pack(kb, &b->words[n * W_]);
                    b->counts[n] = c;
                    ++n;
                }

                {
                    std::lock_guard<std::mutex> lk(m_);
                    b->n = n;
                    if (n) ready_.push_back(b); else free_.push_back(b);
                    if (eof) eof_ = true;
                }
                cv_.notify_all();
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(m_);
            err_ = std::current_exception();
            scheduled_ = false;
            cv_.notify_all();
        }
    }

    SortedDb    sd_;                   // producer side
    TaskPool&   pool_;
    size_t      cap_, W_ = 1, kb_ = 0;
    const uint64_t* cw_ = nullptr;     // consumer side: current k-mer words
    std::vector<uint8_t> prev_;
    bool        have_prev_ = false;

    std::vector<Batch>  store_;
    std::deque<Batch*>  ready_;
    std::vector<Batch*> free_;
    Batch*  cur_ = nullptr;            // consumer-only
    size_t  pos_ = 0;                  // consumer-only

    std::mutex m_;
    std::condition_variable cv_;
    bool eof_ = false, scheduled_ = false, closing_ = false;
    std::exception_ptr err_;
};

// ============================================================
// K-mer text formatting (same table method as kmc_tools' CDumpWriterBase)
// ============================================================
struct KmerFmt {
    uint32_t W = 1, nbytes = 0, first = 4;
    char tab[1024];

    void init(uint32_t k)
    {
        W = (k + 31) / 32;
        nbytes = (k + 3) / 4;
        first = k % 4 ? k % 4 : 4;
        const char c[4] = {'A', 'C', 'G', 'T'};
        for (uint32_t v = 0; v < 256; ++v)
            for (int j = 0; j < 4; ++j) tab[v * 4 + j] = c[(v >> (6 - 2 * j)) & 3];
    }

    char* write(char* p, const uint64_t* w) const
    {
        auto byte = [&](uint32_t q) {
            const uint32_t bp = 8 * (nbytes - 1 - q);
            return static_cast<uint32_t>((w[W - 1 - (bp >> 6)] >> (bp & 63)) & 0xFF);
        };
        const char* t = tab + 4 * byte(0) + 4 - first;
        for (uint32_t i = 0; i < first; ++i) *p++ = *t++;
        for (uint32_t q = 1; q < nbytes; ++q) { std::memcpy(p, tab + 4 * byte(q), 4); p += 4; }
        return p;
    }
};

// ============================================================
// SIMD helpers (AVX2, with scalar fallback)
// ============================================================
#ifdef __AVX2__
// inputs must already be biased by the sign bit (unsigned compare via signed compare)
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
    __m256i m0 = ld(0), m1 = ld(4), m2 = ld(8), m3 = ld(12);   // 4 independent chains
    for (size_t i = 16; i < n; i += 16) {
        m0 = vmin_biased(m0, ld(i));
        m1 = vmin_biased(m1, ld(i + 4));
        m2 = vmin_biased(m2, ld(i + 8));
        m3 = vmin_biased(m3, ld(i + 12));
    }
    m0 = vmin_biased(vmin_biased(m0, m1), vmin_biased(m2, m3));
    alignas(32) uint64_t t[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(t), m0);
    // Undo the bias BEFORE the scalar (unsigned) reduction.
    const uint64_t B = 0x8000000000000000ULL;
    return std::min(std::min(t[0] ^ B, t[1] ^ B), std::min(t[2] ^ B, t[3] ^ B));
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
                            int pin_cpu, int reader_threads, size_t batch_kmers,
                            size_t suf_bytes, size_t lut_entries)
{
    const size_t n = databases.size();
    std::cerr << "Opening " << n << " KMC databases (" << reader_threads
              << " reader threads)..." << std::endl;

    TaskPool pool(reader_threads);   // declared BEFORE readers => destroyed after them,
                                     // and created before pinning so workers stay unpinned
    std::vector<std::unique_ptr<Reader>> readers;
    readers.reserve(n);
    for (const auto& d : databases) {
        std::cerr << "  " << d << std::endl;
        readers.emplace_back(std::make_unique<Reader>(d, pool, batch_kmers, suf_bytes, lut_entries));
    }

    const uint32_t k = readers[0]->k;
    for (auto& r : readers)
        if (r->k != k) throw std::runtime_error("Databases have different k");

    KmerFmt fmt;
    fmt.init(k);

    // ---- output buffers ----
    const size_t max_row  = k + 2 + n * 21;
    const size_t buf_size = std::max<size_t>(1u << 22, max_row * 2);

    {
        std::string h = "kmer";
        for (auto& d : databases) { h += '\t'; h += d; }
        h += '\n';
        sink.write(h.data(), h.size());
    }

    AsyncWriter writer(sink, buf_size);       // writer thread starts BEFORE pinning,
                                              // so it is free to run on other cores
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

    // "\t0\t0\t0..." so zero runs are one memcpy
    std::string zeros;
    zeros.reserve(2 * n);
    for (size_t i = 0; i < n; ++i) zeros += "\t0";

    // ---- selection state: fixed slots, padded to a multiple of 16 ----
    const size_t npad = (n + 15) & ~size_t(15);
    std::vector<uint64_t> keys(npad, UINT64_MAX);
    std::vector<uint8_t>  alive(npad, 0);
    size_t remaining = 0;

    for (uint32_t i = 0; i < n; ++i)
        if (readers[i]->next()) {
            keys[i]  = readers[i]->words()[0];
            alive[i] = 1;
            ++remaining;
        }

    const bool multiword = k > 32;
    std::vector<uint32_t> cand, matches;
    uint64_t prev_mn = 0;
    cand.reserve(n);
    matches.reserve(n);

    uint64_t out_kmers = 0;
    const auto t0 = std::chrono::steady_clock::now();

    while (remaining) {
        // ---- SIMD: smallest first word, and every slot equal to it ----
        const uint64_t mn = min_key(keys.data(), npad);

        if (mn < prev_mn)
            throw std::runtime_error("internal error: k-mers out of order (scan bug)");
        prev_mn = mn;
#ifdef VERIFY_SCAN
        {
            uint64_t ref = UINT64_MAX;
            for (size_t i = 0; i < n; ++i) if (alive[i]) ref = std::min(ref, keys[i]);
            if (ref != mn) throw std::runtime_error("VERIFY_SCAN: SIMD minimum mismatch");
        }
#endif
        collect_equal(keys.data(), npad, mn, cand);

        // dead/padding slots hold UINT64_MAX and can tie only when mn == UINT64_MAX
        cand.erase(std::remove_if(cand.begin(), cand.end(),
                       [&](uint32_t i) { return !alive[i]; }), cand.end());

        if (cand.empty())
            throw std::runtime_error("internal error: no candidate for the minimum key");

        if (!multiword) {
            matches.swap(cand);
        } else {
            // equal first words don't imply equal k-mers: resolve with full compare
            uint32_t best = cand[0];
            for (size_t j = 1; j < cand.size(); ++j)
                if (readers[cand[j]]->less(*readers[best])) best = cand[j];
            matches.clear();
            for (uint32_t i : cand)
                if (readers[i]->equal(*readers[best])) matches.push_back(i);
        }

        // ---- format row ----
        if (buf_end - p < static_cast<ptrdiff_t>(max_row)) flush();

        p = fmt.write(p, readers[matches[0]]->words());

        uint32_t pos = 0;
        for (uint32_t m : matches) {
            const size_t gap = m - pos;
            if (gap) { std::memcpy(p, zeros.data(), 2 * gap); p += 2 * gap; }
            *p++ = '\t';
            const uint64_t c = readers[m]->count;
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

        // ---- advance matched readers and refresh their keys ----
        for (uint32_t m : matches) {
            if (readers[m]->next()) {
                keys[m] = readers[m]->words()[0];
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
                      << std::flush;
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
              << " unique k-mers/s" << std::endl;
}

// ============================================================
// Main
// ============================================================
static void usage(const char* prog)
{
    std::cerr << "Usage:\n  " << prog
              << " [-o FILE] [-z] [-p N] [-l N] [-B KiB] [-c \"CMD\"] [-m CPU]"
                 " [-t N] [-b N] [-S BYTES] [-L N] db1 db2 ...\n"
              << "  -o FILE   output file (default: stdout)\n"
              << "  -z        compress with pigz (implied by .gz extension)\n"
              << "  -p N      pigz threads (default: cores - 1)\n"
              << "  -l N      pigz level 1-9 (default: 4)\n"
              << "  -B KiB    pigz block size in KiB (default: 1024)\n"
              << "  -c CMD    custom compressor command writing to stdout (overrides pigz options)\n"
              << "  -m CPU    pin merge thread to logical CPU\n"
              << "  -t N      database reader threads (default: min(4, cores-2))\n"
              << "  -b N      k-mers per prefetch batch (default: 16384)\n"
              << "  -S BYTES  suffix read buffer per bin (default: 2048)\n"
              << "  -L N      LUT entries buffered per bin (default: 256)\n";
}

int main(int argc, char** argv)
{
    signal(SIGPIPE, SIG_IGN);   // report write errors instead of dying silently

    try {
        std::string out_path;
        bool compress = false;
        int ncpu = static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN));
        int threads = std::max(1, ncpu - 1);
        int level = 4;
        int block_kib = 1024;
        int pin_cpu = -1;
        int reader_threads = std::max(1, std::min(4, ncpu - 2));
        int batch_kmers = 16384;
        int suf_bytes = 2048;
        int lut_entries = 256;
        std::vector<std::string> cmd;

        std::vector<std::string> dbs;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if      (a == "-o" && i + 1 < argc) out_path = argv[++i];
            else if (a == "-z")                 compress = true;
            else if (a == "-p" && i + 1 < argc) threads = std::max(1, std::atoi(argv[++i]));
            else if (a == "-l" && i + 1 < argc) level = std::clamp(std::atoi(argv[++i]), 1, 9);
            else if (a == "-B" && i + 1 < argc) block_kib = std::max(32, std::atoi(argv[++i]));
            else if (a == "-m" && i + 1 < argc) pin_cpu = std::atoi(argv[++i]);
            else if (a == "-t" && i + 1 < argc) reader_threads = std::max(1, std::atoi(argv[++i]));
            else if (a == "-b" && i + 1 < argc) batch_kmers = std::max(1024, std::atoi(argv[++i]));
            else if (a == "-S" && i + 1 < argc) suf_bytes = std::max(64, std::atoi(argv[++i]));
            else if (a == "-L" && i + 1 < argc) lut_entries = std::max(2, std::atoi(argv[++i]));
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
                  << "Output:    " << (out_path.empty() ? "stdout" : out_path) << "\n";
        if (compress) {
            std::cerr << "Compressor:";
            for (const auto& w : cmd) std::cerr << ' ' << w;
            std::cerr << "\n";
        }
        std::cerr << std::flush;

        Sink sink;
        sink.open(out_path, compress, cmd);
        merge_databases(dbs, sink, pin_cpu, reader_threads,
                        static_cast<size_t>(batch_kmers),
                        static_cast<size_t>(suf_bytes),
                        static_cast<size_t>(lut_entries));
        return EXIT_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
