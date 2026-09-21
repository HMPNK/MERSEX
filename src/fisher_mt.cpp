// fisher_mt.cpp - multithreaded port of the AWK fisher_log script
//
// Build:  g++ -O2 -std=c++17 -pthread -o fisher_mt fisher_mt.cpp
// Usage:  fisher_mt --fcount N [--cmin1 X] [--cmin2 Y] [--pcut P] [--alt two|greater|less|point]
//                    [--header] [--threads T] [--batch B] < in.tsv > out.tsv
//
// Input : tab-separated table, column 1 = id, columns 2..fcount+1 = group 1,
//         remaining columns = group 2.
// With --header the first input line is a header: it is not analysed, and is
//         written to the output with the columns pval, grp1_high, grp1_low,
//         grp2_high, grp2_low prepended.
// Output: Fisher exact p \t a \t b \t c \t d \t original_line   (input order preserved)

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <deque>
#include <fcntl.h>
#include <getopt.h>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

struct Params {
    double cmin1 = 1;
    double cmin2 = 0;
    int    alt = 0;            // 0 = two-sided, 1 = less, 2 = greater, 3 = point probability
    long   fcount = -1;
    bool   header = false;
    bool   stats = false;
    bool   use_pcut = false;
    double pcut = 1.0;
    unsigned threads = 0;      // 0 = 80% of hardware concurrency
    size_t batch = 1000000;    // lines per work unit
};

// Raw malloc'd byte buffer: unlike std::vector<char> it is never zero-filled and
// grows with realloc (no copy for large blocks), which matters for huge batches.
struct RawBuf {
    char *p = nullptr;
    size_t n = 0, cap = 0;
    RawBuf() = default;
    RawBuf(const RawBuf &) = delete;
    RawBuf &operator=(const RawBuf &) = delete;
    RawBuf(RawBuf &&o) noexcept : p(o.p), n(o.n), cap(o.cap) { o.p = nullptr; o.n = o.cap = 0; }
    RawBuf &operator=(RawBuf &&o) noexcept {
        if (this != &o) { std::free(p); p = o.p; n = o.n; cap = o.cap; o.p = nullptr; o.n = o.cap = 0; }
        return *this;
    }
    ~RawBuf() { std::free(p); }
    void reserve(size_t c) {
        if (c <= cap) return;
        size_t nc = std::max(c, cap + cap / 2);
        char *q = (char *)std::realloc(p, nc);
        if (!q) { std::fprintf(stderr, "out of memory\n"); std::abort(); }
        p = q; cap = nc;
    }
    void release() { std::free(p); p = nullptr; n = cap = 0; }
};

// A batch = a block of complete, '\n'-terminated lines kept as one raw buffer
// (no per-line allocation), plus the formatted output for those lines.
struct Batch {
    size_t seq = 0;
    RawBuf data;
    std::string out;
};

// Per-thread cumulative log-factorial table (same summation order as the AWK code).
struct LogFact {
    std::vector<double> t{0.0, 0.0};
    double get(size_t n) {
        while (t.size() <= n) t.push_back(t.back() + std::log((double)t.size()));
        return t[n];
    }
};

// True Fisher exact test on the 2x2 table [[a,b],[c,d]].
//   two-sided: sum of P(x) over all tables (same margins) with P(x) <= P(observed)
//              (relative tolerance 1e-7, same convention as R's fisher.test)
//   greater  : one-sided, P(X >= a)  (= hypergeometric test, enrichment in group 1)
//   less     : one-sided, P(X <= a)  (depletion in group 1)
//   point    : point probability P(X = a) only (what the original AWK script computed;
//              NOT a p-value)
static double fisher_p(LogFact &lf, long a, long b, long c, long d, int alt) {
    long r1 = a + b, r2 = c + d, c1 = a + c, c2 = b + d, n = r1 + r2;
    long lo = std::max(0L, c1 - r2), hi = std::min(r1, c1);
    double konst = lf.get(r1) + lf.get(r2) + lf.get(c1) + lf.get(c2) - lf.get(n);
    auto logp = [&](long x) {
        return konst - (lf.get(x) + lf.get(r1 - x) + lf.get(c1 - x) + lf.get(r2 - c1 + x));
    };
    double sum = 0.0;
    if (alt == 3) return std::exp(logp(a));
    if (alt == 1) {
        for (long x = lo; x <= a; x++) sum += std::exp(logp(x));
    } else if (alt == 2) {
        for (long x = a; x <= hi; x++) sum += std::exp(logp(x));
    } else {
        double thr = logp(a) + 1e-7;
        for (long x = lo; x <= hi; x++) {
            double lp = logp(x);
            if (lp <= thr) sum += std::exp(lp);
        }
    }
    return sum > 1.0 ? 1.0 : sum;
}

// Per-thread cache: many rows share the same 2x2 table.
struct PCache {
    std::unordered_map<uint64_t, double> m;
};

// Fast field -> number. Plain integers take a fast path; anything else
// (decimals, exponents, nan/inf, empty, text) falls back to strtod on a copy.
// Non-numeric text becomes 0, as before.
static inline double parse_num(const char *s, const char *e) {
    const char *q = s;
    while (q < e && *q == ' ') q++;
    bool neg = false;
    if (q < e && (*q == '-' || *q == '+')) { neg = (*q == '-'); q++; }
    const char *ds = q;
    long long iv = 0;
    while (q < e && (unsigned)(*q - '0') < 10u) { iv = iv * 10 + (*q - '0'); q++; }
    size_t nd = (size_t)(q - ds);
    if (nd > 0 && nd <= 15 && (q == e || (*q != '.' && *q != 'e' && *q != 'E')))
        return neg ? -(double)iv : (double)iv;
    std::string tmp(s, e);
    return std::strtod(tmp.c_str(), nullptr);
}

static inline void process_line(const char *b, const char *e, const Params &P, LogFact &lf,
                                PCache &pc, std::string &out) {
    long f = 0, f2 = 0, m = 0, m2 = 0;
    long field = 1;                       // 1-based field index, like AWK
    const char *s = b;
    while (true) {
        const char *t = (const char *)std::memchr(s, '\t', (size_t)(e - s));
        const char *fe = t ? t : e;
        if (field >= 2) {
            double v = parse_num(s, fe);
            bool grp1 = field <= P.fcount + 1;
            if (v >= P.cmin1)      { if (grp1) f++;  else m++;  }
            else if (v <= P.cmin2) { if (grp1) f2++; else m2++; }
        }
        if (!t) break;
        s = t + 1;
        field++;
    }

    double p;
    bool cacheable = f < 65536 && f2 < 65536 && m < 65536 && m2 < 65536;
    uint64_t key = cacheable ? ((uint64_t)f << 48) | ((uint64_t)f2 << 32) | ((uint64_t)m << 16) | (uint64_t)m2 : 0;
    auto it = cacheable ? pc.m.find(key) : pc.m.end();
    if (it != pc.m.end()) {
        p = it->second;
    } else {
        p = fisher_p(lf, f, f2, m, m2, P.alt);
        if (cacheable) {
            if (pc.m.size() > 1000000) pc.m.clear();
            pc.m.emplace(key, p);
        }
    }
    if (P.use_pcut && !(p < P.pcut)) return;

    char buf[128];
    int k = std::snprintf(buf, sizeof buf, "%.6g\t%ld\t%ld\t%ld\t%ld\t", p, f, f2, m, m2);
    out.append(buf, k);
    out.append(b, (size_t)(e - b));
    out.push_back('\n');
}

// ---- pipeline state ----
static std::mutex mtx;
static std::condition_variable cv_in, cv_out;
static std::deque<Batch> in_q;
static std::map<size_t, Batch> results;
static bool reading_done = false;
static size_t total_batches = 0;
static size_t max_queue = 8;

using Clock = std::chrono::steady_clock;
static inline double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

// statistics (only reported with --stats)
static double st_read_wait = 0;     // reader: time inside fread (waiting for upstream)
static double st_push_wait = 0;     // reader: time blocked because all workers are busy
static double st_worker_idle = 0;   // sum over workers: time waiting for input
static size_t st_lines = 0, st_bytes = 0;

// Reader: reads big blocks from stdin, cuts them into batches of P.batch lines,
// hands them to the workers. Data is copied as little as possible.
static void reader(const Params &P) {
    const size_t CHUNK = 1 << 22;              // 4 MB per fread
    RawBuf cur;
    size_t start = 0;                          // start of the batch being assembled
    size_t scan = 0;                           // next byte to scan for '\n'
    size_t lines = 0;
    size_t seq = 0;

    auto push = [&](Batch &&bb) {
        auto t0 = Clock::now();
        std::unique_lock<std::mutex> lk(mtx);
        cv_in.wait(lk, [] { return in_q.size() < max_queue; });
        in_q.push_back(std::move(bb));
        cv_in.notify_all();
        st_push_wait += secs(t0, Clock::now());
    };

    while (true) {
        cur.reserve(cur.n + CHUNK);
        auto t0 = Clock::now();
        size_t n = std::fread(cur.p + cur.n, 1, CHUNK, stdin);
        st_read_wait += secs(t0, Clock::now());
        if (n == 0) break;
        cur.n += n;
        st_bytes += n;

        while (true) {
            const char *nl = (const char *)std::memchr(cur.p + scan, '\n', cur.n - scan);
            if (!nl) { scan = cur.n; break; }
            scan = (size_t)(nl - cur.p) + 1;
            if (++lines < P.batch) continue;

            // batch complete: bytes [start, scan)
            st_lines += lines;
            lines = 0;
            size_t bytes = scan - start, tail = cur.n - scan;
            Batch b;
            b.seq = seq++;
            if (start == 0 && bytes >= tail) {
                // big batch: hand over the whole buffer, copy only the small tail
                RawBuf nb;
                nb.reserve(tail + CHUNK);
                if (tail) std::memcpy(nb.p, cur.p + scan, tail);
                nb.n = tail;
                cur.n = scan;
                b.data = std::move(cur);
                cur = std::move(nb);
                start = 0;
                scan = 0;
            } else {
                // small batch: copy it out, keep working in the same buffer
                b.data.reserve(bytes);
                std::memcpy(b.data.p, cur.p + start, bytes);
                b.data.n = bytes;
                start = scan;
            }
            push(std::move(b));
        }
        // drop what has been handed out, keep the unfinished tail
        if (start > 0) {
            std::memmove(cur.p, cur.p + start, cur.n - start);
            cur.n -= start;
            scan -= start;
            start = 0;
        }
    }

    // EOF: whatever is left is a final (possibly unterminated) batch
    if (cur.n > 0) {
        if (cur.p[cur.n - 1] != '\n') { cur.reserve(cur.n + 1); cur.p[cur.n++] = '\n'; }
        for (size_t i = 0; i < cur.n; i++) st_lines += (cur.p[i] == '\n');
        Batch b;
        b.seq = seq++;
        b.data = std::move(cur);
        push(std::move(b));
    }
    {
        std::lock_guard<std::mutex> lk(mtx);
        total_batches = seq;
        reading_done = true;
    }
    cv_in.notify_all();
    cv_out.notify_all();
}

static void worker(const Params &P) {
    LogFact lf;
    PCache pc;
    double idle = 0;
    while (true) {
        Batch b;
        {
            auto t0 = Clock::now();
            std::unique_lock<std::mutex> lk(mtx);
            cv_in.wait(lk, [] { return !in_q.empty() || reading_done; });
            idle += secs(t0, Clock::now());
            if (in_q.empty()) { st_worker_idle += idle; return; }
            b = std::move(in_q.front());
            in_q.pop_front();
            cv_in.notify_all();
        }
        const char *p = b.data.p;
        const char *end = p + b.data.n;
        while (p < end) {
            const char *nl = (const char *)std::memchr(p, '\n', (size_t)(end - p));
            if (!nl) nl = end;
            process_line(p, nl, P, lf, pc, b.out);
            p = nl + 1;
        }
        b.data.release();                      // free the input block early
        {
            std::lock_guard<std::mutex> lk(mtx);
            size_t s = b.seq;
            results.emplace(s, std::move(b));
        }
        cv_out.notify_all();
    }
}

static void usage(const char *prog) {
    std::fprintf(stderr,
        "Usage: %s --fcount N [options] < table.tsv > out.tsv\n"
        "  -f, --fcount N    number of group-1 columns (columns 2..N+1)  [required]\n"
        "  -1, --cmin1 X     value >= X counts as 'high'                 [default 1]\n"
        "  -2, --cmin2 Y     value <= Y counts as 'low'                  [default 0]\n"
        "  -H, --header      first input line is a header: copy it to the output with\n"
        "                    the columns pval, grp1_high, grp1_low, grp2_high, grp2_low added\n"
        "  -p, --pcut P      only output lines with p-value < P\n"
        "  -a, --alt A       test type: two (default), greater, less, point\n"
        "      --two-sided   two-sided Fisher exact test (default)\n"
        "      --one-sided   one-sided test, upper tail P(X>=a) (same as --greater)\n"
        "      --greater     one-sided, upper tail (enrichment in group 1)\n"
        "      --less        one-sided, lower tail (depletion in group 1)\n"
        "      --point       point probability P(X=a) only (not a p-value)\n"
        "      --stats       print a throughput/bottleneck report to stderr at the end\n"
        "  -t, --threads T   worker threads                              [default: 80%% of cores]\n"
        "  -b, --batch B     lines per thread buffer                     [default 1000000]\n",
        prog);
}

int main(int argc, char **argv) {
    Params P;
    static option opts[] = {
        {"fcount", required_argument, nullptr, 'f'},
        {"cmin1", required_argument, nullptr, '1'},
        {"cmin2", required_argument, nullptr, '2'},
        {"pcut", required_argument, nullptr, 'p'},
        {"header", no_argument, nullptr, 'H'},
        {"alt", required_argument, nullptr, 'a'},
        {"two-sided", no_argument, nullptr, 1000},
        {"one-sided", no_argument, nullptr, 1001},
        {"greater", no_argument, nullptr, 1001},
        {"less", no_argument, nullptr, 1002},
        {"point", no_argument, nullptr, 1003},
        {"threads", required_argument, nullptr, 't'},
        {"batch", required_argument, nullptr, 'b'},
        {"stats", no_argument, nullptr, 1004},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "f:1:2:p:Ha:t:b:h", opts, nullptr)) != -1) {
        switch (c) {
        case 'f': P.fcount = std::atol(optarg); break;
        case '1': P.cmin1 = std::atof(optarg); break;
        case '2': P.cmin2 = std::atof(optarg); break;
        case 'p': P.pcut = std::atof(optarg); P.use_pcut = true; break;
        case 'H': P.header = true; break;
        case 'a': {
            std::string v = optarg;
            if (v == "two" || v == "two-sided") P.alt = 0;
            else if (v == "greater" || v == "one-sided") P.alt = 2;
            else if (v == "less") P.alt = 1;
            else if (v == "point") P.alt = 3;
            else { usage(argv[0]); return 1; }
            break;
        }
        case 1000: P.alt = 0; break;
        case 1001: P.alt = 2; break;
        case 1002: P.alt = 1; break;
        case 1003: P.alt = 3; break;
        case 1004: P.stats = true; break;
        case 't': P.threads = (unsigned)std::max(1L, std::atol(optarg)); break;
        case 'b': P.batch = (size_t)std::max(1L, std::atol(optarg)); break;
        default: usage(argv[0]); return c == 'h' ? 0 : 1;
        }
    }
    if (P.fcount < 0) { usage(argv[0]); return 1; }
    if (P.threads == 0) {
        unsigned hc = std::max(1u, std::thread::hardware_concurrency());
        P.threads = std::max(1u, (unsigned)(hc * 0.8));
    }
    max_queue = P.threads;   // bounds memory: at most ~2*threads batches in flight

#ifdef F_SETPIPE_SZ
    // Bigger pipe buffers (default 64 KB) let the upstream decompressor run
    // ahead instead of stalling on every small write. Harmless if not a pipe.
    fcntl(0, F_SETPIPE_SZ, 1 << 20);
    fcntl(1, F_SETPIPE_SZ, 1 << 20);
#endif

    if (P.header) {
        char *hdr = nullptr;
        size_t cap = 0;
        ssize_t len = getline(&hdr, &cap, stdin);
        if (len > 0) {
            if (hdr[len - 1] == '\n') len--;
            std::fputs("pval\tgrp1_high\tgrp1_low\tgrp2_high\tgrp2_low\t", stdout);
            std::fwrite(hdr, 1, (size_t)len, stdout);
            std::fputc('\n', stdout);
        }
        std::free(hdr);
    }

    auto t_start = Clock::now();
    std::thread rd(reader, std::cref(P));
    std::vector<std::thread> pool;
    for (unsigned i = 0; i < P.threads; i++) pool.emplace_back(worker, std::cref(P));

    // Writer (main thread): emit batches in original order.
    size_t next = 0;
    while (true) {
        Batch b;
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv_out.wait(lk, [&] {
                return results.count(next) || (reading_done && next >= total_batches);
            });
            auto it = results.find(next);
            if (it == results.end()) break;
            b = std::move(it->second);
            results.erase(it);
        }
        std::fwrite(b.out.data(), 1, b.out.size(), stdout);
        next++;
    }

    rd.join();
    for (auto &t : pool) t.join();
    std::fflush(stdout);

    if (P.stats) {
        double wall = secs(t_start, Clock::now());
        double reader_busy = std::max(0.0, wall - st_read_wait - st_push_wait);
        std::fprintf(stderr,
            "--- fisher_mt stats ---\n"
            "threads %u | %.0f lines | %.1f MB | %.2f s wall | %.1f MB/s | %.2f M lines/s\n"
            "reader : %.0f%% waiting for input | %.0f%% blocked (workers all busy) | %.0f%% own work\n"
            "workers: %.0f%% idle (waiting for input), on average\n",
            P.threads, (double)st_lines, st_bytes / 1e6, wall, st_bytes / 1e6 / wall,
            st_lines / 1e6 / wall,
            100 * st_read_wait / wall, 100 * st_push_wait / wall, 100 * reader_busy / wall,
            100 * st_worker_idle / (wall * P.threads));
        std::fprintf(stderr,
            "how to read this: workers mostly idle + reader mostly waiting for input -> upstream "
            "(decompressor/pipe) is the limit;\n"
            "                  workers mostly idle + reader 'own work' high -> reader is the limit;\n"
            "                  reader mostly blocked -> workers are the limit (raise -t).\n");
    }
    return 0;
}
