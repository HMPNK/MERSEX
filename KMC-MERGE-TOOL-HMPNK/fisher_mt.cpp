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
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <deque>
#include <getopt.h>
#include <iostream>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

struct Params {
    double cmin1 = 1;
    double cmin2 = 0;
    int    alt = 0;            // 0 = two-sided, 1 = less, 2 = greater, 3 = point probability
    long   fcount = -1;
    bool   header = false;
    bool   use_pcut = false;
    double pcut = 1.0;
    unsigned threads = 0;      // 0 = 80% of hardware concurrency
    size_t batch = 1000000;    // lines per work unit
};

struct Batch {
    size_t seq = 0;
    std::vector<std::string> lines;
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

static void process_line(const std::string &line, const Params &P, LogFact &lf, PCache &pc, std::string &out) {
    long f = 0, f2 = 0, m = 0, m2 = 0;
    long field = 1;                       // 1-based field index, like AWK
    size_t pos = 0, len = line.size();
    while (true) {
        size_t end = line.find('\t', pos);
        if (end == std::string::npos) end = len;
        if (field >= 2) {
            double v = std::strtod(std::string(line, pos, end - pos).c_str(), nullptr);
            bool grp1 = field <= P.fcount + 1;
            if (v >= P.cmin1)      { if (grp1) f++;  else m++;  }
            else if (v <= P.cmin2) { if (grp1) f2++; else m2++; }
        }
        if (end >= len) break;
        pos = end + 1;
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
    out.append(line);
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

static void reader(const Params &P) {
    Batch b;
    size_t seq = 0;
    std::string line;
    auto push = [&](Batch &&bb) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_in.wait(lk, [] { return in_q.size() < max_queue; });
        in_q.push_back(std::move(bb));
        cv_in.notify_all();
    };
    b.seq = seq;
    while (std::getline(std::cin, line)) {
        b.lines.push_back(std::move(line));
        if (b.lines.size() >= P.batch) {
            push(std::move(b));
            b = Batch();
            b.seq = ++seq;
        }
    }
    {
        std::unique_lock<std::mutex> lk(mtx);
        if (!b.lines.empty()) {
            cv_in.wait(lk, [] { return in_q.size() < max_queue; });
            in_q.push_back(std::move(b));
            seq++;
        }
        total_batches = seq;
        reading_done = true;
    }
    cv_in.notify_all();
    cv_out.notify_all();
}

static void worker(const Params &P) {
    LogFact lf;
    PCache pc;
    while (true) {
        Batch b;
        {
            std::unique_lock<std::mutex> lk(mtx);
            cv_in.wait(lk, [] { return !in_q.empty() || reading_done; });
            if (in_q.empty()) return;
            b = std::move(in_q.front());
            in_q.pop_front();
            cv_in.notify_all();
        }
        b.out.reserve(b.lines.size() * 64);
        for (const auto &l : b.lines) process_line(l, P, lf, pc, b.out);
        b.lines.clear();
        b.lines.shrink_to_fit();
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

    std::ios::sync_with_stdio(false);
    if (P.header) {
        std::string hdr;
        if (std::getline(std::cin, hdr)) {
            std::fputs("pval\tgrp1_high\tgrp1_low\tgrp2_high\tgrp2_low\t", stdout);
            std::fwrite(hdr.data(), 1, hdr.size(), stdout);
            std::fputc('\n', stdout);
        }
    }

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
    return 0;
}
