
#include <pthread.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include "vector.hpp"
#include "dataset.hpp"

struct ThreadCtx {
    unsigned tid;
    unsigned T;
    unsigned N;
    const std::vector<Vector>* datasets;
    const std::vector<double>* means;
    const std::vector<double>* invstds;
    const std::vector<unsigned>* lens;
    const std::vector<size_t>* row_offset; // where row i starts in the matrix
    std::vector<double>* results;          // continuous order for matrix
};

// const pointers for kernel operation
static inline double pearson_pair_kernel_idx(
    const Vector& a, const Vector& b,
    unsigned n,
    const double mean_a, const double invstd_a,
    const double mean_b, const double invstd_b,
    const double invN)
{
    double acc = 0.0;
    unsigned k = 0;
    const unsigned n4 = n & ~3u;
    for (; k < n4; k += 4) {
        const double a0 = (a[k+0] - mean_a) * invstd_a; const double b0 = (b[k+0] - mean_b) * invstd_b;
        const double a1 = (a[k+1] - mean_a) * invstd_a; const double b1 = (b[k+1] - mean_b) * invstd_b;
        const double a2 = (a[k+2] - mean_a) * invstd_a; const double b2 = (b[k+2] - mean_b) * invstd_b;
        const double a3 = (a[k+3] - mean_a) * invstd_a; const double b3 = (b[k+3] - mean_b) * invstd_b;
        acc += a0*b0 + a1*b1 + a2*b2 + a3*b3;
    }
    for (; k < n; ++k) {
        const double da = (a[k] - mean_a) * invstd_a;
        const double db = (b[k] - mean_b) * invstd_b;
        acc += da * db;
    }
    return acc * invN;
}

// population mean for the assigned vector
static inline double compute_mean_vec(const Vector& v) {
    const unsigned n = v.get_size();
    if (n == 0) return 0.0;
    double s = 0.0;
    for (unsigned k = 0; k < n; ++k) s += v[k];
    return s / static_cast<double>(n);
}

// splitting the correlation calculation such that the calculations are done only once snice upper and lower triangles are interchangable.
// precomputing each vectors mean and SD to avoid rehashing per pair
static inline double compute_invstd_from_mean_vec(const Vector& v, double mean) {
    const unsigned n = v.get_size();
    if (n == 0) return 0.0;
    double sumsq = 0.0;
    for (unsigned k = 0; k < n; ++k) {
        const double d = v[k] - mean;
        sumsq += d * d;
    }
    const double var = sumsq / static_cast<double>(n);
    if (var <= 0.0) return 0.0;
    return 1.0 / std::sqrt(var);
}
//flatten the matrix order with a single loop assuming constants are already available
static void* worker(void* arg) {
    ThreadCtx* C = reinterpret_cast<ThreadCtx*>(arg);
    const unsigned N   = C->N;
    const unsigned T   = C->T;
    const unsigned tid = C->tid;

    const auto& ds   = *C->datasets;
    const auto& mu   = *C->means;
    const auto& invs = *C->invstds;
    const auto& len  = *C->lens;
    const auto& roff = *C->row_offset;
    auto& out        = *C->results;

    // Give each thread a contiguous block of rows while splitting rows
    const unsigned rows_per = (N + T - 1) / T;
    const unsigned i_begin  = tid * rows_per;
    const unsigned i_end    = std::min(N, i_begin + rows_per);

    for (unsigned i = i_begin; i < i_end; ++i) {
        const unsigned ni = len[i];
        const double mi   = mu[i];
        const double invsi= invs[i];
        const bool ok_i   = (ni > 0 && invsi != 0.0);

        // Starting index in the flattened results for row i as base'
        const size_t base = roff[i];

        if (i + 1 >= N) continue; // last row has length 0

        size_t w = 0; // write cursor within this row's segment
        for (unsigned j = i + 1; j < N; ++j, ++w) {
            double rij = 0.0;

            const unsigned nj   = len[j];
            const double   mj   = mu[j];
            const double   invsj= invs[j];
            const bool     ok_j = (nj > 0 && invsj != 0.0);

            if (ok_i && ok_j && ni == nj) {
                const double invN = 1.0 / static_cast<double>(ni);
                rij = pearson_pair_kernel_idx(ds[i], ds[j], ni, mi, invsi, mj, invsj, invN);
            } else {
                rij = 0.0;
            }

            out[base + w] = rij;
        }
    }
    return nullptr;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [num_threads]\n";
        std::exit(1);
    }

    const char* infile  = argv[1];
    const char* outfile = argv[2];
    int num_threads     = std::atoi(argv[3]);
    if (num_threads <= 0) {
        std::cerr << "num_threads must be > 0\n";
        return 1;
    }
    const unsigned T = static_cast<unsigned>(num_threads);

    // Reading the datasets in order
    std::vector<Vector> datasets = Dataset::read(std::string(infile));
    const unsigned N = static_cast<unsigned>(datasets.size());
    if (N < 2) {
        std::vector<double> empty;
        Dataset::write(empty, std::string(outfile));
        return 0;
    }

    // Precompute per-vector stats once
    std::vector<double>   means(N, 0.0), invstds(N, 0.0);
    std::vector<unsigned> lens(N, 0);
    for (unsigned i = 0; i < N; ++i) {
        const unsigned n = datasets[i].get_size();
        lens[i] = n;
        const double m = compute_mean_vec(datasets[i]);
        means[i]   = m;
        invstds[i] = compute_invstd_from_mean_vec(datasets[i], m);
    }

    // Total pairs and flattened results in the SAME order as sequential
    const size_t total_pairs = static_cast<size_t>(N) * static_cast<size_t>(N - 1) / 2;
    std::vector<double> results(total_pairs, 0.0);

    // Precompute where each row i starts in the flattened array
    std::vector<size_t> row_offset(N, 0);
    for (unsigned i = 1; i < N; ++i) {
        row_offset[i] = row_offset[i-1] + static_cast<size_t>((N - 1) - (i - 1));
    }

    // launching threads
    std::vector<pthread_t> threads(T);
    std::vector<ThreadCtx> ctx(T);
    for (unsigned t = 0; t < T; ++t) {
        ctx[t] = ThreadCtx{t, T, N, &datasets, &means, &invstds, &lens, &row_offset, &results};
        const int rc = pthread_create(&threads[t], nullptr, worker, &ctx[t]);
        if (rc != 0) {
            std::cerr << "pthread_create failed: " << rc << "\n";
            std::exit(1);
        }
    }
    for (unsigned t = 0; t < T; ++t) {
        pthread_join(threads[t], nullptr);
    }

    // writing results for same half of the computed matrix flattened upper-triangle order as sequential
    Dataset::write(results, std::string(outfile));
    return 0;
}
