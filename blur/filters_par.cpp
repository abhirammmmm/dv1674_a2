#include "filters_par.hpp"
#include "matrix.hpp"
#include <pthread.h>
#include <vector>
#include <cmath>
#include <algorithm>

// used the same code of filters.cpp and added pthreads implementation.
namespace Filter {

namespace Gauss {
    constexpr unsigned max_radius{1000};
    constexpr float max_x{1.33f};
    constexpr float pi{3.14159f};

    inline void get_weights(int n, double* weights_out)
    {
        for (int i = 0; i <= n; ++i) {
            double x = static_cast<double>(i) * max_x / static_cast<double>(n);
            weights_out[i] = std::exp(-x * x * pi);
        }
    }
}
//seperating the horizontal operation and vertical operations parellel meaning assigning first nested loop to one worker and the next to another worker.
//passing first nested loop for horizontal to first worker.
struct Pass1Args {
    const Matrix* src;    
    Matrix* scratch;      
    const double* w;      
    int radius;
    unsigned W, H;
    unsigned y0, y1;      
};

// first worker implementation
static void* pass1_worker(void* p) {
    auto* a = static_cast<Pass1Args*>(p);
    const Matrix& src = *a->src;
    Matrix& scratch   = *a->scratch;
    const double* w   = a->w;
    const int R       = a->radius;
    const unsigned W  = a->W;

    for (unsigned y = a->y0; y < a->y1; ++y) {
        for (unsigned x = 0; x < W; ++x) {
            double r = w[0] * src.r(x, y);
            double g = w[0] * src.g(x, y);
            double b = w[0] * src.b(x, y);
            double n = w[0];

            for (int wi = 1; wi <= R; ++wi) {
                const double wc = w[wi];

                if (x >= static_cast<unsigned>(wi)) {
                    const unsigned xl = x - static_cast<unsigned>(wi);
                    r += wc * src.r(xl, y);
                    g += wc * src.g(xl, y);
                    b += wc * src.b(xl, y);
                    n += wc;
                }
                
                const unsigned xr = x + static_cast<unsigned>(wi);
                if (xr < a->W) {
                    r += wc * src.r(xr, y);
                    g += wc * src.g(xr, y);
                    b += wc * src.b(xr, y);
                    n += wc;
                }
            }

            scratch.r(x, y) = static_cast<unsigned char>(r / n);
            scratch.g(x, y) = static_cast<unsigned char>(g / n);
            scratch.b(x, y) = static_cast<unsigned char>(b / n);
        }
    }
    return nullptr;
}

//passing second nested loop for vertical to second worker.
struct Pass2Args {
    const Matrix* scratch; 
    Matrix* dst;           
    const double* w;
    int radius;
    unsigned W, H;
    unsigned y0, y1;
};

// second worker implementation
static void* pass2_worker(void* p) {
    auto* a = static_cast<Pass2Args*>(p);
    const Matrix& scratch = *a->scratch;
    Matrix& dst           = *a->dst;
    const double* w       = a->w;
    const int R           = a->radius;
    const unsigned H      = a->H;

    for (unsigned y = a->y0; y < a->y1; ++y) {
        for (unsigned x = 0; x < a->W; ++x) {
            double r = w[0] * scratch.r(x, y);
            double g = w[0] * scratch.g(x, y);
            double b = w[0] * scratch.b(x, y);
            double n = w[0];

            for (int wi = 1; wi <= R; ++wi) {
                const double wc = w[wi];

                if (y >= static_cast<unsigned>(wi)) {
                    const unsigned yu = y - static_cast<unsigned>(wi);
                    r += wc * scratch.r(x, yu);
                    g += wc * scratch.g(x, yu);
                    b += wc * scratch.b(x, yu);
                    n += wc;
                }

                const unsigned yd = y + static_cast<unsigned>(wi);
                if (yd < H) {
                    r += wc * scratch.r(x, yd);
                    g += wc * scratch.g(x, yd);
                    b += wc * scratch.b(x, yd);
                    n += wc;
                }
            }

            dst.r(x, y) = static_cast<unsigned char>(r / n);
            dst.g(x, y) = static_cast<unsigned char>(g / n);
            dst.b(x, y) = static_cast<unsigned char>(b / n);
        }
    }
    return nullptr;
}

// parellal implementation of blur
Matrix blur_par(Matrix m, int radius_in, int num_threads)
{
    // verifying radius once
    int radius = std::max(0, std::min<int>(radius_in, static_cast<int>(Gauss::max_radius) - 1));

    // Cache size
    const unsigned W = m.get_x_size();
    const unsigned H = m.get_y_size();
    if (W == 0 || H == 0 || radius == 0) {
        return m;
    }

    // verrifying the  thread counts
    if (num_threads < 1) num_threads = 1;
    if ((unsigned)num_threads > H) num_threads = (int)H;

    // Precompute Gaussian weights once
    std::vector<double> w((std::size_t)radius + 1);
    Gauss::get_weights(radius, w.data());

    // will have the output image
    Matrix dst{ m };                  
    Matrix scratch(W, H, m.get_color_max());

    //threads implementation - divinding the image rows in horizontal among the threads
    {
        std::vector<pthread_t> th((std::size_t)num_threads);
        std::vector<Pass1Args> args((std::size_t)num_threads);
        // dividing the total image rows equally among the threads, if not perfect division the remainder rows to first threads.
        unsigned rows_per = H / (unsigned)num_threads;
        unsigned rem = H % (unsigned)num_threads;
        unsigned y = 0;

        for (int t = 0; t < num_threads; ++t) {
            unsigned take = rows_per + (rem ? 1u : 0u);
            if (rem) --rem;

            args[(std::size_t)t] = Pass1Args{
                &dst, &scratch,
                w.data(), radius,
                W, H,
                y, y + take
            };
            //create the thread and start executing the worker 1 with its own range
            pthread_create(&th[(std::size_t)t], nullptr, pass1_worker, &args[(std::size_t)t]);
            y += take;
        }
        for (int t = 0; t < num_threads; ++t) {
            //waiting till thread finishes entire rows and then joining them
            pthread_join(th[(std::size_t)t], nullptr);
        }
    }

    //threads implementation - divinding the image rows in vertically among the threads
    {
        std::vector<pthread_t> th((std::size_t)num_threads);
        std::vector<Pass2Args> args((std::size_t)num_threads);
        // dividing the total image rows equally among the threads, if not perfect division the remainder rows to first threads.

        unsigned rows_per = H / (unsigned)num_threads;
        unsigned rem = H % (unsigned)num_threads;
        unsigned y = 0;

        for (int t = 0; t < num_threads; ++t) {
            unsigned take = rows_per + (rem ? 1u : 0u);
            if (rem) --rem;

            args[(std::size_t)t] = Pass2Args{
                &scratch, &dst,
                w.data(), radius,
                W, H,
                y, y + take
            };
            //create the thread and start executing the worker 2 with its own range
            pthread_create(&th[(std::size_t)t], nullptr, pass2_worker, &args[(std::size_t)t]);
            y += take;
        }
        for (int t = 0; t < num_threads; ++t) {
            //waiting till thread finishes entire rows and then joining them
            pthread_join(th[(std::size_t)t], nullptr);
        }
    }

    return dst;
}

}
