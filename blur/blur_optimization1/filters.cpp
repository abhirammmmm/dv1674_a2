/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace Filter
{
    namespace Gauss
    {
        void get_weights(int n, double* weights_out)
        {
            for (int i = 0; i <= n; ++i)
            {
                double x = static_cast<double>(i) * max_x / static_cast<double>(n);
                weights_out[i] = std::exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        // caching the width and height once by calling the get_x_size() and get_y_size() functions once before the both nested loops to eliminate unnecessary calling of the functions.
        const unsigned W = m.get_x_size();
        const unsigned H = m.get_y_size();

        // Precomputing the Gaussian weights once rather than per pixel which is very expensive due to nested loops. The exp() function calls reduces significantly.
        std::vector<double> w(static_cast<std::size_t>(radius) + 1);
        Gauss::get_weights(radius, w.data());

        Matrix scratch{PPM::max_dimension};

        auto dst = m;

        for (unsigned x = 0; x < W; x++)
        {
            for (unsigned y = 0; y < H; y++)
            {
                double r = w[0] * dst.r(x, y);
                double g = w[0] * dst.g(x, y);
                double b = w[0] * dst.b(x, y);
                double n = w[0];

                for (int wi = 1; wi <= radius; ++wi)
                {
                    const double wc = w[wi];

                    // static casting
                    if (x >= static_cast<unsigned>(wi)) {
                        const unsigned xl = x - static_cast<unsigned>(wi);
                        r += wc * dst.r(xl, y);
                        g += wc * dst.g(xl, y);
                        b += wc * dst.b(xl, y);
                        n += wc;
                    }
                    // static casting
                    const unsigned xr = x + static_cast<unsigned>(wi);
                    if (xr < W) {
                        r += wc * dst.r(xr, y);
                        g += wc * dst.g(xr, y);
                        b += wc * dst.b(xr, y);
                        n += wc;
                    }
                }

                scratch.r(x, y) = static_cast<unsigned char>(r / n);
                scratch.g(x, y) = static_cast<unsigned char>(g / n);
                scratch.b(x, y) = static_cast<unsigned char>(b / n);
            }
        }

        for (unsigned x = 0; x < W; x++)
        {
            for (unsigned y = 0; y < H; y++)
            {
                double r = w[0] * scratch.r(x, y);
                double g = w[0] * scratch.g(x, y);
                double b = w[0] * scratch.b(x, y);
                double n = w[0];

                for (int wi = 1; wi <= radius; ++wi)
                {
                    const double wc = w[wi];

                    // static casting
                    if (y >= static_cast<unsigned>(wi)) {
                        const unsigned yu = y - static_cast<unsigned>(wi);
                        r += wc * scratch.r(x, yu);
                        g += wc * scratch.g(x, yu);
                        b += wc * scratch.b(x, yu);
                        n += wc;
                    }
                    // static casting
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

        return dst;
    }

}
