/*
Author: David Holmqvist <daae19@student.bth.se>
*/
#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

namespace Filter {

    namespace Gauss {
        void get_weights(int n, double* weights_out) {
            for (int i = 0; i <= n; ++i) {
                double x = static_cast<double>(i) * max_x / n;
                weights_out[i] = std::exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius) {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};

        // Cache as signed to avoid unsigned promotions
        const int W = static_cast<int>(dst.get_x_size());
        const int H = static_cast<int>(dst.get_y_size());

        // Precompute weights once, reuse in both passes
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

        // ---- Horizontal pass (x outer, y inner to match baseline) ----
        for (int x = 0; x < W; ++x) {
            for (int y = 0; y < H; ++y) {
                double r = w[0] * dst.r(x, y);
                double g = w[0] * dst.g(x, y);
                double b = w[0] * dst.b(x, y);
                double n = w[0];

                for (int wi = 1; wi <= radius; ++wi) {
                    const double wc = w[wi];

                    int x2 = x - wi;
                    if (x2 >= 0) {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < W) {
                        r += wc * dst.r(x2, y);
                        g += wc * dst.g(x2, y);
                        b += wc * dst.b(x2, y);
                        n += wc;
                    }
                }
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        // ---- Vertical pass (x outer, y inner to match baseline) ----
        for (int x = 0; x < W; ++x) {
            for (int y = 0; y < H; ++y) {
                double r = w[0] * scratch.r(x, y);
                double g = w[0] * scratch.g(x, y);
                double b = w[0] * scratch.b(x, y);
                double n = w[0];

                for (int wi = 1; wi <= radius; ++wi) {
                    const double wc = w[wi];

                    int y2 = y - wi;
                    if (y2 >= 0) {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < H) {
                        r += wc * scratch.r(x, y2);
                        g += wc * scratch.g(x, y2);
                        b += wc * scratch.b(x, y2);
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }

        return dst;
    }
}
