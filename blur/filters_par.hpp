#pragma once
#include "matrix.hpp"

namespace Filter {
    // calling the parallel blur function.
    Matrix blur_par(Matrix m, int radius, int num_threads);
}
