/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
// Large object is passed by const reference and avoids copying the entire vector everytime
std::vector<double> correlation_coefficients(const std::vector<Vector>& datasets);
double pearson(const Vector& vec1, const Vector& vec2);
};

#endif
