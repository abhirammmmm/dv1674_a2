/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

namespace Analysis {

std::vector<double> correlation_coefficients(const std::vector<Vector>& datasets)
{
    std::vector<double> result {};


    // Reserve capacity for number of unique pairs before handto reduce computations.
    auto n = datasets.size();
    if (n > 1) result.reserve(n * (n - 1) / 2);

    for (auto sample1 { 0u }; sample1 < datasets.size() - 1; sample1++) {
        for (auto sample2 { sample1 + 1u }; sample2 < datasets.size(); ++sample2) {
            auto corr { pearson(datasets[sample1], datasets[sample2]) };
            result.push_back(corr);
        }
    }

    return result;
}

double pearson(const Vector& vec1, const Vector& vec2)
{
    auto x_mean { vec1.mean() };
    auto y_mean { vec2.mean() };

    auto x_mm { vec1 - x_mean };
    auto y_mm { vec2 - y_mean };

    auto x_mag { x_mm.magnitude() };
    auto y_mag { y_mm.magnitude() };

     // To protect against division by zero
    if (x_mag == 0.0 || y_mag == 0.0) {
        return 0.0;
    }
    
    auto x_mm_over_x_mag { x_mm / x_mag };
    auto y_mm_over_y_mag { y_mm / y_mag };

    auto r { x_mm_over_x_mag.dot(y_mm_over_y_mag) };

    return std::max(std::min(r, 1.0), -1.0);
}
};
