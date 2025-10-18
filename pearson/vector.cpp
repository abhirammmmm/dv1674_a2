/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

Vector::Vector()
    : size{0}, data{nullptr}
{
}

Vector::~Vector()
{
    if (data)
    {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size{size}, data{new double[size]}
{
}

Vector::Vector(unsigned size, double *data)
    : size{size}, data{data}
{
}

Vector::Vector(const Vector &other)
    : Vector{other.size}
{
//    for (auto i{0}; i < size; i++)
//    {
//        data[i] = other.data[i];
//    }

// Faster than initial element by element copy inside the loop.
    std::memcpy(data, other.data, size(double) * size);
}

// Reassign resources forcefully using the move constructor which is cheaper.
Vector::Vector(Vector&& other) noexcept
    : size(other.size), data(other.data)
{
    other.size = 0;
    other.data = nullptr;
}

// Move definition for assignment operator such that it frees previous unused resources.
Vector& Vector::operator=(Vector&& other) noexcept
{
    if (this == &other) return *this;

    delete[] data; // free current
    size = other.size;
    data = other.data;

    other.size = 0;
    other.data = nullptr;
    return *this;
}

unsigned Vector::get_size() const
{
    return size;
}

double *Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double &Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    double sum{0};

    for (auto i{0}; i < size; i++)
    {
        sum += data[i];
    }

    return sum / static_cast<double>(size);
}

double Vector::magnitude() const
{
    auto dot_prod{dot(*this)};
    return std::sqrt(dot_prod);
}

Vector Vector::operator/(double div) const
{
//    auto result{*this};
    Vector result{size};
    for (auto i{0}; i < size; i++)
    {
        //result[i] /= div;
        result[i] = data[i] / div;

    }

    return result;
}

Vector Vector::operator-(double sub) const
{
//    auto result{*this};
    Vector result{size};
    for (auto i{0}; i < size; i++)
    {
        //result[i] -= sub;
        result[i] = data[i] - sub;
    }

    return result;
}

double Vector::dot(const Vector& rhs) const
{
    double result{0};

    for (auto i{0}; i < size; i++)
    {
        result += data[i] * rhs[i];
    }

    return result;
}
