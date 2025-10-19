/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>    // for std::memcpy
#include <algorithm>  // for std::copy_n
#include <numeric>    // optional

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
    if (size > 0 && other.data) {
        // using memcpy as it is faster than initial element by element copy inside the loop.
        std::memcpy(data, other.data, sizeof(double) * size);
    }
    //for (auto i{0}; i < size; i++)
    //{
      //  data[i] = other.data[i];
    //}
}


// Reassign resources forcefully using the move constructor which is cheaper
Vector::Vector(Vector&& other) noexcept
    : size(other.size), data(other.data)
{
    other.size = 0;
    other.data = nullptr;
}

// Move definition for assignment operator such that it frees previous unused resources.
Vector& Vector::operator=(const Vector& other)
{
    if (this == &other) return *this;

    if (size != other.size) {
        delete[] data;
        size = other.size;
        data = (size > 0) ? new double[size] : nullptr;
    }
    if (size > 0 && other.data) {
        std::memcpy(data, other.data, sizeof(double) * size);
    }
    return *this;
}

// Move cosntructor assignment
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
//unigned not auto
    for (unsigned i=0; i < size; i++)
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
//Vector Vector::operator/(double div) 

Vector Vector::operator/(double div) const
{
    //auto result{*this};
    Vector result{size};
    for (auto i{0u}; i < size; i++)
    {
        //result[i] /= div;
        result[i] = data[i] / div;
    }

    return result;
}
//Vector Vector::operator-(double sub)

Vector Vector::operator-(double sub) const
{
    //auto result{*this};
    Vector result{size};
    for (auto i{0u}; i < size; i++)
    {
        //result[i] -= sub
        result[i] = data[i] - sub;
            
    }

    return result;
}
//double Vector::dot(constVector rhs) const
double Vector::dot(const Vector& rhs) const
{
    double result{0};

    for (auto i{0u}; i < size; i++)
    {
        result += data[i] * rhs[i];
    }

    return result;
}
