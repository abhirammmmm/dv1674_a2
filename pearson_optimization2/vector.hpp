/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#if !defined(VECTOR_HPP)
#define VECTOR_HPP
#include <cstddef>

class Vector {
private:
    unsigned size;
    double* data;

public:
    Vector();
    Vector(unsigned size);
    Vector(unsigned size, double* data);
    Vector(const Vector& other);
    // Defining new move constructor to avoid multiple double copy operations
    Vector(Vector&& other) noexcept;
    // Defining "=" operator overload to cut runtime with lower deep copies and lower heap allocations. 
    Vector& operator=(const Vector& other);       // copy assign
    Vector& operator=(Vector&& other) noexcept;   // move assign
    ~Vector();

    double magnitude() const;
    double mean() const;
    double normalize() const;
    // adding const reference for passing vector
    double dot(const Vector& rhs) const;

    unsigned get_size() const;
    double* get_data();
//adding const
    Vector operator/(double div) const;
    Vector operator-(double sub) const;
    double operator[](unsigned i) const;
    double& operator[](unsigned i);
};

#endif
