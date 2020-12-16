/*
 * double.hpp
 *
 *  Created on: Sep 9, 2016
 *      Author: akozlins
 */

#ifndef MU3E_UTIL_DOUBLE_HPP_
#define MU3E_UTIL_DOUBLE_HPP_

#include "float.hh"

struct double2 {
    double x, y;
};

#include <CLHEP/Vector/ThreeVector.h>

struct double3 {
    double x, y, z;

    operator CLHEP::Hep3Vector () const {
        return { x, y, z };
    }
};


inline
double2 make_double2(const double3& d3) {
    return { d3.x, d3.y };
}


inline
double3 make_double3(const float3& f3) {
    return { (double)f3.x, (double)f3.y, (double)f3.z };
}


inline
float3 make_float3(const double3& d3) {
    return { (float)d3.x, (float)d3.y, (float)d3.z };
}

//
// add
//


inline
double2 operator + (const double2& a, const double2& b) {
    return { a.x + b.x, a.y + b.y };
}


inline
double2& operator += (double2& a, const double2& b) {
    a.x += b.x; a.y += b.y;
    return a;
}


inline
double3 operator + (const double3& a, const double3& b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}


inline
double3& operator += (double3& a, const double3& b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

//
// sub
//


inline
double2 operator - (const double2& a, const double2& b) {
    return { a.x - b.x, a.y - b.y };
}


inline
double2& operator -= (double2& a, const double2& b) {
    a.x -= b.x; a.y -= b.y;
    return a;
}


inline
double3 operator - (const double3& a, const double3& b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}


inline
double3& operator -= (double3& a, const double3& b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

//
// mul
//


inline
double2 operator * (const double2& a, double b) {
    return { a.x * b, a.y * b };
}


inline
double2 operator * (double b, const double2& a) {
    return { a.x * b, a.y * b };
}


inline
double3 operator * (const double3& a, double b) {
    return { a.x * b, a.y * b, a.z * b };
}


inline
double3 operator * (double b, const double3& a) {
    return { a.x * b, a.y * b, a.z * b };
}

//
// dot
//


inline
double dot(const double2& a, const double2& b) {
    return a.x * b.x + a.y * b.y;
}


inline
double dot(const double3& a, const double3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

#endif /* MU3E_UTIL_DOUBLE_HPP_ */
