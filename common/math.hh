/*
 * math.hpp
 *
 *  Created on: Feb 16, 2018
 *      Author: akozlins
 */

#ifndef MU3E_UTIL_MATH_HPP_
#define MU3E_UTIL_MATH_HPP_

//#include "../log.hpp"

#include <cassert>
#include <cfloat>
#include <cmath> // std::abs(float)
#include <cstdlib> // std::abs(int)
#include <cstdint> // int32_t
#include <cstring> // memcpy

namespace mu3e { namespace util {

template < class T >
const T& clamp(const T& v, const T& lo, const T& hi) {
    assert(!(hi < lo));
    if(v < lo) return lo;
    if(hi < v) return hi;
    return v;
}

template < class T >
T pown(T x, int n) {
    assert(n >= 0);
    T r = 1;
    while(true) {
        if(n & 1) r *= x;
        n >>= 1;
        if(n == 0) break;
        x *= x;
    }
    return r;
}

/**
 * Calculate '1 - x / tan(x)' using Pade approximant [1].
 *
 * Pade [4,4] = 1 - x / tan(x) = x^2 * (315 - 14 * x^2) / (945 - (105 - x^2) * x^2)
 * Error: 0 < |x| < 1 => |ulp| <= 7
 *
 * [1] 'https://en.wikipedia.org/wiki/Pade_approximant'
 */
template < typename float_t >
float_t xtanx(float_t x) {
    assert(std::abs(x) < float_t(M_PI_2));

    // Taylor: 1 - x / tan(x) = x^2 / 3 + x^4 / 45 + x^6 * 2 / 945 + O(X^8)
    // Pade [4,2]: x2 * (35 - x2) / (105 - 10 * x2)
    // Pade [4,4]: x2 * (315 - 14 * x2) / (945 - (105 - x2) * x2)
    // Pade [6,4]: x2 * (165 - (189 - x2) * x2 / 21) / (495 - (60 - x2) * x2)

    float_t x2 = x * x;
    return x2 * (315 - 14 * x2) / (945 - (105 - x2) * x2);
}

template < typename float_t >
float_t dphi(float_t phi1, float_t phi2) {
    float_t dphi = phi2 - phi1;
//    while(std::abs(dphi) > float_t(M_PI)) dphi -= std::copysign(float_t(2*M_PI), dphi);
    if(std::abs(dphi) > float_t(M_PI)) dphi -= std::copysign(float_t(2*M_PI), dphi);
    return dphi;
}

/**
 * Returns angle between sides A and B of the triangle ABC.
 */
template < typename float_t >
float_t angleAB(float_t A, float_t B, float_t C) {
//    float_t dot = (A*A + B*B - C*C) / (2 * A * B);
    float_t dot = (A + B + C) * (A + B - C) / (2 * A * B) - 1;
    if(dot > +1) return 0;
    if(dot < -1) return float_t(M_PI);
    return std::acos(dot);
}

} } // namespace mu3e::util

#endif /* MU3E_UTIL_MATH_HPP_ */
