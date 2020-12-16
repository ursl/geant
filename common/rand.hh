/*
 * rand.hpp
 *
 *  Created on: May 31, 2017
 *      Author: akozlins
 */

#ifndef MU3E_UTIL_RAND_HPP_
#define MU3E_UTIL_RAND_HPP_

#include "double.hh"

#include <CLHEP/Random/RandExponential.h>
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandPoisson.h>

namespace mu3e { namespace util {

inline
double2 rand_u2d() {
    double x, y, r;
    do {
        x = 1 - 2 * CLHEP::RandFlat::shoot();
        y = 1 - 2 * CLHEP::RandFlat::shoot();
        r = x*x + y*y;
    } while(r > 1);
    return { x, y };
}

/**
 * Random point on unit sphere.
 *
 * cpu: 2 uniform + sqrt + sincos
 */
//inline
//CLHEP::Hep3Vector rand_u3d() {
//    double phi = M_PI * CLHEP::RandFlat::shoot(-1.0, 1.0);
//    double costh = CLHEP::RandFlat::shoot(-1.0, 1.0);
//    double sinth = std::sqrt(1 - costh * costh);
//    return { sinth * std::cos(phi), sinth * std::sin(phi), costh };
//}

/**
 * Random point on unit sphere.
 *
 * G.Marsaglia, http://projecteuclid.org/euclid.aoms/1177692644
 * cpu: 2.55 uniform + sqrt
 */
inline
double3 rand_u3d() {
    double x, y, s;
    do {
        x = 1 - 2 * CLHEP::RandFlat::shoot();
        y = 1 - 2 * CLHEP::RandFlat::shoot();
        s = x*x + y*y;
    } while(s > 1);
    double a = 2 * std::sqrt(1 - s);
    return { x * a, y * a, 1 - 2 * s };
}

} } // namespace mu3e::util

#endif /* MU3E_UTIL_RAND_HPP_ */
