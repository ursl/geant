/*
 * float.hpp
 *
 *  Created on: Nov 11, 2016
 *      Author: akozlins
 */

#ifndef MU3E_UTIL_FLOAT_HPP_
#define MU3E_UTIL_FLOAT_HPP_

#include "math.hh"

//
// float2
//

struct float2 {
    float x, y;

    float rt() const { return std::hypot(x, y); }
    float phi() const { return std::atan2(y, x); }
};

inline
float float2_rt(const float2& f) {
    return std::hypot(f.x, f.y);
}

//
// float3
//

struct float3 {
    float x, y, z;

    float r() const { return std::sqrt(x*x + y*y + z*z); }
    float rt() const { return std::hypot(x, y); }
    float phi() const { return std::atan2(y, x); }
    float lam() const { return std::atan(z / rt()); }
};

inline
float2 make_float2(float x, float y) {
    return { x, y };
}

inline
float2 make_float2(double x, double y) {
    return { (float)x, (float)y };
}

inline
float2 make_float2(const float3& f3) {
    return { f3.x, f3.y };
}

inline
float3 make_float3(float x, float y, float z) {
    return { x, y, z };
}

inline
float3 make_float3(double x, double y, double z) {
    return { (float)x, (float)y, (float)z };
}

inline
float float3_phi(const float3& f) {
    return std::atan2(f.y, f.x);
}

//
// float4
//

struct float4 {
    float x, y, z, w;
};

//
// add
//

inline
float2 operator + (const float2& a, const float2& b) {
    return { a.x + b.x, a.y + b.y };
}

inline
float2& operator += (float2& a, const float2& b) {
    a.x += b.x; a.y += b.y;
    return a;
}

inline
float3 operator + (const float3& a, const float3& b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

inline
float3& operator += (float3& a, const float3& b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

inline
float4 operator + (const float4& a, const float4& b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

inline
float4& operator += (float4& a, const float4& b) {
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
    return a;
}

//
// sub
//

inline
float2 operator - (const float2& a, const float2& b) {
    return { a.x - b.x, a.y - b.y };
}

inline
float2& operator -= (float2& a, const float2& b) {
    a.x -= b.x; a.y -= b.y;
    return a;
}

inline
float3 operator - (const float3& a, const float3& b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

inline
float3& operator -= (float3& a, const float3& b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

inline
float4 operator - (const float4& a, const float4& b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
}

inline
float4& operator -= (float4& a, const float4& b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z; a.w -= b.w;
    return a;
}

//
// mul
//

inline
float2 operator * (const float2& a, float b) {
    return { a.x * b, a.y * b };
}

inline
float2 operator * (float b, const float2& a) {
    return { a.x * b, a.y * b };
}

inline
float3 operator * (const float3& a, float b) {
    return { a.x * b, a.y * b, a.z * b };
}

inline
float3 operator * (float b, const float3& a) {
    return { a.x * b, a.y * b, a.z * b };
}

inline
float4 operator * (const float4& a, float b) {
    return { a.x * b, a.y * b, a.z * b, a.w * b };
}

inline
float4 operator * (float b, const float4& a) {
    return { a.x * b, a.y * b, a.z * b, a.w * b };
}

//
// div
//

inline
float2 operator / (const float2& a, float b) {
    return { a.x / b, a.y / b };
}

inline
float3 operator / (const float3& a, float b) {
    return { a.x / b, a.y / b, a.z / b };
}

inline
float4 operator / (const float4& a, float b) {
    return { a.x / b, a.y / b, a.z / b, a.w / b };
}

//
// dot
//

inline
float dot(const float2& a, const float2& b) {
    return a.x * b.x + a.y * b.y;
}

inline
float dot(const float3& a, const float3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline
float dot(const float4& a, const float4& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

//
// cross
//

inline
float3 cross(const float3& a, const float3& b) {
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}

#endif /* MU3E_UTIL_FLOAT_HPP_ */
