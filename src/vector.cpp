/*
 *   Vector.cpp
 *
 *   Author: ROBOTIS
 *
 */

#include <math.h>
#include "vector.h"

using namespace darwin;


vector::vector() : x(0), y(0), z(0) {}


vector::vector(double x, double y, double z)
        : x(x), y(y), z(z) {}


vector::vector(const point3d& pt1, const point3d& pt2) {
    x = pt2.x - pt1.x;
    y = pt2.y - pt1.y;
    z = pt2.z - pt1.z;
}


vector::vector(const vector& vector) {
    x = vector.x;
    y = vector.y;
    z = vector.z;
}


vector::~vector() {
}


double vector::length() const {
    return sqrt(x * x + y * y + z * z);
}


void vector::normalize() {
    double length = this->length();
    x = x / length;
    y = y / length;
    z = z / length;
}


double vector::dot(const vector& vec) const {
    return (x * vec.x + y * vec.y + z * vec.z);
}


vector vector::cross(const vector& vec) const {
    vector res;
    res.x = y * vec.z - z * vec.y;
    res.y = z * vec.x - x * vec.z;
    res.z = x * vec.y - y * vec.x;
    return res;
}


double vector::angle_between(const vector& vec) const {
    return acos((x * vec.x + y * vec.y + z * vec.z) / (this->length() * vec.length())) * (180.0 / 3.141592);
}


double vector::angle_between(const vector& vec, const vector& axis) const {
    double angle = this->angle_between(vec);
    vector cross = this->cross(vec);
    if (cross.dot(axis) < 0.0)
        angle *= -1.0;

    return angle;
}


vector& vector::operator=(const vector& vector) {
    x = vector.x;
    y = vector.y;
    z = vector.z;
    return *this;
}


vector& vector::operator+=(const vector& vector) {
    x += vector.x;
    y += vector.y;
    z += vector.z;
    return *this;
}


vector& vector::operator-=(const vector& vector) {
    x -= vector.x;
    y -= vector.y;
    z -= vector.z;
    return *this;
}


vector& vector::operator+=(const double value) {
    x += value;
    y += value;
    z += value;
    return *this;
}


vector& vector::operator-=(const double value) {
    x -= value;
    y -= value;
    z -= value;
    return *this;
}


vector& vector::operator*=(const double value) {
    x *= value;
    y *= value;
    z *= value;
    return *this;
}


vector& vector::operator/=(const double value) {
    x /= value;
    y /= value;
    z /= value;
    return *this;
}


vector vector::operator+(const vector& vec) {
    return vector(x + vec.x, y + vec.y, z + vec.z);
}


vector vector::operator-(const vector& vec) {
    return vector(x - vec.x, y - vec.y, z - vec.z);
}


vector vector::operator+(const double value) {
    return vector(x + value, y + value, z + value);
}


vector vector::operator-(const double value) {
    return vector(x - value, y - value, z - value);
}


vector vector::operator*(const double value) {
    return vector(x * value, y * value, z * value);
}


vector vector::operator/(const double value) {
    return vector(x / value, y / value, z / value);
}
