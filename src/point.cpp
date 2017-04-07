/*
 *   Point.cpp
 *   represents a point in 2D
 *   Author: ROBOTIS
 *
 */

#include <math.h>
#include "point.h"

using namespace darwin;


point2d::point2d()
        : x(0), y(0) {}


point2d::point2d(double x, double y)
        : x(x), y(y) {}


point2d::~point2d() {
}


/*returns the euclidian distance between pt1 and pt2*/
double point2d::distance(const point2d& pt1, const point2d& pt2) {
    double x = pt1.x - pt2.x;
    double y = pt1.y - pt2.y;
    return sqrt(x * x + y * y);
}


point2d& point2d::operator+=(const point2d& point) {
    x += point.x;
    y += point.y;
    return *this;
}


point2d& point2d::operator-=(const point2d& point) {
    x -= point.x;
    y -= point.y;
    return *this;
}


point2d& point2d::operator+=(double value) {
    x += value;
    y += value;
    return *this;
}


point2d& point2d::operator-=(double value) {
    x -= value;
    y -= value;
    return *this;
}


point2d& point2d::operator*=(double value) {
    x *= value;
    y *= value;
    return *this;
}


point2d& point2d::operator/=(double value) {
    x /= value;
    y /= value;
    return *this;
}


point2d point2d::operator+(const point2d& point) const {
    return point2d(x + point.x, y + point.y);
}


point2d point2d::operator-(const point2d& point) const {
    return point2d(x - point.x, y - point.y);
}


point2d point2d::operator+(double value) const {
    return point2d(x + value, y + value);
}


point2d point2d::operator-(double value) const {
    return point2d(x - value, y - value);
}


point2d point2d::operator*(double value) const {
    return point2d(x * value, y * value);
}


point2d point2d::operator/(double value) const {
    return point2d(x / value, y / value);
}


point3d::point3d()
        : x(0), y(0), z(0) {}


point3d::point3d(double x, double y, double z)
        : x(x), y(y), z(z) {}


point3d::~point3d() {
}


double point3d::distance(const point3d& pt1, const point3d& pt2) {
    double x = pt1.x - pt2.x;
    double y = pt1.y - pt2.y;
    double z = pt1.z - pt2.z;
    return sqrt(x * x + y * y + z * z);
}


point3d& point3d::operator+=(const point3d& point) {
    x += point.x;
    y += point.y;
    z += point.z;
    return *this;
}


point3d& point3d::operator-=(const point3d& point) {
    x -= point.x;
    y -= point.y;
    z -= point.z;
    return *this;
}


point3d& point3d::operator+=(double value) {
    x += value;
    y += value;
    z += value;
    return *this;
}


point3d& point3d::operator-=(double value) {
    x -= value;
    y -= value;
    z -= value;
    return *this;
}


point3d& point3d::operator*=(double value) {
    x *= value;
    y *= value;
    z *= value;
    return *this;
}


point3d& point3d::operator/=(double value) {
    x /= value;
    y /= value;
    z /= value;
    return *this;
}


point3d point3d::operator+(const point3d& point) const {
    return point3d(x + point.x, y + point.y, z + point.z);
}


point3d point3d::operator-(const point3d& point) const {
    return point3d(x - point.x, y - point.y, z - point.z);
}


point3d point3d::operator+(double value) const {
    return point3d(x + value, y + value, z + value);
}


point3d point3d::operator-(double value) const {
    return point3d(x - value, y - value, z - value);
}


point3d point3d::operator*(double value) const {
    return point3d(x * value, y * value, z * value);
}


point3d point3d::operator/(double value) const {
    return point3d(x / value, y / value, z / value);
}
