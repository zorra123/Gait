/*
 *   Vector.h
 *
 *   Author: ROBOTIS
 *
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "point.h"

namespace darwin {
    class vector {
    private:

    protected:

    public:
        double x;
        double y;
        double z;

        vector();

        vector(double x, double y, double z);

        vector(const point3d& pt1, const point3d& pt2);

        vector(const vector& vector);

        ~vector();

        double length() const;

        void normalize();

        double dot(const vector& vec) const;

        vector cross(const vector& vec) const;

        double angle_between(const vector& vec) const;

        double angle_between(const vector& vec, const vector& axis) const;

        vector& operator=(const vector& vector);

        vector& operator+=(const vector& vector);

        vector& operator-=(const vector& vector);

        vector& operator+=(const double value);

        vector& operator-=(const double value);

        vector& operator*=(const double value);

        vector& operator/=(const double value);

        vector operator+(const vector& vec);

        vector operator-(const vector& vec);

        vector operator+(const double value);

        vector operator-(const double value);

        vector operator*(const double value);

        vector operator/(const double value);
    };
}

#endif