/*
 *   Matrix.h
 *
 *   Author: ROBOTIS
 *
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "vector.h"
#include "point.h"


namespace darwin {
    class matrix {
    public:
        enum {
            m00 = 0,
            m01,
            m02,
            m03,
            m10,
            m11,
            m12,
            m13,
            m20,
            m21,
            m22,
            m23,
            m30,
            m31,
            m32,
            m33,
            MAXNUM_ELEMENT
        };

        double m[MAXNUM_ELEMENT]; // Element

        matrix();

        matrix(const matrix& mat);

        ~matrix();

        void identity();

        bool inverse();

        void scale(const vector& scale);

        void rotate(double angle, const vector& axis);

        void translate(const vector& offset);

        point3d transform(const point3d& point) const;

        vector transform(const vector& vec) const;

        void set_transform(const point3d& point, const vector& angle);

        matrix& operator=(const matrix& mat);

        matrix& operator*=(const matrix& mat);

        matrix operator*(const matrix& mat);
    };
}

#endif