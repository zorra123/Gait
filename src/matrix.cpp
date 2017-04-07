/*
 *   Matrix.cpp
 *   Matrix3D is the class of matrix of size 4X4.
 *   Author: ROBOTIS
 *
 */

#include <math.h>
#include "matrix.h"

using namespace darwin;


/**
create the identity matrix (of size 4X4)
*/
matrix::matrix() {
    this->identity();
}


/**
TODO: no clear about the side effects
*/
matrix::matrix(const matrix& mat) {
    *this = mat;
}


matrix::~matrix() {
}


/*
set the current matrix to identity
*/
void matrix::identity() {
    m[m00] = 1;
    m[m01] = 0;
    m[m02] = 0;
    m[m03] = 0;
    m[m10] = 0;
    m[m11] = 1;
    m[m12] = 0;
    m[m13] = 0;
    m[m20] = 0;
    m[m21] = 0;
    m[m22] = 1;
    m[m23] = 0;
    m[m30] = 0;
    m[m31] = 0;
    m[m32] = 0;
    m[m33] = 1;
}


/*
compute the inverse of the matrix
returns true iff the matrix is invertible
*/
bool matrix::inverse() {
    matrix src, dst, tmp;
    double det;

    /* transpose matrix */
    for (int i = 0; i < 4; i++) {
        src.m[i] = m[i * 4];
        src.m[i + 4] = m[i * 4 + 1];
        src.m[i + 8] = m[i * 4 + 2];
        src.m[i + 12] = m[i * 4 + 3];
    }

    /* calculate pairs for first 8 elements (cofactors) */
    tmp.m[0] = src.m[10] * src.m[15];
    tmp.m[1] = src.m[11] * src.m[14];
    tmp.m[2] = src.m[9] * src.m[15];
    tmp.m[3] = src.m[11] * src.m[13];
    tmp.m[4] = src.m[9] * src.m[14];
    tmp.m[5] = src.m[10] * src.m[13];
    tmp.m[6] = src.m[8] * src.m[15];
    tmp.m[7] = src.m[11] * src.m[12];
    tmp.m[8] = src.m[8] * src.m[14];
    tmp.m[9] = src.m[10] * src.m[12];
    tmp.m[10] = src.m[8] * src.m[13];
    tmp.m[11] = src.m[9] * src.m[12];
    /* calculate first 8 elements (cofactors) */
    dst.m[0] = (tmp.m[0] * src.m[5] + tmp.m[3] * src.m[6] + tmp.m[4] * src.m[7]) -
               (tmp.m[1] * src.m[5] + tmp.m[2] * src.m[6] + tmp.m[5] * src.m[7]);
    dst.m[1] = (tmp.m[1] * src.m[4] + tmp.m[6] * src.m[6] + tmp.m[9] * src.m[7]) -
               (tmp.m[0] * src.m[4] + tmp.m[7] * src.m[6] + tmp.m[8] * src.m[7]);
    dst.m[2] = (tmp.m[2] * src.m[4] + tmp.m[7] * src.m[5] + tmp.m[10] * src.m[7]) -
               (tmp.m[3] * src.m[4] + tmp.m[6] * src.m[5] + tmp.m[11] * src.m[7]);
    dst.m[3] = (tmp.m[5] * src.m[4] + tmp.m[8] * src.m[5] + tmp.m[11] * src.m[6]) -
               (tmp.m[4] * src.m[4] + tmp.m[9] * src.m[5] + tmp.m[10] * src.m[6]);
    dst.m[4] = (tmp.m[1] * src.m[1] + tmp.m[2] * src.m[2] + tmp.m[5] * src.m[3]) -
               (tmp.m[0] * src.m[1] + tmp.m[3] * src.m[2] + tmp.m[4] * src.m[3]);
    dst.m[5] = (tmp.m[0] * src.m[0] + tmp.m[7] * src.m[2] + tmp.m[8] * src.m[3]) -
               (tmp.m[1] * src.m[0] + tmp.m[6] * src.m[2] + tmp.m[9] * src.m[3]);
    dst.m[6] = (tmp.m[3] * src.m[0] + tmp.m[6] * src.m[1] + tmp.m[11] * src.m[3]) -
               (tmp.m[2] * src.m[0] + tmp.m[7] * src.m[1] + tmp.m[10] * src.m[3]);
    dst.m[7] = (tmp.m[4] * src.m[0] + tmp.m[9] * src.m[1] + tmp.m[10] * src.m[2]) -
               (tmp.m[5] * src.m[0] + tmp.m[8] * src.m[1] + tmp.m[11] * src.m[2]);
    /* calculate pairs for second 8 elements (cofactors) */
    tmp.m[0] = src.m[2] * src.m[7];
    tmp.m[1] = src.m[3] * src.m[6];
    tmp.m[2] = src.m[1] * src.m[7];
    tmp.m[3] = src.m[3] * src.m[5];
    tmp.m[4] = src.m[1] * src.m[6];
    tmp.m[5] = src.m[2] * src.m[5];
    //Streaming SIMD Extensions - inverse of 4x4 matrix 8
    tmp.m[6] = src.m[0] * src.m[7];
    tmp.m[7] = src.m[3] * src.m[4];
    tmp.m[8] = src.m[0] * src.m[6];
    tmp.m[9] = src.m[2] * src.m[4];
    tmp.m[10] = src.m[0] * src.m[5];
    tmp.m[11] = src.m[1] * src.m[4];
    /* calculate second 8 elements (cofactors) */
    dst.m[8] = (tmp.m[0] * src.m[13] + tmp.m[3] * src.m[14] + tmp.m[4] * src.m[15]) -
               (tmp.m[1] * src.m[13] + tmp.m[2] * src.m[14] + tmp.m[5] * src.m[15]);
    dst.m[9] = (tmp.m[1] * src.m[12] + tmp.m[6] * src.m[14] + tmp.m[9] * src.m[15]) -
               (tmp.m[0] * src.m[12] + tmp.m[7] * src.m[14] + tmp.m[8] * src.m[15]);
    dst.m[10] = (tmp.m[2] * src.m[12] + tmp.m[7] * src.m[13] + tmp.m[10] * src.m[15]) -
                (tmp.m[3] * src.m[12] + tmp.m[6] * src.m[13] + tmp.m[11] * src.m[15]);
    dst.m[11] = (tmp.m[5] * src.m[12] + tmp.m[8] * src.m[13] + tmp.m[11] * src.m[14]) -
                (tmp.m[4] * src.m[12] + tmp.m[9] * src.m[13] + tmp.m[10] * src.m[14]);
    dst.m[12] = (tmp.m[2] * src.m[10] + tmp.m[5] * src.m[11] + tmp.m[1] * src.m[9]) -
                (tmp.m[4] * src.m[11] + tmp.m[0] * src.m[9] + tmp.m[3] * src.m[10]);
    dst.m[13] = (tmp.m[8] * src.m[11] + tmp.m[0] * src.m[8] + tmp.m[7] * src.m[10]) -
                (tmp.m[6] * src.m[10] + tmp.m[9] * src.m[11] + tmp.m[1] * src.m[8]);
    dst.m[14] = (tmp.m[6] * src.m[9] + tmp.m[11] * src.m[11] + tmp.m[3] * src.m[8]) -
                (tmp.m[10] * src.m[11] + tmp.m[2] * src.m[8] + tmp.m[7] * src.m[9]);
    dst.m[15] = (tmp.m[10] * src.m[10] + tmp.m[4] * src.m[8] + tmp.m[9] * src.m[9]) -
                (tmp.m[8] * src.m[9] + tmp.m[11] * src.m[10] + tmp.m[5] * src.m[8]);
    /* calculate determinant */
    det = src.m[0] * dst.m[0] + src.m[1] * dst.m[1] + src.m[2] * dst.m[2] + src.m[3] * dst.m[3];
    /* calculate matrix inverse */
    if (det == 0) {
        det = 0;
        return false;
    } else {
        det = 1 / det;
    }

    for (int i = 0; i < MAXNUM_ELEMENT; i++)
        m[i] = dst.m[i] * det;

    return true;
}


/*
compute the matrix
scale.x         0        0        0
0            scale.y     0        0
0               0      scale.z    0
0               0        0        1
*/
void matrix::scale(const vector& scale) {
    matrix mat;
    mat.m[m00] = scale.x;
    mat.m[m11] = scale.y;
    mat.m[m22] = scale.z;

    *this *= mat;
}


/*construct the matrix corresponding to a rotation of angle around the axis*/
void matrix::rotate(double angle, const vector& axis) {
    double rad = angle * 3.141592 / 180.0;
    double C = cos(rad);
    double S = sin(rad);
    matrix mat;

    mat.m[m00] = C + axis.x * axis.x * (1 - C);
    mat.m[m01] = axis.x * axis.y * (1 - C) - axis.z * S;
    mat.m[m02] = axis.x * axis.z * (1 - C) + axis.y * S;
    mat.m[m10] = axis.x * axis.y * (1 - C) + axis.z * S;
    mat.m[m11] = C + axis.y * axis.y * (1 - C);
    mat.m[m12] = axis.y * axis.z * (1 - C) - axis.x * S;
    mat.m[m20] = axis.x * axis.z * (1 - C) - axis.y * S;
    mat.m[m21] = axis.y * axis.z * (1 - C) + axis.x * S;
    mat.m[m22] = C + axis.z * axis.z * (1 - C);

    *this *= mat;
}


/*compute the matrix corresponding to a translation of vector offset*/
void matrix::translate(const vector& offset) {
    matrix mat;
    mat.m[m03] = offset.x;
    mat.m[m13] = offset.y;
    mat.m[m23] = offset.z;
    *this *= mat;
}


/*make the product of the current matrix with the point*/
point3d matrix::transform(const point3d& point) const {
    point3d result;
    result.x = m[m00] * point.x + m[m01] * point.y + m[m02] * point.z + m[m03];
    result.y = m[m10] * point.x + m[m11] * point.y + m[m12] * point.z + m[m13];
    result.z = m[m20] * point.x + m[m21] * point.y + m[m22] * point.z + m[m23];
    return result;
}


/*make the product of the current matrix with the point*/
vector matrix::transform(const vector& vec) const {
    vector result;
    result.x = m[m00] * vec.x + m[m01] * vec.y + m[m02] * vec.z + m[m03];
    result.y = m[m10] * vec.x + m[m11] * vec.y + m[m12] * vec.z + m[m13];
    result.z = m[m20] * vec.x + m[m21] * vec.y + m[m22] * vec.z + m[m23];
    return result;
}


/*?*/
void matrix::set_transform(const point3d& point, const vector& angle) {
    double Cx = cos(angle.x * 3.141592 / 180.0);
    double Cy = cos(angle.y * 3.141592 / 180.0);
    double Cz = cos(angle.z * 3.141592 / 180.0);
    double Sx = sin(angle.x * 3.141592 / 180.0);
    double Sy = sin(angle.y * 3.141592 / 180.0);
    double Sz = sin(angle.z * 3.141592 / 180.0);

    identity();
    m[0] = Cz * Cy;
    m[1] = Cz * Sy * Sx - Sz * Cx;
    m[2] = Cz * Sy * Cx + Sz * Sx;
    m[3] = point.x;
    m[4] = Sz * Cy;
    m[5] = Sz * Sy * Sx + Cz * Cx;
    m[6] = Sz * Sy * Cx - Cz * Sx;
    m[7] = point.y;
    m[8] = -Sy;
    m[9] = Cy * Sx;
    m[10] = Cy * Cx;
    m[11] = point.z;
}


matrix& matrix::operator=(const matrix& mat) {
    for (int i = 0; i < MAXNUM_ELEMENT; i++)
        m[i] = mat.m[i];
    return *this;
}


matrix& matrix::operator*=(const matrix& mat) {
    matrix tmp = *this * mat;
    *this = tmp;
    return *this;
}


matrix matrix::operator*(const matrix& mat) {
    matrix result;

    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
            for (int c = 0; c < 4; c++) {
                result.m[j * 4 + i] += m[j * 4 + c] * mat.m[c * 4 + i];
            }
        }
    }

    return result;
}
