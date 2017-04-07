/*
 *   Point.h
 *
 *   Author: ROBOTIS
 *
 */

#ifndef _POINT_H_
#define _POINT_H_


namespace darwin {
    class point2d {
    public:
        double x;
        double y;

        point2d();

        point2d(double x, double y);

        point2d(const point2d& point) = default;

        ~point2d();

        /*compute the euclidean distance between pt1 and pt2*/
        static double distance(const point2d& pt1, const point2d& pt2);

        point2d& operator=(const point2d& point) = default;

        point2d& operator+=(const point2d& point);

        point2d& operator-=(const point2d& point);

        point2d& operator+=(double value);

        point2d& operator-=(double value);

        point2d& operator*=(double value);

        point2d& operator/=(double value);

        point2d operator+(const point2d& point) const;

        point2d operator-(const point2d& point) const;

        point2d operator+(double value) const;

        point2d operator-(double value) const;

        point2d operator*(double value) const;

        point2d operator/(double value) const;
    };

    class point3d {
    public:
        double x;
        double y;
        double z;

        point3d();

        point3d(double x, double y, double z);

        point3d(const point3d& point) = default;

        ~point3d();

        /*compute the euclidean distance between pt1 and pt2*/
        static double distance(const point3d& pt1, const point3d& pt2);

        point3d& operator=(const point3d& point) = default;

        point3d& operator+=(const point3d& point);

        point3d& operator-=(const point3d& point);

        point3d& operator+=(double value);

        point3d& operator-=(double value);

        point3d& operator*=(double value);

        point3d& operator/=(double value);

        point3d operator+(const point3d& point) const;

        point3d operator-(const point3d& point) const;

        point3d operator+(double value) const;

        point3d operator-(double value) const;

        point3d operator*(double value) const;

        point3d operator/(double value) const;
    };
}

#endif
