/******************************************************************************
 * Copyright (c) 2009-2013 Artur Molchanov <artur.molchanov@gmail.com>        *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 ******************************************************************************/

#include "Coords3D.h"

#include <cmath>

namespace {
    const double FLOAT_TOL = 0.00001;

    bool cmp_float(double x, double y) {
        return fabs(x - y) < FLOAT_TOL;
    }
}

const Coords3D Coords3D::operator +(const Coords3D& v) const {

    Coords3D temp = {x + v.x, y + v.y, z + v.z};
    return temp;
}


double Coords3D::resultPrev = 0;
float Coords3D::zPrev = 0;
float Coords3D::yPrev = 0;
float Coords3D::xPrev = 0;

double Coords3D::sqr() const {
    if (cmp_float(xPrev, x) && cmp_float(yPrev, y) && cmp_float(zPrev, z))
        return resultPrev;

    xPrev = x;
    yPrev = y;
    zPrev = z;
    resultPrev = pow(xPrev, 2) + pow(yPrev, 2) + pow(zPrev, 2);

    return resultPrev;
}

double Coords3D::length() const {
    return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

double Coords3D::operator*(const Coords3D& v) const {
    return x * v.x + y * v.y + z * v.z;
}

const Coords3D Coords3D::operator-(const Coords3D& v) const {
    Coords3D result = {x - v.x, y - v.y, z - v.z};
    return result;
}

bool Coords3D::operator==(const Coords3D& v) const {
    return cmp_float(x, v.x) && cmp_float(y, v.y) && cmp_float(z, v.z);
}

Coords3D::Coords3D(float x, float y, float z) : x(x), y(y), z(z) {

}

Coords3D::Coords3D() {
}

void Coords3D::normalize() {
    const float len = static_cast<float> (sqrt(pow(x, 2) + pow(y, 2) +
        pow(z, 2)));

    x /= len;
    y /= len;
    z /= len;
}

Coords3D operator *(const Coords3D& v, int n) {
    Coords3D temp = {n * v.x, n * v.y, n * v.z};
    return temp;
}

Coords3D operator *(int n, const Coords3D& v) {
    return v * n;
}
