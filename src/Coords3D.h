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

#pragma once

struct Coords3D {
    float x;
    float y;
    float z;

    Coords3D();
    Coords3D(float x, float y, float z);

    const Coords3D operator +(const Coords3D& v) const;

    //Разность координат двух точек
    const Coords3D operator -(const Coords3D& v) const;

    //Скалярное произведение векторов
    double operator *(const Coords3D& v) const;

    bool operator ==(const Coords3D& v) const;

    //Возведение вектора в квадрат
    double sqr() const;
    double length() const; //Вычисляет длину вектора

private:
    static float xPrev;
    static float yPrev;
    static float zPrev;

    static double resultPrev;
}; 

Coords3D operator *(const Coords3D& v, int n);

Coords3D operator *(int n, const Coords3D& v);
