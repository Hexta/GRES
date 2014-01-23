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

#include "functions.h"

#include "Coords3D.h"

#include <QTime>
#include <QDebug>

#include <random>
#include <functional>

using std::swap;
using std::vector;

namespace {

const int NUMBER_OF_ATOMS_IN_CELL = 8;

Coords3D atomTypes[] = {
    {1.0, 1.0, 1.0},
    {1.0, 0.5, 0.5},
    {0.5, 1.0, 0.5},
    {0.5, 0.5, 1.0},
    {0.25, 0.25, 0.25},
    {0.75, 0.75, 0.25},
    {0.75, 0.25, 0.75},
    {0.25, 0.75, 0.75}
};

struct length {
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float distance;
};

typedef vector<length>lengthes;

struct rectangle {
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float x3;
    float y3;
    float z3;
    float x4;
    float y4;
    float z4;
    float S;
};

typedef vector<rectangle>rectangles;

double
distance(const double& x1, const double& y1, const double& z1, const double& x2,
const double& y2, const double& z2) {
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2));
}

bool rect_comp(const rectangle &r1, const rectangle &r2) {

    const double S1 = distance(r1.x1, r1.y1, r1.z1, r1.x2, r1.y2, r1.z2) *
        distance(r1.x3, r1.y3, r1.z3, r1.x4, r1.y4, r1.z4);
    const double S2 = distance(r2.x1, r2.y1, r2.z1, r2.x2, r2.y2, r2.z2) *
        distance(r2.x3, r2.y3, r2.z3, r2.x4, r2.y4, r2.z4);
    return S1 < S2;
}

} // namespace

bool cmp_float(double x, double y) {
    const double FLOAT_TOL = 0.000001;
    return fabs(x - y) < FLOAT_TOL;
}

Cell findCell(int h, int k, int l, float &xs, float &ys, float &zs,
    Coords3D &vX, Coords3D &vY, Coords3D &vZ)
{
    const int SIZE_X = 5;
    const int SIZE_Y = 5;
    const int SIZE_Z = 5;
    Atoms allAtoms, atomsP1; // атомы, лежащие на плоскостях
    lengthes ls;
    rectangles rectanglesP1; //прямоугольники
    //создадим кристалл из нескольких ячеек
    for (int z = 0; z < 9; ++z)
        for (int y = 0; y < 9; ++y)
            for (int x = 0; x < 9; ++x)
                for (int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a) {
                    Coords3D atom = {x + atomTypes[a].x, y + atomTypes[a].y,
                                     z + atomTypes[a].z};
                    allAtoms.push_back(atom);
                }
    //Найдем свободный член в уравнении секущей плоскости hx+ky+lz-C=0
    int C = (h * (SIZE_X - 1) + k * (SIZE_Y - 1) + l * (SIZE_Z - 1)) / 2 + 1;
    //Найдем атомы, лежащие на плоскости №1

    atomsP1.reserve(SIZE_Z * SIZE_Y * SIZE_X * NUMBER_OF_ATOMS_IN_CELL);
    for (int z = 0; z < SIZE_Z; ++z)
        for (int y = 0; y < SIZE_Y; ++y)
            for (int x = 0; x < SIZE_X; ++x)
                for (int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a)
                    if (cmp_float(h * (x + atomTypes[a].x) + k * (y + atomTypes[a].y) + l * (z + atomTypes[a].z), C)) {
                        Coords3D atom = {x + atomTypes[a].x, y + atomTypes[a].y, z + atomTypes[a].z};
                        atomsP1.push_back(atom);
                    }

    //Найдем прямоугольники, лежащие на плоскости №1
    const size_t atomsP1Size = atomsP1.size();
    ls.reserve(atomsP1Size * atomsP1Size / 2);
    for (unsigned int i = 0; i < atomsP1Size; ++i)
        for (unsigned int j = i + 1; j < atomsP1Size; ++j) {
            const auto atomP1A = atomsP1[i];
            const auto atomP1B = atomsP1[j];
            float l1 = distance(atomP1A.x, atomP1A.y, atomP1A.z,
                                atomP1B.x, atomP1B.y, atomP1B.z);
            length L = {atomP1A.x, atomP1A.y, atomP1A.z,
                        atomP1B.x, atomP1B.y, atomP1B.z, l1};
            ls.push_back(L);
        }
    float kx = -1;
    float ky = -10;
    float kz = -100;

    for (unsigned int i = 0; i < ls.size(); ++i)
        for (unsigned int j = i + 1; j < ls.size(); ++j) {
            float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, l1, l2;
            x1 = ls[i].x1;
            y1 = ls[i].y1;
            z1 = ls[i].z1;
            x2 = ls[i].x2;
            y2 = ls[i].y2;
            z2 = ls[i].z2;
            x3 = ls[j].x1;
            y3 = ls[j].y1;
            z3 = ls[j].z1;
            x4 = ls[j].x2;
            y4 = ls[j].y2;
            z4 = ls[j].z2;
            l1 = ls[i].distance;
            l2 = ls[j].distance;

            if (cmp_float(l1, l2)) {
                if (!cmp_float(x4, x3))
                    kx = (x2 - x1) / (x4 - x3);
                else if (cmp_float(x2, x1))
                    kx = 1;

                if (!cmp_float(y4, y3))
                    ky = (y2 - y1) / (y4 - y3);
                else if (cmp_float(y2, y1))
                    ky = 1;

                if (!cmp_float(z4, z3))
                    kz = (z2 - z1) / (z4 - z3);
                else if (cmp_float(z2, z1))
                    kz = 1;

                if ((cmp_float(kx, ky) && cmp_float(kx, kz)) &&
                    cmp_float((x2 - x1)*(x4 - x2) + (y2 - y1)*(y4 - y2) + (z2 - z1)*(z4 - z2), 0)) {
                    rectangle rect = {x1, y1, z1, x2, y2, z2,
                                      x3, y3, z3, x4, y4, z4, l1};
                    rectanglesP1.push_back(rect);

                    //удаляем дубли
                    for (auto v = ls.begin() + j + 1; v < ls.end(); ++v) {
                        const auto& v_x1 = v->x1;
                        const auto& v_x2 = v->x2;
                        const auto& v_y1 = v->y1;
                        const auto& v_y2 = v->y2;
                        const auto& v_z1 = v->z1;
                        const auto& v_z2 = v->z2;

                        if (((cmp_float(v_x1, x1) && cmp_float(v_x2, x3)
                              && cmp_float(v_y1, y1) && cmp_float(v_y2, y3)
                              && cmp_float(v_z1, z1) && cmp_float(v_z2, z3))
                             || (cmp_float(v_x1, x2) && cmp_float(v_x2, x4)
                                 && cmp_float(v_y1, y2) && cmp_float(v_y2, y4)
                                 && cmp_float(v_z1, z2) && cmp_float(v_z2, z4)))) {
                            const size_t n = v - ls.begin();
                            ls.erase(v);
                            v = ls.begin() + n - 1;
                        }
                    }
                }
            }
        }
    stable_sort(rectanglesP1.begin(), rectanglesP1.end(), rect_comp);

    Cells allCells;

    for (auto &rectangle : rectanglesP1) {
        float x1 = rectangle.x1, y1 = rectangle.y1, z1 = rectangle.z1;
        float x2 = rectangle.x2, y2 = rectangle.y2, z2 = rectangle.z2;
        float x3 = rectangle.x3, y3 = rectangle.y3, z3 = rectangle.z3;
        float x4 = rectangle.x4, y4 = rectangle.y4, z4 = rectangle.z4;

        for (float n = 0.5; n < 5; n += 0.5) {
            int atoms = 0;

            for (auto &atom : allAtoms) {
                if ((cmp_float(atom.x, x1 + n * h)
                     && cmp_float(atom.y, y1 + n * k)
                     && cmp_float(atom.z, z1 + n * l))
                    || (cmp_float(atom.x, x2 + n * h)
                        && cmp_float(atom.y, y2 + n * k)
                        && cmp_float(atom.z, z2 + n * l))
                    || (cmp_float(atom.x, x3 + n * h) && cmp_float(atom.y, y3 + n * k)
                        && cmp_float(atom.z, z3 + n * l))
                    || (cmp_float(atom.x, x4 + n * h)
                        && cmp_float(atom.y, y4 + n * k)
                        && cmp_float(atom.z, z4 + n * l)))
                    ++atoms;
            }

            if (atoms == 4) {
                Coords3D Vx, Vy, Vz;

                Coords3D const P1(x3, y3, z3);
                Coords3D const P2(x4, y4, z4);
                Coords3D const P3(x1, y1, z1);
                Coords3D const P4(x3 + n*h, y3 + n*k, z3 + n * l);

                Vz = P4 - P1;
                Vy = P3 - P1;
                Vx = P2 - P1;
                Atoms cellAtoms = Cell(allAtoms, Vx, Vy, Vz, P1).atoms;

                // а в конец списка атомов запишем векторы координат и координаты начала координат :)
                cellAtoms.push_back(P1);
                cellAtoms.push_back(Vx);
                cellAtoms.push_back(Vy);
                cellAtoms.push_back(Vz);
                allCells.push_back(cellAtoms);
            }
        }
    }
    stable_sort(allCells.begin(), allCells.end());

    for (auto &cell : allCells) {
        const size_t cell_size = cell.size();

        //считаем векторы координат
        const Coords3D &P1 = *(cell.atoms.end() - 4);
        const Coords3D &Vx = *(cell.atoms.end() - 3);
        const Coords3D &Vy = *(cell.atoms.end() - 2);
        const Coords3D &Vz = *(cell.atoms.end() - 1);

        if (cell + Vx == allAtoms //транслируем по OX
            && cell + Vy == allAtoms //транслируем по OY
            && cell + Vz == allAtoms //транслируем по OZ
            && cell + (-1) * Vz == allAtoms //транслируем по -OZ
            && cell + (-1) * Vy == allAtoms //транслируем по -OY
            && cell + (-1) * Vx == allAtoms //транслируем по -OX
            && Cell(allAtoms, Vx, Vy, Vz, P1 + Vx).size() == cell_size - 4
            && Cell(allAtoms, Vx, Vy, Vz, P1 + Vy).size() == cell_size - 4
            && Cell(allAtoms, Vx, Vy, Vz, P1 + Vz).size() == cell_size - 4
            && Cell(allAtoms, Vx, Vy, Vz, P1 + -1 * Vx).size() == cell_size - 4
            && Cell(allAtoms, Vx, Vy, Vz, P1 + -1 * Vy).size() == cell_size - 4
            && Cell(allAtoms, Vx, Vy, Vz, P1 + -1 * Vz).size() == cell_size - 4) {

            auto tmpCell = cell;
            allCells.clear();
            allCells.push_back(tmpCell);
            vX = Vx;
            vY = Vy;
            vZ = Vz;
            break;
        }
    }

    auto firstCell = allCells[0];

    firstCell.moveCoords(*(firstCell.atoms.end() - 4),
        *(firstCell.atoms.end() - 3),
        *(firstCell.atoms.end() - 2),
        *(firstCell.atoms.end() - 1));

    xs = firstCell.getXSize();
    ys = firstCell.getYSize();
    zs = firstCell.getZSize();

    firstCell.optimize();

    return firstCell;
}
