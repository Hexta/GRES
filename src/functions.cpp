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

#include <QTime>
#include <QDebug>

#include <random>
#include <functional>

using std::swap;
using std::vector;

namespace {
float P1 = 1.0f;
float P2 = 0.030f;
float P3 = 0.00010f;

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

void addLayer(Surface3D &surface, const AllNeighbors& sosedi, int sX, int sY, int sZ)
{
    Surface2D surfaceXY;
    surfaceXY.reserve(sY);

    for (int y = 0; y < sY; ++y) {
        Surface1D surfaceX;
        surfaceX.reserve(sX);

        for (int x = 0; x < sX; ++x) {
            CellInfo cell;
            cell.reserve(sosedi.size());
            for (auto &sosed : sosedi) {
                Neighbors neighbs;
                char numberNeighbs = 0; //Число первых соседей
                for (int nb = 0; nb < 4; ++nb) {
                    if (x + sosed[nb].x >= 0 && y + sosed[nb].y >= 0 && x + sosed[nb].x < sX && y + sosed[nb].y < sY) {
                        ++numberNeighbs;
                        AtomType neighb = {x + sosed[nb].x, y + sosed[nb].y, sZ + sosed[nb].z, sosed[nb].type, false};
                        neighbs.push_back(neighb);
                    }
                }
                AtomInfo atom = {neighbs, numberNeighbs, !numberNeighbs};
                cell.push_back(atom);
            }
            surfaceX.push_back(cell);
        }
        surfaceXY.push_back(surfaceX);
    }

    surface.push_back(surfaceXY);
}

Cell atomsInBox(const Atoms &atoms, const Coords3D &Vx, const Coords3D &Vy,
           const Coords3D &Vz, const Coords3D &P1)
{
    Cell cellAts;

    for (auto &atom : atoms) {
        const float x = atom.x;
        const float y = atom.y;
        const float z = atom.z;
        Coords3D X = {x, y, z};
        Coords3D V = X - P1;

        double k = (V * Vz) / Vz.sqr();
        if (k >= 0.0 && k <= 1.0) {
            k = (V * Vy) / Vy.sqr();
            if (k >= 0.0 && k <= 1.0) {
                k = (V * Vx) / Vx.sqr();
                if (k >= 0.0 && k <= 1.0)
                    cellAts.atoms.push_back(X); //Записываем атомы в ячейке
            }
        }
    }

    return cellAts;
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
                Atoms cellAtoms;
                Coords3D Vx, Vy, Vz, P1;
                P1.x = x3;
                P1.y = y3;
                P1.z = z3;
                Coords3D P2 = {x4, y4, z4}, P3 = {x1, y1, z1};
                Coords3D P4 = {x3 + n*h, y3 + n*k, z3 + n * l};

                Vz = P4 - P1;
                Vy = P3 - P1;
                Vx = P2 - P1;
                cellAtoms = atomsInBox(allAtoms, Vx, Vy, Vz, P1).atoms;

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
        //транслируем по OX
        bool happy = cell + Vx == allAtoms;

        //транслируем по OY
        if (happy)
            happy = cell + Vy == allAtoms;
        //транслируем по OZ
        if (happy)
            happy = cell + Vz == allAtoms;

        //транслируем по -OZ
        if (happy)
            happy = cell + (-1) * Vz == allAtoms;
        //транслируем по -OY
        if (happy)
            happy = cell + (-1) * Vy == allAtoms;
        //транслируем по -OX
        if (happy)
            happy = cell + (-1) * Vx == allAtoms;

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + Vx);
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + Vy);
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + Vz);
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + -1 * Vx);
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + -1 * Vy);
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            Cell translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                P1 + -1 * Vz);
            happy = (cell_size - 4 == translCell.size());
        }
        if (happy) {
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

bool selAtom(Surface3D &surface, vector<AtomType> &surfAtoms, AllNeighbors &neighbs,
        int z_min, Cell &tA, const vector<bool> &mask, const float *rates)
{
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = static_cast<int> (surface[z_min].size());
    const int xC = static_cast<int> (surface[z_min][0].size());
    bool result = false;
    bool maskON = !mask.empty();

    static std::mt19937 rdevice(std::random_device {}
                                ());
    static std::uniform_real_distribution<float> dis;
    static auto rd = bind(dis, rdevice);

    int i = static_cast<int> (std::uniform_int_distribution<size_t> {0, surfAtoms.size() - 1}
                              (rdevice));

    auto &surfAtom = surfAtoms[i];
    int x = surfAtom.x;
    int y = surfAtom.y;
    int z = surfAtom.z;
    unsigned char a = surfAtom.type;
    auto &surfaceZY = surface[z][y];
    if ((maskON && ((z == 0
                     && !mask[ (y - 2)*(surfaceZY.size() - 4) + x - 2 ])
                    || (z + tA.atoms[a].z >= 0.5)))
        || !maskON) {
        float randN = rd();
        int bonds = surfaceZY[x][a].fNbCount;
        if (!bonds) {
            swap(surfAtoms[i], surfAtoms.back());
            surfAtoms.pop_back();
            swap(surfAtoms[i], surfAtoms.back());
            --i;
            return false;
        }

        if ((bonds == 1 && randN < P1)
            || (bonds == 2 && randN < P2)
            || (bonds == 3 && randN < P3)) {
            surface.delAtom(surfAtoms, x, y, z, a, i);
            /*
            удалим и соостветсвенный краевой атом
             */
            if (x == 4 || x == 5) {
                surface.delAtom(surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                surface.delAtom(surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                surface.delAtom(surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                surface.delAtom(surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }
            if (z >= surface.size() - 3) {
                addLayer(surface, neighbs, xC, yC, static_cast<int> (surface.size()));
            }
        }
    }

    return result;
}

bool selAtomCA(Surface3D &surface, vector<AtomType> &surfAtoms, int z_min,
          Cell &tA, vector<bool> &mask, float* rates)
{
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = static_cast<int> (surface[z_min].size());
    const int xC = static_cast<int> (surface[z_min][0].size());
    bool result = false;
    bool maskON = !mask.empty();

    size_t surfAtomsSize = surfAtoms.size();

    static std::mt19937 rdevice(std::random_device {}
                                ());
    static std::uniform_real_distribution<float> dis;
    static auto rd = bind(dis, rdevice);

    for (unsigned int i = 0; i < surfAtomsSize; ++i) {
        const auto &surfAtom = surfAtoms[i];
        int x = surfAtom.x;
        int y = surfAtom.y;
        unsigned int z = surfAtom.z;
        unsigned short a = surfAtom.type;
        const auto &tAz = tA.atoms[a].z;
        auto &surfaceZY = surface[z][y];
        if ((maskON && (((z + tAz < 0.5) && !mask[ (y - 2)*(surfaceZY.size() - 4) + x - 2 ])
                        || (z + tAz >= 0.5)))
            || !maskON) {
            float randN = rd();
            int bonds = surfaceZY[x][a].fNbCount;
            if (!bonds) {
                surfAtoms.erase(surfAtoms.begin() + i--);
                --surfAtomsSize;
                continue;
            }

            if ((bonds == 1 && randN < P1)
                || (bonds == 2 && randN < P2)
                || (bonds == 3 && randN < P3)) {
                surfAtoms[i].toDel = true;
                if (z == surface.size() - 3)
                    result = true;
            }
        }
    }

    int i = 0;
    auto surfAtomsEnd = surfAtoms.end();
    for (auto surfAtomIter = surfAtoms.begin(); surfAtomIter != surfAtomsEnd;
         ++surfAtomIter) {
        if (surfAtomIter->toDel) {
            int x = surfAtomIter->x;
            int y = surfAtomIter->y;
            int z = surfAtomIter->z;
            unsigned char a = surfAtomIter->type;

            surface.delAtom(surfAtoms, x, y, z, a, i);

            if (x == 4 || x == 5) {
                surface.delAtom(surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                surface.delAtom(surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                surface.delAtom(surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                surface.delAtom(surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }

            surfAtomIter = surfAtoms.begin() + --i;
            surfAtomsEnd = surfAtoms.end();
        }
        ++i;
    }
    return result;
}
