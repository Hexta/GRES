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

#include "Atoms.h"
#include "Coords3D.h"
#include "Cell.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

using std::vector;

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
typedef vector<AtomInfo> cell;
typedef vector<cell> surface1D;
typedef vector<surface1D> surface2D;
typedef vector<surface2D> surface3D;

bool rect_comp(const rectangle &r1, const rectangle &r2);

double distance(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2);

bool compareTranslCell(const Cell &cellAtoms, const Cell &allAtoms, const Coords3D &V); //Совпадение атомов после трансляции
void coordsMove(Cell &ca, const Coords3D &O, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz); // смещение координат в ячейку

void recallNeighbours(surface3D &surface, vector<AtomType> &surfAtoms, int x, int y, int z, int type);

bool selAtom(surface3D &surface, vector<AtomType> &surfAtoms, AllNeighbors &neighbs, int z_min, Cell &tA, const vector<bool> &mask, const float *rates);
bool selAtomCA(surface3D &surface, vector<AtomType> &surfAtoms, int z_min, Cell &tA, vector<bool> &mask, float* rates);

Cell findCell(int h, int k, int l, float &xs, float &ys, float &zs, Coords3D &Vx, Coords3D &Vy, Coords3D &Vz); //Поиск эл ячейки
void addLayer(surface3D &surface, const AllNeighbors& sosedi, int sX, int sY, int sZ);

Cell atomsInBox(const Atoms &atoms, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz, const Coords3D &P1);
void findZmin(const surface3D &surface, int &zm);
void optimizeSurface(surface3D &surface, int z_min);
void delAtom(surface3D &surface, vector<AtomType> &surfAtoms, int x, int y, int z, int type, int surfAtN);

bool cmp_float(double x, double y);
