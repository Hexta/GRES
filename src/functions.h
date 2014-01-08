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

#include "Coords3D.h"

#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

using std::vector;

struct atomType {
    int x; //сдвиг ячейки по OX (-1;0;+1)
    int y;
    int z;
    unsigned char type; //тип атома (1 -- 8)
    bool toDel;
};
typedef vector<atomType> soseds; //соседи одного атома
typedef vector<soseds> allSoseds; //соседи всех атомов

struct atomInfo {
    soseds neighbours;
    char fNbCount;
    bool deleted;
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
typedef vector<Coords3D> cellAtom;
typedef vector<Coords3D> atomsCoords;
typedef vector<atomsCoords> cells;
typedef vector<length>lengthes;
typedef vector<rectangle>rectangles;
typedef vector<atomInfo> cell;
typedef vector<cell> surface1D;
typedef vector<surface1D> surface2D;
typedef vector<surface2D> surface3D;

void findSoseds(allSoseds &allSosedi, const atomsCoords &atom_Types, float xs,
	float ys, float zs);
bool rect_comp(const rectangle &r1, const rectangle &r2);
bool cells_comp(const atomsCoords &c1, const atomsCoords &c2);
double distance(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2);

bool compareTranslCell(const atomsCoords &cellAtoms, const atomsCoords &allAtoms, const Coords3D &V); //Совпадение атомов после трансляции
void coordsMove(atomsCoords &ca, const Coords3D &O, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz); // смещение координат в ячейку
void cellOptim(atomsCoords &ca, float&, float&, float&); //чистка от общих атомов

void recallNeighbours(surface3D &surface, vector<atomType> &surfAtoms, int x, int y, int z, int type);
bool selAtom(surface3D &surface, vector<atomType> &surfAtoms, allSoseds &neighbs, int z_min, atomsCoords &tA, const vector<bool> &mask, const float *rates);
bool selAtomCA(surface3D &surface, vector<atomType> &surfAtoms, int z_min, atomsCoords &tA, vector<bool> &mask, float* rates);

atomsCoords findCell(int h, int k, int l, float &xs, float &ys, float &zs, Coords3D &Vx, Coords3D &Vy, Coords3D &Vz); //Поиск эл ячейки
void addLayer(surface3D &surface, const allSoseds& sosedi, int sX, int sY, int sZ);

atomsCoords atomsInBox(const atomsCoords &atoms, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz, const Coords3D &P1);
void findZmin(const surface3D &surface, int &zm);
void optimizeSurface(surface3D &surface, int z_min);
void delAtom(surface3D &surface, vector<atomType> &surfAtoms, int x, int y, int z, int type, int surfAtN);
bool operator==(const atomType &a1, const atomType &a2);

bool cmp_float(double x, double y);
