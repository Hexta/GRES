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
#include "Surface3D.h"

#include <vector>
#include <cstdlib>

using std::vector;

bool compareTranslCell(const Cell &cellAtoms, const Cell &allAtoms, const Coords3D &V); //Совпадение атомов после трансляции
void coordsMove(Cell &ca, const Coords3D &O, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz); // смещение координат в ячейку

void recallNeighbors(Surface3D &surface, vector<AtomType> &surfAtoms, int x, int y, int z, int type);

bool selAtom(Surface3D &surface, vector<AtomType> &surfAtoms, AllNeighbors &neighbs, int z_min, Cell &tA, const vector<bool> &mask, const float *rates);
bool selAtomCA(Surface3D &surface, vector<AtomType> &surfAtoms, int z_min, Cell &tA, vector<bool> &mask, float* rates);

Cell findCell(int h, int k, int l, float &xs, float &ys, float &zs, Coords3D &Vx, Coords3D &Vy, Coords3D &Vz); //Поиск эл ячейки
void addLayer(Surface3D &surface, const AllNeighbors& sosedi, int sX, int sY, int sZ);

Cell atomsInBox(const Atoms &atoms, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz, const Coords3D &P1);

void delAtom(Surface3D &surface, vector<AtomType> &surfAtoms, int x, int y, int z, int type, int surfAtN);

bool cmp_float(double x, double y);
