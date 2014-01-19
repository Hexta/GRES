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
#include "Cell.h"
#include "Surface3D.h"

#include <vector>
#include <cstdlib>

static float P1 = 1.0f;
static float P2 = 0.030f;
static float P3 = 0.00010f;

struct Coords3D;

Cell findCell(int h, int k, int l, float &xs, float &ys, float &zs, Coords3D &Vx, Coords3D &Vy, Coords3D &Vz); //Поиск эл ячейки
Cell atomsInBox(const Atoms &atoms, const Coords3D &Vx, const Coords3D &Vy, const Coords3D &Vz, const Coords3D &P1);

bool cmp_float(double x, double y);
