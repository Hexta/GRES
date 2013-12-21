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

#ifndef GEOMETRY
#define GEOMETRY

#include "functions.h"
#include "dataTypes.h"

void createAtomsAndBondes(surface3D &surface, vector<atomType>&,
        atomsCoords &cellAts, float xs_, float ys_, float zs_, int z_min,
        float scaling, vector<atomName> &atNames_, Bonds &outBonds);

void createSphere(GLdouble radius, GLint slices, GLint stacks, int &vSize1,
        int &vSize2, int &vSize3);
void normalize(float v[3]);
void normalize(float v[3], coords3D&);
coords3D normalize(coords3D);
void norm(coords3D &in);
coords3D normcrossprod(coords3D in1, coords3D in2);

#endif
