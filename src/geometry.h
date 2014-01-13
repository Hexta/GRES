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

#include "functions.h"

#ifdef _WIN32 
#include <windows.h>
#include <wingdi.h>
#endif

#include <GL/glu.h>

#include <vector>

using std::vector;

struct atomName {
    int name;
    int xC;
    int yC;
    int zC;
    float x;
    float y;
    float z;
    unsigned char type;
    int fNbCount;
};

typedef vector<atomName>AtomsNames;

struct Bond {
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
};

typedef vector<Bond> Bonds;

void createAtomsAndBondes(Surface3D &surface, const vector<AtomType>& surfAtoms,
        const Cell &cellAts, float xs_, float ys_, float zs_, int z_min,
        float scaling, vector<atomName> &atNames_, Bonds &outBonds);

void createSphere(GLdouble radius, GLint slices, GLint stacks, int &vSize1,
        int &vSize2, int &vSize3);
void normalize(float v[3]);
void normalize(float v[3], Coords3D&);
Coords3D normalize(const Coords3D& in);
void norm(Coords3D &in);
Coords3D normcrossprod(const Coords3D& in1, const Coords3D& in2);
