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

#include <cstddef>
#include <vector>
#include <memory>

struct AtomType {
    int x; //сдвиг €чейки по OX (-1;0;+1)
    int y;
    int z;
    unsigned char type; //тип атома (1 -- 8)
    bool toDel;
};

bool operator==(const AtomType &a1, const AtomType &a2);

typedef std::vector<AtomType> AtomTypes;
 // atom's neighbors
typedef std::vector<AtomType> Neighbors;
 // all neighbors
typedef std::vector<Neighbors> AllNeighbors;

struct AtomInfo {
    Neighbors neighbors;
    char fNbCount;
    bool deleted;
};

typedef std::vector<Coords3D> Atoms;
typedef std::shared_ptr<Atoms> AtomsPtr;
typedef std::vector<AtomInfo> CellInfo;
