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
    AtomType() {
    }

    AtomType(int x, int y, int z, uint8_t type, bool deleteFlag) :
        x(x),
        y(y),
        z(z),
        type(type),
        toDel(deleteFlag) {
    }

    explicit AtomType(Coords3D const& coords) :
        coords(coords) {
    }

    Coords3D coords;

    int x; //сдвиг €чейки по OX (-1;0;+1)
    int y;
    int z;

    // atom type (1 - 8)
    uint8_t type;
    bool toDel;
};

bool operator==(const AtomType &a1, const AtomType &a2);

typedef std::vector<AtomType> AtomTypes;
 // atom's neighbors
typedef std::vector<AtomType> Neighbors;
 // all neighbors
typedef std::vector<Neighbors> AllNeighbors;

struct AtomInfo {
    AtomInfo() {
    }

    explicit AtomInfo(AtomType const& type) :
        type(type) {
    }

    AtomType type;
    std::vector<AtomType> neighbors;
    char fNbCount;
    bool deleted;
};

class Atoms
{
public:
    typedef std::vector<AtomInfo>::iterator Iterator;
    typedef std::vector<AtomInfo>::const_iterator ConstIterator;

public:
    Atoms();
    Atoms(Atoms const& other);
    Atoms(Atoms&& other);
    Atoms(std::size_t, AtomInfo const& atom);

    Atoms(const Atoms &atomsIn, const Coords3D &Vx, const Coords3D &Vy,
        const Coords3D &Vz, const Coords3D &P1);

    ~Atoms();

    Atoms& operator=(Atoms&& other);
    Atoms& operator=(Atoms const& other);

    void reserve(std::size_t size);
    void push_back(AtomInfo const& atom);

    Iterator begin();
    ConstIterator begin() const;

    Iterator end();
    ConstIterator end() const;

    std::size_t size() const;
    AtomInfo& operator[](size_t n);

    void clear();
    bool empty() const;

    bool checkContains(Atoms const& other);

private:
    class Private;
    std::unique_ptr<Private> m_impl;
};

class AtomsHelper {
public:
    static Atoms createAllCellAtoms(Coords3DList const& atomTypes);

private:
    static Atoms allCellAtoms;
};
