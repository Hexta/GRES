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

#include <cstddef>
#include <memory>
#include <vector>

typedef std::vector<CellInfo> Surface1D;
typedef std::vector<Surface1D> Surface2D;

class Cell;

class Surface3D {
public:
    Surface3D();
    Surface3D(std::vector<Surface2D> const& surfaces2D);
    ~Surface3D();

    Surface2D& operator[](size_t n);
    Surface2D const& operator[](size_t n) const;
    size_t size() const;
    void reserve(size_t n);

    int findZmin(int zBegin);
    void optimize(int zMin);
    void clear();
    void push_back(Surface2D const& surface2D);
    void recallNeighbors(int x, int y, int z, int type);

    const AtomTypes& getSurfaceAtoms() const;
    void delAtom(int x, int y, int z, int type, int surfAtN);

    bool selAtom(AllNeighbors const& neighbs, size_t z_min, Cell &tA,
        const std::vector<bool> &mask, const float *rates);

    bool selAtomCA(int z_min, Cell &tA, std::vector<bool> const& mask,
        float* rates);

    void addLayer(const AllNeighbors& totalNeighbors, int sX, int sY, int sZ);
    void rebuildSurfaceAtoms();

private:
    AtomTypes initSurfaceAtoms() const;
    std::vector<Surface2D> m_surfaces2D;
    AtomTypes m_surfaceAtoms;
};

typedef std::shared_ptr<Surface3D> Surface3DPtr;
