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

#include "Surface3D.h"

#include <algorithm>


Surface3D::Surface3D() {
}

Surface3D::Surface3D(std::vector<Surface2D> const& surfaces2D) : m_surfaces2D(surfaces2D) {
}


Surface3D::~Surface3D() {
}

Surface2D& Surface3D::operator[](size_t n) {
    return m_surfaces2D[n];
}

Surface2D const& Surface3D::operator[](size_t n) const {
    return m_surfaces2D[n];
}

size_t Surface3D::size() const {
    return m_surfaces2D.size();
}

void Surface3D::reserve(size_t n) {
    m_surfaces2D.reserve(n);
}

int Surface3D::findZmin(int zBegin) {
    for (auto z = m_surfaces2D.begin() + zBegin; z != m_surfaces2D.end(); ++z)
    for (auto y = z->begin() + 2; y != z->end() - 1; ++y) {
        const auto yEnd = y->end();
        for (auto x = y->begin() + 2; x != yEnd - 1; ++x)

        for (auto &atom : *x)
        if (!atom.deleted) {
            return static_cast<int> (z - m_surfaces2D.begin());
        }
    }

    return zBegin;
}

void Surface3D::optimize(int zMin) {
    if (zMin > 1)
        m_surfaces2D[zMin - 2].clear();
    m_surfaces2D.shrink_to_fit();
}

void Surface3D::clear() {
    m_surfaces2D.clear();
}

void Surface3D::push_back(Surface2D const& surface2D) {
    m_surfaces2D.push_back(surface2D);
}

void Surface3D::recallNeighbors(int x, int y, int z, int type,
    AtomTypes& surfAtoms) {
    for (auto &neighb : m_surfaces2D[z][y][x][type].neighbors) {
        const int xNb = neighb.x;
        const int yNb = neighb.y;
        const int zNb = neighb.z;
        const unsigned char typeNb = neighb.type;

        AtomInfo &neihgbAtomInfo = m_surfaces2D[zNb][yNb][xNb][typeNb];

        size_t i = 0;
        for (auto &neighb_int : neihgbAtomInfo.neighbors) {
            if (neighb_int.x == x && neighb_int.y == y
                && neighb_int.z == z
                && neighb_int.type == type) {
                neihgbAtomInfo.neighbors.erase(neihgbAtomInfo.neighbors.begin() + i);
                neihgbAtomInfo.fNbCount -= 1;
                break;
            }
            ++i;
        }
        if (!neihgbAtomInfo.fNbCount) {
            neihgbAtomInfo.deleted = true;
            neihgbAtomInfo.neighbors.clear();
        }

        if ((neihgbAtomInfo.fNbCount == 3) &&
            ((xNb > 1 && xNb < m_surfaces2D[zNb][yNb].size() - 2)
            && (yNb > 1 && yNb < m_surfaces2D[zNb].size() - 2))) {

            AtomType aT = {xNb, yNb, zNb, typeNb, false};
            m_surfaceAtoms.push_back(aT);
            surfAtoms.push_back(aT);
        }
    }
}

const AtomTypes& Surface3D::getSurfaceAtoms() const {
    return m_surfaceAtoms;
}

void Surface3D::delAtom(AtomTypes& surfAtoms, int x, int y, int z, int type,
    int surfAtN) {
    if (surfAtN != -1) {

        std::swap(surfAtoms[surfAtN], surfAtoms.back());
        surfAtoms.pop_back();
        std::swap(surfAtoms[surfAtN], surfAtoms.back());
    }

    recallNeighbors(x, y, z, type, surfAtoms);

    m_surfaces2D[z][y][x][type].fNbCount = 0;
    m_surfaces2D[z][y][x][type].deleted = true;
    m_surfaces2D[z][y][x][type].neighbors.clear();
}
