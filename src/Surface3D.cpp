/******************************************************************************
* Copyright (c) 2009-2016 Artur Molchanov <artur.molchanov@gmail.com>        *
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
#include "Cell.h"

#include <algorithm>
#include <random>
#include <functional>
#include <iostream>

Surface3D::Surface3D(Cell const& cell) :
    m_cell(cell) {
}

Surface3D::Surface3D(std::vector<Surface2D> const& surfaces2D) :
    m_surfaces2D(surfaces2D),
    m_surfaceAtoms(initSurfaceAtoms()) {
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
    for (auto z = m_surfaces2D.cbegin() + zBegin; z != m_surfaces2D.cend(); ++z)
        for (auto y = z->cbegin() + 2; y != z->cend() - 1; ++y) {
            auto const yEnd = y->cend();

            for (auto x = y->cbegin() + 2; x != yEnd - 1; ++x)

                for (auto const& atom : x->getAtoms())
                    if (!atom.deleted) {
                        return static_cast<int>(z - m_surfaces2D.cbegin());
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
    m_surfaceAtoms.clear();
}

void Surface3D::push_back(Surface2D const& surface2D) {
    m_surfaces2D.push_back(surface2D);
}

void Surface3D::recallNeighbors(int x, int y, int z, int type) {
    for (auto const& neighb : m_surfaces2D[z][y][x].getAtoms()[type].neighbors) {
        const int xNb = neighb.x;
        const int yNb = neighb.y;
        const int zNb = neighb.z;
        const unsigned char typeNb = neighb.type;

        AtomInfo& neihgbAtomInfo = m_surfaces2D[zNb][yNb][xNb].getAtoms()[typeNb];

        size_t i = 0;
        for (auto& neighb_int : neihgbAtomInfo.neighbors) {
            if (neighb_int.x == x && neighb_int.y == y
                && neighb_int.z == z
                && neighb_int.type == type) {
                neihgbAtomInfo.neighbors.erase(neihgbAtomInfo.neighbors.begin() + i);
                neihgbAtomInfo.firstNeighborsCount -= 1;
                break;
            }
            ++i;
        }
        if (!neihgbAtomInfo.firstNeighborsCount) {
            neihgbAtomInfo.deleted = true;
            neihgbAtomInfo.neighbors.clear();
        }

        if ((neihgbAtomInfo.firstNeighborsCount == 3)
            && ((xNb > 1 && xNb < static_cast<decltype(xNb)>(m_surfaces2D[zNb][yNb].size()) - 2)
                && (yNb > 1 && yNb < static_cast<decltype(yNb)>(m_surfaces2D[zNb].size()) - 2))) {
            AtomType aT = {xNb, yNb, zNb, typeNb, false};
            m_surfaceAtoms.push_back(aT);
        }
    }
}

const AtomTypes& Surface3D::getSurfaceAtoms() const {
    return m_surfaceAtoms;
}

void Surface3D::deleteAtom(int x, int y, int z, int type) {
    recallNeighbors(x, y, z, type);

    m_surfaces2D[z][y][x].getAtoms()[type].firstNeighborsCount = 0;
    m_surfaces2D[z][y][x].getAtoms()[type].deleted = true;
    m_surfaces2D[z][y][x].getAtoms()[type].neighbors.clear();
}

void Surface3D::deleteAtom(int x, int y, int z, int type, std::size_t atomNumberAtSurface) {
    std::swap(m_surfaceAtoms[atomNumberAtSurface], m_surfaceAtoms.back());
    m_surfaceAtoms.pop_back();
    std::swap(m_surfaceAtoms[atomNumberAtSurface], m_surfaceAtoms.back());

    deleteAtom(x, y, z, type);
}

bool Surface3D::isMaskedAtom(AtomType const& atom, std::vector<bool> const& mask) const {
    if (mask.empty())
    {
        return false;
    }

    auto const y = atom.y;
    auto const x = atom.x;
    auto const z = atom.z;
    auto const a = atom.type;

    auto const& surfaceZY = m_surfaces2D[z][y];

    return !((z == 0 && !mask[(y - 2) * (surfaceZY.size() - 4) + x - 2])
        || (z + m_cell.getAtoms()[a].type.coords.z >= 0.5));
}

template<typename ResultType>
std::function<ResultType()> getUniformRandomCalculator() {
    std::mt19937 randomProviderDevice(std::random_device{} ());
    std::uniform_real_distribution<float> uniformDistribution;

    return bind(uniformDistribution, randomProviderDevice);
}

bool shouldDeleteAtom(int bondsCount, float const * rates) {
    static auto calculateRandom = getUniformRandomCalculator<float>();

    auto const& P1 = rates[0];
    auto const& P2 = rates[1];
    auto const& P3 = rates[2];

    float const randN = calculateRandom();

    return (bondsCount == 1 && randN < P1)
        || (bondsCount == 2 && randN < P2)
        || (bondsCount == 3 && randN < P3);

}

bool Surface3D::deleteRandomAtomKmc(AllNeighbors const& neighbors,
    size_t z_min,
    const std::vector<bool>& mask,
    const float* rates) {
    const int yC = static_cast<int>(m_surfaces2D[z_min].size());
    const int xC = static_cast<int>(m_surfaces2D[z_min][0].size());
    bool result = false;

    static std::mt19937 randomProviderDevice(std::random_device {} ());

    auto const i = std::uniform_int_distribution<size_t>{0, m_surfaceAtoms.size() - 1} (randomProviderDevice);

    auto& surfAtom = m_surfaceAtoms[i];
    int x = surfAtom.x;
    int y = surfAtom.y;
    int z = surfAtom.z;
    unsigned char a = surfAtom.type;
    auto& surfaceZY = m_surfaces2D[z][y];

    if (isMaskedAtom(surfAtom, mask)) {
        return result;
    }

    int const bonds = surfaceZY[x].getAtoms()[a].firstNeighborsCount;

    if (bonds <= 0) {
        std::swap(m_surfaceAtoms[i], m_surfaceAtoms.back());
        m_surfaceAtoms.pop_back();
        std::swap(m_surfaceAtoms[i], m_surfaceAtoms.back());
        return false;
    }

    if (!shouldDeleteAtom(bonds, rates)) {
        return result;
    }

    deleteAtom(x, y, z, a, i);

    // Delete appropriate edge atom too
    if (x == 4 || x == 5) {
        deleteAtom(5 - x, y, z, a);
    }

    if (x == xC - 6 || x == xC - 5) {
        deleteAtom(2 * (xC - 1) - 5 - x, y, z, a);
    }

    if (y == 4 || y == 5) {
        deleteAtom(x, 5 - y, z, a);
    }

    if (y == yC - 6 || y == yC - 5) {
        deleteAtom(x, 2 * yC - 7 - y, z, a);
    }

    if (z >= static_cast<decltype(z)>(size()) - 3) {
        addLayer(neighbors, xC, yC, static_cast<int>(size()));
    }
    
    return result;
}

void Surface3D::addLayer(const AllNeighbors& totalNeighbors, int sX, int sY, int sZ) {
    Surface2D surfaceXY;
    surfaceXY.reserve(sY);

    for (int y = 0; y < sY; ++y) {
        Surface1D surfaceX;
        surfaceX.reserve(sX);

        for (int x = 0; x < sX; ++x) {
            Atoms atoms;
            atoms.reserve(totalNeighbors.size());

            for (auto const& neighbors : totalNeighbors) {
                Neighbors neighbs;
                // First neighbors count
                char numberNeighbs = 0;

                for (int nb = 0; nb < 4; ++nb) {
                    if (x + neighbors[nb].x >= 0
                        && y + neighbors[nb].y >= 0
                        && x + neighbors[nb].x < sX
                        && y + neighbors[nb].y < sY) {
                        ++numberNeighbs;
                        AtomType neighb = {x + neighbors[nb].x,
                                           y + neighbors[nb].y, sZ + neighbors[nb].z,
                                           neighbors[nb].type, false};
                        neighbs.push_back(neighb);
                    }
                }

                AtomInfo atom;
                atom.neighbors = neighbs;
                atom.firstNeighborsCount = numberNeighbs;
                atom.deleted = numberNeighbs == 0;
                atoms.push_back(atom);
            }

            surfaceX.push_back(Cell(atoms));
        }

        surfaceXY.push_back(surfaceX);
    }

    push_back(surfaceXY);
}

bool Surface3D::deleteRandomAtomCa(int z_min, std::vector<bool> const& mask, float* rates) {
    const int yC = static_cast<int>(m_surfaces2D[z_min].size());
    const int xC = static_cast<int>(m_surfaces2D[z_min][0].size());
    bool result = false;

    auto surfAtomsSize = m_surfaceAtoms.size();

    for (unsigned int i = 0; i < surfAtomsSize; ++i) {
        auto const& surfAtom = m_surfaceAtoms[i];
        auto const x = surfAtom.x;
        auto const y = surfAtom.y;
        auto const z = surfAtom.z;
        auto const a = surfAtom.type;
        auto const& surfaceZY = m_surfaces2D[z][y];

        if (isMaskedAtom(surfAtom, mask)) {
            continue;
        }

        auto const bonds = surfaceZY[x].getAtoms()[a].firstNeighborsCount;

        if (bonds == 0) {
            m_surfaceAtoms.erase(m_surfaceAtoms.begin() + i--);
            --surfAtomsSize;

            continue;
        }

        if (shouldDeleteAtom(bonds, rates))
        {
            m_surfaceAtoms[i].shouldBeDeleted = true;

            if (z == size() - 3) {
                result = true;
            }
        }
    }

    std::size_t i = 0;
    auto surfAtomsEnd = m_surfaceAtoms.end();

    for (auto surfAtomIter = m_surfaceAtoms.begin(); surfAtomIter != surfAtomsEnd;
         ++surfAtomIter) {

        if (!surfAtomIter->shouldBeDeleted) {
            ++i;

            continue;
        }

        int x = surfAtomIter->x;
        int y = surfAtomIter->y;
        int z = surfAtomIter->z;
        unsigned char a = surfAtomIter->type;

        deleteAtom(x, y, z, a, i);

        if (x == 4 || x == 5) {
            deleteAtom(5 - x, y, z, a);
        }
        if (x == xC - 6 || x == xC - 5) {
            deleteAtom(2 * (xC - 1) - 5 - x, y, z, a);
        }
        if (y == 4 || y == 5) {
            deleteAtom(x, 5 - y, z, a);
        }
        if (y == yC - 6 || y == yC - 5) {
            deleteAtom(x, 2 * yC - 7 - y, z, a);
        }

        if (i == 0) {
            i = 1;
        }

        if (i >= m_surfaceAtoms.size()) {
            surfAtomIter = m_surfaceAtoms.begin() + (m_surfaceAtoms.size() - 1);
        }
        else {
            surfAtomIter = m_surfaceAtoms.begin() + --i;
        }

        surfAtomsEnd = m_surfaceAtoms.end();

        ++i;
    }
    return result;
}

AtomTypes Surface3D::initSurfaceAtoms() const {
    AtomTypes surfaceAtoms;
    int z = 0;
    for (auto const& surface2D : m_surfaces2D) {
        int y = 0;
        for (auto const& surface1D : surface2D) {
            int x = 0;
            for (auto const& cell : surface1D) {
                unsigned char a = 0;
                for (auto const& atom : cell.getAtoms()) {
                    if (x > 1 && x < static_cast<int>(surface1D.size()) - 2
                        && y > 1 && y < static_cast<int>(surface2D[z].size()) - 2
                        && atom.neighbors.size() < 4) {
                        AtomType atomInfo = {x, y, z, a, false};
                        surfaceAtoms.push_back(atomInfo);
                    }
                    ++a;
                }
                ++x;
            }
            ++y;
        }
        ++z;
    }
    return surfaceAtoms;
}

void Surface3D::rebuildSurfaceAtoms() {
    m_surfaceAtoms = initSurfaceAtoms();
}
