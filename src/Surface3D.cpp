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
#include "Cell.h"

#include <algorithm>
#include <random>
#include <functional>

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

bool Surface3D::selAtom(AtomTypes& surfAtoms, AllNeighbors &neighbs, int z_min,
    Cell &tA, const std::vector<bool> &mask, const float *rates, float& P1,
    float& P2, float& P3) {
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = static_cast<int> (m_surfaces2D[z_min].size());
    const int xC = static_cast<int> (m_surfaces2D[z_min][0].size());
    bool result = false;
    bool maskON = !mask.empty();

    static std::mt19937 rdevice(std::random_device{}
    ());
    static std::uniform_real_distribution<float> dis;
    static auto rd = bind(dis, rdevice);

    int i = static_cast<int> (std::uniform_int_distribution<size_t> {0, surfAtoms.size() - 1}
    (rdevice));

    auto &surfAtom = surfAtoms[i];
    int x = surfAtom.x;
    int y = surfAtom.y;
    int z = surfAtom.z;
    unsigned char a = surfAtom.type;
    auto &surfaceZY = m_surfaces2D[z][y];
    if ((maskON && ((z == 0
        && !mask[(y - 2)*(surfaceZY.size() - 4) + x - 2])
        || (z + tA.atoms[a].z >= 0.5)))
        || !maskON) {
        float randN = rd();
        int bonds = surfaceZY[x][a].fNbCount;
        if (!bonds) {
            std::swap(surfAtoms[i], surfAtoms.back());
            surfAtoms.pop_back();
            std::swap(surfAtoms[i], surfAtoms.back());
            --i;
            return false;
        }

        if ((bonds == 1 && randN < P1)
            || (bonds == 2 && randN < P2)
            || (bonds == 3 && randN < P3)) {
            delAtom(surfAtoms, x, y, z, a, i);
            /*
            удалим и соостветсвенный краевой атом
            */
            if (x == 4 || x == 5) {
                delAtom(surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                delAtom(surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                delAtom(surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                delAtom(surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }
            if (z >= size() - 3) {
                addLayer(neighbs, xC, yC, static_cast<int>(size()));
            }
        }
    }

    return result;
}

void Surface3D::addLayer(const AllNeighbors& sosedi, int sX, int sY, int sZ) {
    Surface2D surfaceXY;
    surfaceXY.reserve(sY);

    for (int y = 0; y < sY; ++y) {
        Surface1D surfaceX;
        surfaceX.reserve(sX);

        for (int x = 0; x < sX; ++x) {
            CellInfo cell;
            cell.reserve(sosedi.size());
            for (auto &sosed : sosedi) {
                Neighbors neighbs;
                char numberNeighbs = 0; //Число первых соседей
                for (int nb = 0; nb < 4; ++nb) {
                    if (x + sosed[nb].x >= 0 && y + sosed[nb].y >= 0 && x + sosed[nb].x < sX && y + sosed[nb].y < sY) {
                        ++numberNeighbs;
                        AtomType neighb = {x + sosed[nb].x, y + sosed[nb].y, sZ + sosed[nb].z, sosed[nb].type, false};
                        neighbs.push_back(neighb);
                    }
                }
                AtomInfo atom = {neighbs, numberNeighbs, !numberNeighbs};
                cell.push_back(atom);
            }
            surfaceX.push_back(cell);
        }
        surfaceXY.push_back(surfaceX);
    }

    push_back(surfaceXY);
}

bool Surface3D::selAtomCA(AtomTypes& surfAtoms, int z_min, Cell &tA, std::vector<bool> &mask, float* rates, float& P1, float& P2, float& P3) {
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = static_cast<int> (m_surfaces2D[z_min].size());
    const int xC = static_cast<int> (m_surfaces2D[z_min][0].size());
    bool result = false;
    bool maskON = !mask.empty();

    size_t surfAtomsSize = surfAtoms.size();

    static std::mt19937 rdevice(std::random_device{}
    ());
    static std::uniform_real_distribution<float> dis;
    static auto rd = bind(dis, rdevice);

    for (unsigned int i = 0; i < surfAtomsSize; ++i) {
        const auto &surfAtom = surfAtoms[i];
        int x = surfAtom.x;
        int y = surfAtom.y;
        unsigned int z = surfAtom.z;
        unsigned short a = surfAtom.type;
        const auto &tAz = tA.atoms[a].z;
        auto &surfaceZY = m_surfaces2D[z][y];
        if ((maskON && (((z + tAz < 0.5) && !mask[(y - 2)*(surfaceZY.size() - 4) + x - 2])
            || (z + tAz >= 0.5)))
            || !maskON) {
            float randN = rd();
            int bonds = surfaceZY[x][a].fNbCount;
            if (!bonds) {
                surfAtoms.erase(surfAtoms.begin() + i--);
                --surfAtomsSize;
                continue;
            }

            if ((bonds == 1 && randN < P1)
                || (bonds == 2 && randN < P2)
                || (bonds == 3 && randN < P3)) {
                surfAtoms[i].toDel = true;
                if (z == size() - 3)
                    result = true;
            }
        }
    }

    int i = 0;
    auto surfAtomsEnd = surfAtoms.end();
    for (auto surfAtomIter = surfAtoms.begin(); surfAtomIter != surfAtomsEnd;
        ++surfAtomIter) {
        if (surfAtomIter->toDel) {
            int x = surfAtomIter->x;
            int y = surfAtomIter->y;
            int z = surfAtomIter->z;
            unsigned char a = surfAtomIter->type;

            delAtom(surfAtoms, x, y, z, a, i);

            if (x == 4 || x == 5) {
                delAtom(surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                delAtom(surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                delAtom(surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                delAtom(surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }

            surfAtomIter = surfAtoms.begin() + --i;
            surfAtomsEnd = surfAtoms.end();
        }
        ++i;
    }
    return result;
}
