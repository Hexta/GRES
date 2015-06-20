/******************************************************************************
 * Copyright (c) 2009-2014 Artur Molchanov <artur.molchanov@gmail.com>        *
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

#include "Cell.h"

#include "Coords3D.h"

#include "functions.h"

#include <cmath>
#include <algorithm>

using std::pow;

namespace {
float distance(const double& x1,
    const double& y1,
    const double& z1,
    const double& x2,
    const double& y2,
    const double& z2) {
    auto const result = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) +
        pow((z2 - z1), 2));

    return static_cast<float>(result);
}

float distance(AtomInfo const& left, AtomInfo const& right) {
    auto const& leftCoords = left.type.coords;
    auto const& rightCoords = right.type.coords;

    return distance(leftCoords.x, leftCoords.y, leftCoords.z,
        rightCoords.x, rightCoords.y, rightCoords.z);
}

struct Length {
    Length(Coords3D const& left, Coords3D const& right) :
        x1(left.x),
        y1(left.y),
        z1(left.z),
        x2(right.x),
        y2(right.y),
        z2(right.z),
        l(distance(x1, y1, z1, x2, y2, z2)) {
    }

    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float l;
};

typedef std::vector<Length> Lengths;

struct rectangle {
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float x3;
    float y3;
    float z3;
    float x4;
    float y4;
    float z4;
    float S;
};

typedef std::vector<rectangle> rectangles;

// Silicon cell
Coords3DList atomTypes = {
    {1.0, 1.0, 1.0},
    {1.0, 0.5, 0.5},
    {0.5, 1.0, 0.5},
    {0.5, 0.5, 1.0},
    {0.25, 0.25, 0.25},
    {0.75, 0.75, 0.25},
    {0.75, 0.25, 0.75},
    {0.25, 0.75, 0.75}
};

bool rect_comp(const rectangle& r1, const rectangle& r2) {
    const double S1 = distance(r1.x1, r1.y1, r1.z1, r1.x2, r1.y2, r1.z2) *
                      distance(r1.x3, r1.y3, r1.z3, r1.x4, r1.y4, r1.z4);
    const double S2 = distance(r2.x1, r2.y1, r2.z1, r2.x2, r2.y2, r2.z2) *
                      distance(r2.x3, r2.y3, r2.z3, r2.x4, r2.y4, r2.z4);
    return S1 < S2;
}
}

class Cell::Private {
public:
    Private() {
    }

    Private(Atoms const& atoms) :
        m_atoms(atoms) {
    }

    Coords3D getp1() const {
        return m_p1;
    }

    Coords3D getVx() const {
        return m_vX;
    }

    Coords3D getVy() const {
        return m_vY;
    }

    Coords3D getVz() const {
        return m_vZ;
    }

    void setp1(Coords3D const& p1) {
        m_p1 = p1;
    }

    void setVx(Coords3D const& vX) {
        m_vX = vX;
    }

    void setVy(Coords3D const& vY) {
        m_vY = vY;
    }

    void setVz(Coords3D const& vZ) {
        m_vZ = vZ;
    }

    Atoms getAtoms() const {
        return m_atoms;
    }

    Atoms& getAtoms() {
        return m_atoms;
    }

    size_t size() const {
        return m_atoms.size();
    }

    const AllNeighbors findNeighbors(float xs, float ys, float zs) {
        const auto numberOfAtomsInCell = size();

        AllNeighbors totalNeighbors;
        for (size_t type0 = 0; type0 < numberOfAtomsInCell; ++type0) {
            Neighbors neighbors;

            auto const& coords0 = m_atoms[type0].type.coords;

            double const x0 = coords0.x;
            double const y0 = coords0.y;
            double const z0 = coords0.z;

            // shift of the cell #1
            for (int xc1 = -1; xc1 < 2; ++xc1)
                for (int yc1 = -1; yc1 < 2; ++yc1)
                    for (int zc1 = -1; zc1 < 2; ++zc1)
                        for (unsigned char type1 = 0;
                             type1 < numberOfAtomsInCell; ++type1)
                            if (!(type1 == type0 && !xc1 && !yc1 && !zc1)) {
                                auto const& coords1 = m_atoms[type1].type.coords;

                                double const x1 = coords1.x + xc1 * xs;
                                double const y1 = coords1.y + yc1 * ys;
                                double const z1 = coords1.z + zc1 * zs;

                                // square of the bond #1 length
                                if (fabs(pow(x1 - x0,
                                            2) + pow(y1 - y0, 2) + pow(z1 - z0, 2) - 3.0 / 16.0) <= 0.001) {
                                    AtomType s1(xc1, yc1, zc1, type1, false);
                                    neighbors.push_back(s1);
                                }
                            }
            totalNeighbors.push_back(neighbors);
        }

        return totalNeighbors;
    }

    void optimize() {
        Atoms newAtoms;

        for (auto const& atom : m_atoms) {
            const auto& coords = atom.type.coords;
            const auto x = coords.x;
            const auto y = coords.y;
            const auto z = coords.z;

            if (x > 0.01 && y > 0.01 && z > 0.01) {
                newAtoms.push_back(atom);
            }
        }

        m_atoms = newAtoms;
    }

    void moveCoords() {
        for (auto atom_coords_it = m_atoms.begin();
            atom_coords_it < m_atoms.end(); ++atom_coords_it) {

            auto& coords = atom_coords_it->type.coords;

            const auto& vec = coords - m_p1;
            const double x = (vec * m_vX) / m_vX.length();
            const double y = (vec * m_vY) / m_vY.length();
            const double z = (vec * m_vZ) / m_vZ.length();

            coords.x = static_cast<float>(x);
            coords.y = static_cast<float>(y);
            coords.z = static_cast<float>(z);
        }
    }

    void addAtom(AtomInfo const& atom) {
        m_atoms.push_back(atom);
    }

private:
    Atoms m_atoms;

    Coords3D m_p1;

    Coords3D m_vX;
    Coords3D m_vY;
    Coords3D m_vZ;
};

Cell::Cell() :
    m_impl(new Private) {
}

Cell::Cell(const Atoms& atoms) :
    m_impl(new Private(atoms)) {
}

Cell::Cell(int h, int k, int l) :
    m_impl(new Private) {
    const int SIZE_X = 5;
    const int SIZE_Y = 5;
    const int SIZE_Z = 5;

    // Найдем свободный член в уравнении секущей плоскости hx+ky+lz-C=0
    int C = (h * (SIZE_X - 1) + k * (SIZE_Y - 1) + l * (SIZE_Z - 1)) / 2 + 1;

    // atoms located on the surface
    Atoms atomsP1;

    // Find atoms located on plane #1
    atomsP1.reserve(SIZE_Z * SIZE_Y * SIZE_X * atomTypes.size());
    for (int z = 0; z < SIZE_Z; ++z)
        for (int y = 0; y < SIZE_Y; ++y)
            for (int x = 0; x < SIZE_X; ++x)
                for (auto const& atomType : atomTypes)
                    if (cmp_float(h *
                            (x + atomType.x) + k * (y + atomType.y) + l * (z + atomType.z), C)) {
                        Coords3D atom = {x + atomType.x, y + atomType.y, z + atomType.z};
                        atomsP1.push_back(AtomInfo(AtomType(atom)));
                    }

    // Find rectangles located on the plane #1
    const size_t atomsP1Size = atomsP1.size();

    Lengths lengths;
    lengths.reserve(atomsP1Size * atomsP1Size / 2);

    for (size_t i = 0; i < atomsP1Size; ++i)
        for (size_t j = i + 1; j < atomsP1Size; ++j) {
            auto const& atomP1A = atomsP1[i];
            auto const& atomP1B = atomsP1[j];
            Length L(atomP1A.type.coords, atomP1B.type.coords);
            lengths.push_back(L);
        }

    float kx = -1;
    float ky = -10;
    float kz = -100;

    rectangles rectanglesP1;

    for (auto i = lengths.begin(); i != lengths.end(); ++i)
        for (auto j = i + 1; j != lengths.end(); ++j) {
            auto const& length1 = *i;
            auto const& length2 = *j;

            auto x1 = length1.x1;
            auto y1 = length1.y1;
            auto z1 = length1.z1;
            auto x2 = length1.x2;
            auto y2 = length1.y2;
            auto z2 = length1.z2;

            auto x3 = length2.x1;
            auto y3 = length2.y1;
            auto z3 = length2.z1;
            auto x4 = length2.x2;
            auto y4 = length2.y2;
            auto z4 = length2.z2;

            auto l1 = length1.l;
            auto l2 = length2.l;

            if (cmp_float(l1, l2)) {
                if (!cmp_float(x4, x3))
                    kx = (x2 - x1) / (x4 - x3);
                else if (cmp_float(x2, x1))
                    kx = 1;

                if (!cmp_float(y4, y3))
                    ky = (y2 - y1) / (y4 - y3);
                else if (cmp_float(y2, y1))
                    ky = 1;

                if (!cmp_float(z4, z3))
                    kz = (z2 - z1) / (z4 - z3);
                else if (cmp_float(z2, z1))
                    kz = 1;

                if ((cmp_float(kx, ky) && cmp_float(kx, kz))
                    && cmp_float((x2 - x1) * (x4 - x2) + (y2 - y1) * (y4 - y2) + (z2 - z1) * (z4 - z2), 0)) {
                    rectangle rect = {x1, y1, z1, x2, y2, z2,
                                      x3, y3, z3, x4, y4, z4, l1};
                    rectanglesP1.push_back(rect);

                    // remove duplicates rectangles edges
                    for (auto v = j + 1; v != lengths.end(); ++v) {
                        const auto& v_x1 = v->x1;
                        const auto& v_x2 = v->x2;
                        const auto& v_y1 = v->y1;
                        const auto& v_y2 = v->y2;
                        const auto& v_z1 = v->z1;
                        const auto& v_z2 = v->z2;

                        if (((cmp_float(v_x1, x1) && cmp_float(v_x2, x3)
                              && cmp_float(v_y1, y1) && cmp_float(v_y2, y3)
                              && cmp_float(v_z1, z1) && cmp_float(v_z2, z3))

                             || (cmp_float(v_x1, x2) && cmp_float(v_x2, x4)
                                 && cmp_float(v_y1, y2) && cmp_float(v_y2, y4)
                                 && cmp_float(v_z1, z2) && cmp_float(v_z2, z4)))

                            || (cmp_float(v_x1, x3) && cmp_float(v_x2, x1)
                                && cmp_float(v_y1, y3) && cmp_float(v_y2, y1)
                                && cmp_float(v_z1, z3) && cmp_float(v_z2, z1))

                            || (cmp_float(v_x1, x4) && cmp_float(v_x2, x2)
                                && cmp_float(v_y1, y4) && cmp_float(v_y2, y2)
                                && cmp_float(v_z1, z4) && cmp_float(v_z2, z2))

                            ) {
                            const size_t n = v - lengths.begin();

                            std::swap(lengths[n], lengths.back());
                            lengths.pop_back();

                            v = lengths.begin() + n - 1;

                            break;
                        }
                    }
                }
            }
        }

    std::stable_sort(rectanglesP1.begin(), rectanglesP1.end(), rect_comp);

    Cells allCells;
    Atoms allAtoms = AtomsHelper::createAllCellAtoms(Coords3DList(atomTypes));

    for (auto const& rectangle : rectanglesP1) {
        auto const x1 = rectangle.x1; auto const y1 = rectangle.y1; auto const z1 = rectangle.z1;
        auto const x2 = rectangle.x2; auto const y2 = rectangle.y2; auto const z2 = rectangle.z2;
        auto const x3 = rectangle.x3; auto const y3 = rectangle.y3; auto const z3 = rectangle.z3;
        auto const x4 = rectangle.x4; auto const y4 = rectangle.y4; auto const z4 = rectangle.z4;

        for (float n = 0.5; n < 5; n += 0.5) {
            int atomsOnRectangleCount = 0;

            for (auto const& atom : allAtoms) {
                auto const& atomCoords = atom.type.coords;

                auto const atomCoordX = atomCoords.x - n * h;
                auto const atomCoordY = atomCoords.y - n * k;
                auto const atomCoordZ = atomCoords.z - n * l;

                if ((cmp_float(atomCoordX, x1)
                     && cmp_float(atomCoordY, y1)
                     && cmp_float(atomCoordZ, z1))

                    || (cmp_float(atomCoordX, x2)
                        && cmp_float(atomCoordY, y2)
                        && cmp_float(atomCoordZ, z2))

                    || (cmp_float(atomCoordX, x3)
                        && cmp_float(atomCoordY, y3)
                        && cmp_float(atomCoordZ, z3))

                    || (cmp_float(atomCoordX, x4)
                        && cmp_float(atomCoordY, y4)
                        && cmp_float(atomCoordZ, z4))) {
                    ++atomsOnRectangleCount;

                    if (atomsOnRectangleCount == 4) {
                        break;
                    }
                }
            }

            if (atomsOnRectangleCount == 4) {
                Coords3D Vx, Vy, Vz;

                Coords3D const P1(x3, y3, z3);
                Coords3D const P2(x4, y4, z4);
                Coords3D const P3(x1, y1, z1);
                Coords3D const P4(x3 + n * h, y3 + n * k, z3 + n * l);

                Vz = P4 - P1;
                Vy = P3 - P1;
                Vx = P2 - P1;

                Atoms cellAtoms(allAtoms, Vx, Vy, Vz, P1);

                Cell cell(cellAtoms);
                cell.setp1(P1);
                cell.setVx(Vx);
                cell.setVy(Vy);
                cell.setVz(Vz);

                allCells.push_back(cell);
            }
        }
    }

    stable_sort(allCells.begin(), allCells.end());

    for (auto const& cell : allCells) {
        const size_t cellSize = cell.size();

        // calculate coordination vectors
        const Coords3D P1 = cell.getp1();
        const Coords3D Vx = cell.getVx();
        const Coords3D Vy = cell.getVy();
        const Coords3D Vz = cell.getVz();

        if (allAtoms.checkContains((cell + Vx).getAtoms()) // транслируем по OX
            && allAtoms.checkContains((cell + Vy).getAtoms()) // транслируем по OY
            && allAtoms.checkContains((cell + Vz).getAtoms()) // транслируем по OZ
            && allAtoms.checkContains((cell + (-1) * Vz).getAtoms()) // транслируем по -OZ
            && allAtoms.checkContains((cell + (-1) * Vy).getAtoms()) // транслируем по -OY
            && allAtoms.checkContains((cell + (-1) * Vx).getAtoms()) // транслируем по -OX

            && Atoms(allAtoms, Vx, Vy, Vz, P1 + Vx).size() == cellSize
            && Atoms(allAtoms, Vx, Vy, Vz, P1 + Vy).size() == cellSize
            && Atoms(allAtoms, Vx, Vy, Vz, P1 + Vz).size() == cellSize
            && Atoms(allAtoms, Vx, Vy, Vz, P1 + -1 * Vx).size() == cellSize
            && Atoms(allAtoms, Vx, Vy, Vz, P1 + -1 * Vy).size() == cellSize
            && Atoms(allAtoms, Vx, Vy, Vz, P1 + -1 * Vz).size() == cellSize) {
            *this = cell;

            break;
        }
    }

    m_impl->moveCoords();
}

Cell::Cell(const Cell& other) :
    m_impl(new Private(*other.m_impl)) {
}

Cell::Cell(Cell&& other) :
    m_impl(std::move(other.m_impl)) {
}

size_t Cell::size() const {
    return m_impl->size();
}

bool Cell::operator<(const Cell& cell) const {
    return size() < cell.size();
}

const AllNeighbors Cell::findNeighbors(float xs, float ys, float zs) {
    return m_impl->findNeighbors(xs, ys, zs);
}

void Cell::optimize() {
    m_impl->optimize();
}

double Cell::getXSize() const {
    return m_impl->getVx().length();
}

double Cell::getYSize() const {
    return m_impl->getVy().length();
}

double Cell::getZSize() const {
    return m_impl->getVz().length();
}

Coords3D Cell::getVx() const {
    return m_impl->getVx();
}

Coords3D Cell::getVy() const {
    return m_impl->getVy();
}

Coords3D Cell::getVz() const {
    return m_impl->getVz();
}

Coords3D Cell::getp1() const {
    return m_impl->getp1();
}

void Cell::setp1(Coords3D const& p1) {
    m_impl->setp1(p1);
}

void Cell::setVx(Coords3D const& vX) {
    m_impl->setVx(vX);
}

void Cell::setVy(Coords3D const& vY) {
    m_impl->setVy(vY);
}

void Cell::setVz(Coords3D const& vZ) {
    m_impl->setVz(vZ);
}

Cell& Cell::operator=(const Cell& other) {
    m_impl.reset(new Private(*other.m_impl));

    return *this;
}

Cell& Cell::operator=(Cell&& other) {
    m_impl = std::move(other.m_impl);

    return *this;
}

Cell::~Cell() {
}

Atoms Cell::getAtoms() const {
    return m_impl->getAtoms();
}

Atoms& Cell::getAtoms() {
    return m_impl->getAtoms();
}

void Cell::addAtom(AtomInfo const& atom) {
    m_impl->addAtom(atom);
}

const Cell operator+(Cell const& cell, Coords3D const& v) {
    Cell resultCell;

    Atoms const atoms = cell.getAtoms();

    for (auto cell_atom_it = atoms.begin();
         cell_atom_it < atoms.end(); ++cell_atom_it) {
        resultCell.addAtom(AtomInfo(AtomType(cell_atom_it->type.coords + v)));
    }

    return resultCell;
}
