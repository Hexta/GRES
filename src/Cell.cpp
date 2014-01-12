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

#include "Cell.h"

#include "Coords3D.h"

#include <cmath>

using std::pow;

Cell::Cell() {

}

Cell::Cell(const std::vector<Coords3D>& atoms) : atoms(atoms) {

}

bool Cell::operator==(const Cell& cell) const {
    return atoms.size() < cell.atoms.size();
}

size_t Cell::size() const {
    return atoms.size();
}

bool Cell::operator<(const Cell& cell) const {
    return atoms.size() < cell.atoms.size();
}

const AllNeighbors Cell::findSoseds(float xs, float ys, float zs) {
    const int NUMBER_OF_ATOMS_IN_CELL = static_cast<int> (atoms.size());
    AllNeighbors allSosedi;
    for (int type0 = 0; type0 < NUMBER_OF_ATOMS_IN_CELL; ++type0) {
        Neighbors sosedi;

        double x0 = atoms[type0].x;
        double y0 = atoms[type0].y;
        double z0 = atoms[type0].z;

        for (int xc1 = -1; xc1 < 2; ++xc1)//смещение соседней €чейки є1
        for (int yc1 = -1; yc1 < 2; ++yc1)
        for (int zc1 = -1; zc1 < 2; ++zc1)
        for (unsigned char type1 = 0;
            type1 < NUMBER_OF_ATOMS_IN_CELL; ++type1)
        if (!(type1 == type0 && !xc1 && !yc1 && !zc1)) {
            double x1 = atoms[type1].x + xc1*xs;
            double y1 = atoms[type1].y + yc1*ys;
            double z1 = atoms[type1].z + zc1*zs;
            if (fabs(pow(x1 - x0, 2) + pow(y1 - y0, 2) + pow(z1 - z0, 2) - 3.0 / 16.0) <= 0.001) //квадрат длины св€зи #1
            {
                AtomType s1 = {xc1, yc1, zc1, type1, false};
                sosedi.push_back(s1);
            }
        }
        allSosedi.push_back(sosedi);
    }

    return allSosedi;
}

void Cell::optimize() {
    Atoms newAtoms;

    for (auto i = atoms.begin(); i < atoms.end() - 4; ++i) {
        const float x = i->x;
        const float y = i->y;
        const float z = i->z;

        if (x > 0.01 && y > 0.01 && z > 0.01) {
            Coords3D atom = {x, y, z};
            newAtoms.push_back(atom);
        }
    }

    atoms = newAtoms;
}

double Cell::getXSize() const
{
    return (*(atoms.end() - 3)).length();
}

double Cell::getYSize() const
{
    return (*(atoms.end() - 2)).length();
}

double Cell::getZSize() const
{
    return (*(atoms.end() - 1)).length();
}
