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
#include "Coords3D.h"

#include <cstddef>
#include <vector>
#include <deque>
#include <memory>

class Cell {
public:
    Cell();
    explicit Cell(const Atoms& atoms);
    Cell(const Cell& other);
    Cell(Cell&& other);

    Cell(int h, int k, int l);

    ~Cell();

    bool operator<(const Cell& cell) const;

    Cell& operator=(const Cell& other);
    Cell& operator=(Cell&& other);

    size_t size() const;

    const AllNeighbors findNeighbors(float xs, float ys, float zs);

    // remove shared atoms
    void optimize();

    double getXSize() const;
    double getYSize() const;
    double getZSize() const;

    Atoms getAtoms() const;
    Atoms& getAtoms();

    void addAtom(AtomInfo const& atom);

    // смещение координат в €чейку
    void moveCoords(const Coords3D &O, const Coords3D &Vx, const Coords3D &Vy,
        const Coords3D &Vz);

    Coords3D getp1() const;
    Coords3D getVx() const;
    Coords3D getVy() const;
    Coords3D getVz() const;

    void setp1(Coords3D const& p1);
    void setVx(Coords3D const& vX);
    void setVy(Coords3D const& vX);
    void setVz(Coords3D const& vX);

private:
    class Private;
    std::unique_ptr<Private> m_impl;
};

const Cell operator+(Cell const& cell, Coords3D const& v);

typedef std::deque<Cell> Cells;
