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

#include "Atoms.h"

bool operator==(const AtomType& a1, const AtomType& a2) {
    return ((a1.x == a2.x) && (a1.y == a2.y) && (a1.z == a2.z) && (a1.type == a2.type));
}

class Atoms::Private {
public:
    Private() :
        m_atomInfos() {
    }

    Private(std::size_t n, AtomInfo const& atom) :
        m_atomInfos(n, atom) {
    }

    Private(Private const& other) :
        m_atomInfos(other.m_atomInfos) {
    }

    ~Private() {
    }

    void reserve(std::size_t size) {
        m_atomInfos.reserve(size);
    }

    void push_back(AtomInfo const& atom) {
        m_atomInfos.push_back(atom);
    }

    Iterator begin() {
        return m_atomInfos.begin();
    }

    ConstIterator begin() const {
        return m_atomInfos.begin();
    }

    Iterator end() {
        return m_atomInfos.end();
    }

    ConstIterator end() const {
        return m_atomInfos.end();
    }

    std::size_t size() const {
        return m_atomInfos.size();
    }

    AtomInfo& operator[](size_t n) {
        return m_atomInfos[n];
    }

    void clear() {
        m_atomInfos.clear();
    }

    bool empty() const {
        return m_atomInfos.empty();
    }

    bool checkContains(Private const& other) {
        bool isEual = false;
        for (auto otherAtomsIt = other.m_atomInfos.begin(); otherAtomsIt != other.m_atomInfos.end();
             ++otherAtomsIt) {
            isEual = false;

            for (auto& atom : m_atomInfos) {
                if (otherAtomsIt->type.coords == atom.type.coords) {
                    isEual = true;
                    break;
                }
            }

            if (!isEual)
                break;
        }

        return isEual;
    }

private:
    typedef std::vector<AtomInfo> AtomInfos;
    AtomInfos m_atomInfos;
};

Atoms::Atoms() :
    m_impl(new Private) {
}

Atoms::Atoms(Atoms&& other) :
    m_impl(std::move(other.m_impl)) {
}

Atoms::Atoms(Atoms const& other) :
    m_impl(new Private(*other.m_impl)) {
}

Atoms::Atoms(std::size_t n, AtomInfo const& atom) :
    m_impl(new Private(n, atom)) {
}

Atoms::Atoms(const Atoms& atomsIn,
    const Coords3D& Vx,
    const Coords3D& Vy,
    const Coords3D& Vz,
    const Coords3D& P1) :
    m_impl(new Private) {
    double const VzSqr = Vz.sqr();
    double const VySqr = Vy.sqr();
    double const VxSqr = Vx.sqr();
    for (auto const& atom : atomsIn) {
        auto const V = atom.type.coords - P1;

        double k = (V * Vz) / VzSqr;
        if (k >= 0.0 && k <= 1.0) {
            k = (V * Vy) / VySqr;
            if (k >= 0.0 && k <= 1.0) {
                k = (V * Vx) / VxSqr;
                if (k >= 0.0 && k <= 1.0)
                    push_back(atom); // place atoms to cells
            }
        }
    }
}

Atoms::~Atoms() {
}

Atoms& Atoms::operator=(Atoms&& other) {
    m_impl = std::move(other.m_impl);
    return *this;
}

Atoms& Atoms::operator=(Atoms const& other) {
    m_impl.reset(new Private(*other.m_impl));
    return *this;
}

void Atoms::reserve(std::size_t size) {
    m_impl->reserve(size);
}

void Atoms::push_back(AtomInfo const& atom) {
    m_impl->push_back(atom);
}

Atoms::Iterator Atoms::begin() {
    return m_impl->begin();
}

Atoms::ConstIterator Atoms::begin() const {
    return m_impl->begin();
}

Atoms::Iterator Atoms::end() {
    return m_impl->end();
}

Atoms::ConstIterator Atoms::end() const {
    return m_impl->end();
}

std::size_t Atoms::size() const {
    return m_impl->size();
}

AtomInfo& Atoms::operator[](size_t n) {
    return (*m_impl)[n];
}

void Atoms::clear() {
    m_impl->clear();
}

bool Atoms::empty() const {
    return m_impl->empty();
}

bool Atoms::checkContains(Atoms const& other) {
    return m_impl->checkContains(*other.m_impl);
}

Atoms AtomsHelper::allCellAtoms;

Atoms AtomsHelper::createAllCellAtoms(Coords3DList const& atomTypes) {
    if (!allCellAtoms.empty()) {
        return allCellAtoms;
    }

    auto const cellSize = atomTypes.size();

    allCellAtoms.reserve(9 * 9 * 9 * cellSize);

    // create crystal from cells
    for (int z = 0; z < 9; ++z)
        for (int y = 0; y < 9; ++y)
            for (int x = 0; x < 9; ++x)
                for (size_t a = 0; a < cellSize; ++a) {
                    Coords3D atom = {x + atomTypes[a].x, y + atomTypes[a].y,
                                     z + atomTypes[a].z};
                    allCellAtoms.push_back(AtomInfo(AtomType(atom)));
                }

    return allCellAtoms;
}
