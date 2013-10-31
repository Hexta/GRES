#include "functions.h"
#include <QTime>
#include "QDebug"

float P1 = 1.0;
float P2 = 0.030;
float P3 = 0.00010;

const int NUMBER_OF_ATOMS_IN_CELL = 8;

coords3D atomTypes[] = {
    {1.0, 1.0, 1.0},
    {1.0, 0.5, 0.5},
    {0.5, 1.0, 0.5},
    {0.5, 0.5, 1.0},
    {0.25, 0.25, 0.25},
    {0.75, 0.75, 0.25},
    {0.75, 0.25, 0.75},
    {0.25, 0.75, 0.75}
};

void
addLayer(surface3D &surface, allSoseds &sosedi, int sX, int sY, int sZ) {
    surface2D surfaceXY;
    surfaceXY.reserve(sY);

    for (int y = 0; y < sY; ++y) {
        surface1D surfaceX;
        surfaceX.reserve(sX);

        for (int x = 0; x < sX; ++x) {
            cell cell;
            cell.reserve(sosedi.size());
            for (auto &sosed : sosedi) {
                soseds neighbs;
                char numberNeighbs = 0; //Число первых соседей
                for (int nb = 0; nb < 4; ++nb) {
                    if (x + sosed[nb].x >= 0 && y + sosed[nb].y >= 0 && x + sosed[nb].x < sX && y + sosed[nb].y < sY) {
                        ++numberNeighbs;
                        atomType neighb = {x + sosed[nb].x, y + sosed[nb].y, sZ + sosed[nb].z, sosed[nb].type, false};
                        neighbs.push_back(neighb);
                    }
                }
                atomInfo atom = {neighbs, numberNeighbs, !numberNeighbs};
                cell.push_back(atom);
            }
            surfaceX.push_back(cell);
        }
        surfaceXY.push_back(surfaceX);
    }
    surface.push_back(surfaceXY);
}

atomsCoords
atomsInBox(atomsCoords &atoms, const coords3D &Vx, const coords3D &Vy,
           const coords3D &Vz, const coords3D &P1)
/*
определяет атомы лежащие внутри прямоугольного параллелипипеда
 */ {
    atomsCoords cellAts;

    for (auto &atom : atoms) {
        const float x = atom.x;
        const float y = atom.y;
        const float z = atom.z;
        coords3D X = {x, y, z};
        coords3D V = pointShifting(P1, X);

        double k = ScalarMult(V, Vz) / VectorQuad(Vz);
        if (k >= 0.0 && k <= 1.0) {
            k = ScalarMult(V, Vy) / VectorQuad(Vy);
            if (k >= 0.0 && k <= 1.0) {
                k = ScalarMult(V, Vx) / VectorQuad(Vx);
                if (k >= 0.0 && k <= 1.0)
                    cellAts.push_back(X); //Записываем атомы в ячейке
            }
        }
    }
    return cellAts;
}

void
cellOptim(atomsCoords &ca, float &xs, float &ys, float &zs) {
    atomsCoords newCell;

    for (auto i = ca.begin(); i < ca.end() - 4; ++i) {
        const float x = i->x;
        const float y = i->y;
        const float z = i->z;

        if (x > 0.01 && y > 0.01 && z > 0.01) {
            coords3D atom = {x, y, z};
            newCell.push_back(atom);
        }
    }

    xs = distance(*(ca.end() - 3));
    ys = distance(*(ca.end() - 2));
    zs = distance(*(ca.end() - 1));
    ca = newCell;
}

bool
cells_comp(const atomsCoords &c1, const atomsCoords &c2) {
    return c1.size() < c2.size();
}

bool
compareTranslCell(const atomsCoords &cellAtoms, const atomsCoords &allAtoms,
                  const coords3D &V) {
    bool happy = false;
    for (auto cell_atom_it = cellAtoms.begin();
         cell_atom_it < cellAtoms.end() - 4; ++cell_atom_it) {

        happy = false;
        coords3D At = pointShift(*cell_atom_it, V);

        for (auto &atom : allAtoms) {
            if (coords3Dcompare(At, atom)) {
                happy = true;
                break;
            }
        }
        if (!happy)
            break;
    }

    return happy;
}

bool
coords3Dcompare(const coords3D &coords1, const coords3D &coords2) {
    //сравнение координат
    return coords1.x == coords2.x && coords1.y == coords2.y && coords1.z == coords2.z;
}

void
coordsMove(atomsCoords &ca, const coords3D &O, const coords3D &Vx,
           const coords3D &Vy, const coords3D &Vz) {
    for (auto atom_coords_it = ca.begin(); atom_coords_it < ca.end() - 4; ++atom_coords_it) {
        const double x = ScalarMult(pointShifting(O, *atom_coords_it), Vx) / distance(Vx);
        const double y = ScalarMult(pointShifting(O, *atom_coords_it), Vy) / distance(Vy);
        const double z = ScalarMult(pointShifting(O, *atom_coords_it), Vz) / distance(Vz);

        atom_coords_it->x = x;
        atom_coords_it->y = y;
        atom_coords_it->z = z;
    }
}

double
distance(const double& x1, const double& y1, const double& z1, const double& x2,
         const double& y2, const double& z2) {
    return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2));
}

double
distance(const coords3D &V) {
    return sqrt(pow(V.x, 2) + pow(V.y, 2) + pow(V.z, 2));
}

coords3D
pointShift(const coords3D &A, const coords3D &V) {
    //сдвигает точку A на вектор V
    coords3D result = {A.x + V.x, A.y + V.y, A.z + V.z};
    return result;
}

coords3D
pointShifting(const coords3D &P1, const coords3D &P2) {
    //Разность координат двух точек
    coords3D result = {P2.x - P1.x, P2.y - P1.y, P2.z - P1.z};
    return result;
}

atomsCoords
findCell(int h, int k, int l, float &Xsize, float &Ysize, float &Zsize,
         coords3D &vX, coords3D &vY, coords3D &vZ) {
    const int SIZE_X = 5;
    const int SIZE_Y = 5;
    const int SIZE_Z = 5;
    atomsCoords allAtoms, atomsP1; // атомы, лежащие на плоскостях
    lengthes ls;
    rectangles rectanglesP1; //прямоугольники
    //создадим кристалл из нескольких ячеек
    for (int z = 0; z < 9; ++z)
        for (int y = 0; y < 9; ++y)
            for (int x = 0; x < 9; ++x)
                for (int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a) {
                    coords3D atom = {x + atomTypes[a].x, y + atomTypes[a].y,
                                     z + atomTypes[a].z};
                    allAtoms.push_back(atom);
                }
    //Найдем свободный член в уравнении секущей плоскости hx+ky+lz-C=0
    int C = 0.5 * (h * (SIZE_X - 1) + k * (SIZE_Y - 1) + l * (SIZE_Z - 1)) + 1;
    //Найдем атомы, лежащие на плоскости №1

    atomsP1.reserve(SIZE_Z * SIZE_Y * SIZE_X * NUMBER_OF_ATOMS_IN_CELL);
    for (int z = 0; z < SIZE_Z; ++z)
        for (int y = 0; y < SIZE_Y; ++y)
            for (int x = 0; x < SIZE_X; ++x)
                for (int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a)
                    if (h * (x + atomTypes[a].x) + k * (y + atomTypes[a].y) + l * (z + atomTypes[a].z) - C == 0.0) {
                        coords3D atom = {x + atomTypes[a].x, y + atomTypes[a].y, z + atomTypes[a].z};
                        atomsP1.push_back(atom);
                    }

    //Найдем прямоугольники, лежащие на плоскости №1
    const size_t atomsP1Size = atomsP1.size();
    ls.reserve(0.5 * atomsP1Size * atomsP1Size);
    for (unsigned int i = 0; i < atomsP1Size; ++i)
        for (unsigned int j = i + 1; j < atomsP1Size; ++j) {
            const auto atomP1A = atomsP1[i];
            const auto atomP1B = atomsP1[j];
            float l1 = distance(atomP1A.x, atomP1A.y, atomP1A.z,
                                atomP1B.x, atomP1B.y, atomP1B.z);
            length L = {atomP1A.x, atomP1A.y, atomP1A.z,
                        atomP1B.x, atomP1B.y, atomP1B.z, l1};
            ls.push_back(L);
        }
    float kx = -1;
    float ky = -10;
    float kz = -100;

    for (unsigned int i = 0; i < ls.size(); ++i)
        for (unsigned int j = i + 1; j < ls.size(); ++j) {
            float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, l1, l2;
            x1 = ls[i].x1;
            y1 = ls[i].y1;
            z1 = ls[i].z1;
            x2 = ls[i].x2;
            y2 = ls[i].y2;
            z2 = ls[i].z2;
            x3 = ls[j].x1;
            y3 = ls[j].y1;
            z3 = ls[j].z1;
            x4 = ls[j].x2;
            y4 = ls[j].y2;
            z4 = ls[j].z2;
            l1 = ls[i].length;
            l2 = ls[j].length;
            if (fabs(l1 - l2) < 0.01) {
                if (fabs(x4 - x3) > 0.01)
                    kx = (x2 - x1) / (x4 - x3);
                else if (fabs(x2 - x1) < 0.01)
                    kx = 1;

                if (fabs(y4 - y3) > 0.01)
                    ky = (y2 - y1) / (y4 - y3);
                else if (fabs(y2 - y1) < 0.01)
                    ky = 1;

                if (fabs(z4 - z3) > 0.01)
                    kz = (z2 - z1) / (z4 - z3);
                else if (fabs(z2 - z1) < 0.01)
                    kz = 1;
                if ((kx == ky && kx == kz) &&
                    !((x2 - x1)*(x4 - x2) + (y2 - y1)*(y4 - y2) + (z2 - z1)*(z4 - z2))) {
                    rectangle rect = {x1, y1, z1, x2, y2, z2,
                                      x3, y3, z3, x4, y4, z4, l1};
                    rectanglesP1.push_back(rect);

                    //удаляем дубли
                    for (auto v = ls.begin() + j + 1; v < ls.end(); ++v)
                        if ((v->x1 == x1 && v->x2 == x3 && v->y1 == y1
                             && v->y2 == y3 && v->z1 == z1
                             && v->z2 == z3)
                            || (v->x1 == x2 && v->x2 == x4
                                && v->y1 == y2 && v->y2 == y4
                                && v->z1 == z2 && v->z2 == z4)) {
                            const size_t n = v - ls.begin();
                            ls.erase(v);
                            v = ls.begin() + n - 1;
                        }
                }
            }
        }
    stable_sort(rectanglesP1.begin(), rectanglesP1.end(), rect_comp);

    cells allCells;

    for (auto &rectangle : rectanglesP1) {
        float x1 = rectangle.x1, y1 = rectangle.y1, z1 = rectangle.z1;
        float x2 = rectangle.x2, y2 = rectangle.y2, z2 = rectangle.z2;
        float x3 = rectangle.x3, y3 = rectangle.y3, z3 = rectangle.z3;
        float x4 = rectangle.x4, y4 = rectangle.y4, z4 = rectangle.z4;

        for (float n = 0.5; n < 5; n += 0.5) {
            int atoms = 0;

            for (auto &atom : allAtoms) {
                if ((atom.x == x1 + n * h && atom.y == y1 + n * k && atom.z == z1 + n * l)
                    || (atom.x == x2 + n * h && atom.y == y2 + n * k && atom.z == z2 + n * l)
                    || (atom.x == x3 + n * h && atom.y == y3 + n * k && atom.z == z3 + n * l)
                    || (atom.x == x4 + n * h && atom.y == y4 + n * k && atom.z == z4 + n * l))
                    ++atoms;
            }

            if (atoms == 4) {
                atomsCoords cellAtoms;
                coords3D Vx, Vy, Vz, P1;
                P1.x = x3;
                P1.y = y3;
                P1.z = z3;
                coords3D P2 = {x4, y4, z4}, P3 = {x1, y1, z1};
                coords3D P4 = {x3 + n*h, y3 + n*k, z3 + n * l};

                Vz = pointShifting(P1, P4);
                Vy = pointShifting(P1, P3);
                Vx = pointShifting(P1, P2);
                cellAtoms = atomsInBox(allAtoms, Vx, Vy, Vz, P1);

                // а в конец списка атомов запишем векторы координат и координаты начала координат :)
                cellAtoms.push_back(P1);
                cellAtoms.push_back(Vx);
                cellAtoms.push_back(Vy);
                cellAtoms.push_back(Vz);
                allCells.push_back(cellAtoms);
            }
        }
    }
    stable_sort(allCells.begin(), allCells.end(), cells_comp);

    for (auto &cell : allCells) {
        const size_t cell_size = cell.size();

        //считаем векторы координат
        const coords3D &P1 = *(cell.end() - 4);
        const coords3D &Vx = *(cell.end() - 3);
        const coords3D &Vy = *(cell.end() - 2);
        const coords3D &Vz = *(cell.end() - 1);
        //транслируем по OX
        bool happy = compareTranslCell(cell, allAtoms, Vx);

        //транслируем по OY
        if (happy)
            happy = compareTranslCell(cell, allAtoms, Vy);
        //транслируем по OZ
        if (happy)
            happy = compareTranslCell(cell, allAtoms, Vz);

        //транслируем по -OZ
        if (happy)
            happy = compareTranslCell(cell, allAtoms, -1 * Vz);
        //транслируем по -OY
        if (happy)
            happy = compareTranslCell(cell, allAtoms, -1 * Vy);
        //транслируем по -OX
        if (happy)
            happy = compareTranslCell(cell, allAtoms, -1 * Vx);

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, Vx));
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, Vy));
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, Vz));
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, -1 * Vx));
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, -1 * Vy));
            happy = (cell_size - 4 == translCell.size());
        }

        if (happy) {
            atomsCoords translCell = atomsInBox(allAtoms, Vx, Vy, Vz,
                                                pointShift(P1, -1 * Vz));
            happy = (cell_size - 4 == translCell.size());
        }
        if (happy) {
            allCells.clear();
            allCells.push_back(cell);
            vX = Vx;
            vY = Vy;
            vZ = Vz;
            break;
        }
    }

    auto &firstCell = allCells[0];

    coordsMove(firstCell, *(firstCell.end() - 4), *(firstCell.end() - 3),
               *(firstCell.end() - 2), *(firstCell.end() - 1));

    cellOptim(firstCell, Xsize, Ysize, Zsize);

    return firstCell;
}

void
findSoseds(allSoseds &allSosedi, atomsCoords &atom_Types, float xs, float ys, float zs) {
    const int NUMBER_OF_ATOMS_IN_CELL = atom_Types.size();
    allSosedi.clear();
    for (int type0 = 0; type0 < NUMBER_OF_ATOMS_IN_CELL; ++type0) {
        soseds sosedi;

        double x0 = atom_Types[type0].x;
        double y0 = atom_Types[type0].y;
        double z0 = atom_Types[type0].z;

        for (int xc1 = -1; xc1 < 2; ++xc1)//смещение соседней ячейки №1
            for (int yc1 = -1; yc1 < 2; ++yc1)
                for (int zc1 = -1; zc1 < 2; ++zc1)
                    for (unsigned char type1 = 0;
                         type1 < NUMBER_OF_ATOMS_IN_CELL; ++type1)
                        if (!(type1 == type0 && !xc1 && !yc1 && !zc1)) {
                            double x1 = atom_Types[type1].x + xc1*xs;
                            double y1 = atom_Types[type1].y + yc1*ys;
                            double z1 = atom_Types[type1].z + zc1*zs;
                            if (fabs(pow(x1 - x0, 2) + pow(y1 - y0, 2) + pow(z1 - z0, 2) - 3.0 / 16.0) <= 0.001) //квадрат длины связи #1
                            {
                                atomType s1 = {xc1, yc1, zc1, type1, false};
                                sosedi.push_back(s1);
                            }
                        }
        allSosedi.push_back(sosedi);
    }
}

void
findZmin(const surface3D &surface, int &zm) {
    for (auto z = surface.begin() + zm; z != surface.end(); ++z)
        for (auto y = z->begin() + 2; y != z->end() - 1; ++y) {
            const auto yEnd = y->end();
            for (auto x = y->begin() + 2; x != yEnd - 1; ++x)

                for (auto &atom : *x)
                    if (!atom.deleted) {
                        zm = z - surface.begin();
                        return;
                    }
        }
}

void
recallNeighbours(surface3D &surface, vector<atomType> &surfAtoms, int x, int y,
                 int z, int type) {
    for (auto &neighb : surface[z][y][x][type].neighbours) {
        const int xNb = neighb.x;
        const int yNb = neighb.y;
        const int zNb = neighb.z;
        const unsigned char typeNb = neighb.type;

        atomInfo &neihgbAtomInfo = surface[zNb][yNb][xNb][typeNb];

        size_t i = 0;
        for (auto &neighb_int : neihgbAtomInfo.neighbours) {
            if (neighb_int.x == x && neighb_int.y == y
                && neighb_int.z == z
                && neighb_int.type == type) {
                neihgbAtomInfo.neighbours.erase(neihgbAtomInfo.neighbours.begin() + i);
                neihgbAtomInfo.fNbCount -= 1;
                break;
            }
            ++i;
        }
        if (!neihgbAtomInfo.fNbCount) {
            neihgbAtomInfo.deleted = true;
            neihgbAtomInfo.neighbours.clear();
        }

        if ((neihgbAtomInfo.fNbCount == 3) &&
            ((xNb > 1 && xNb < surface[zNb][yNb].size() - 2)
             && (yNb > 1 && yNb < surface[zNb].size() - 2))) {
            atomType aT = {xNb, yNb, zNb, typeNb, false};
            surfAtoms.push_back(aT);
        }
    }
}

bool
rect_comp(const rectangle &r1, const rectangle &r2) {
    const double S1 = distance(r1.x1, r1.y1, r1.z1, r1.x2, r1.y2, r1.z2) *
            distance(r1.x3, r1.y3, r1.z3, r1.x4, r1.y4, r1.z4);
    const double S2 = distance(r2.x1, r2.y1, r2.z1, r2.x2, r2.y2, r2.z2) *
            distance(r2.x3, r2.y3, r2.z3, r2.x4, r2.y4, r2.z4);
    return S1 < S2;
}

double
ScalarMult(coords3D V1, coords3D V2) {
    //Скалярное произведение векторов
    return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
}

bool
selAtom(surface3D &surface, vector<atomType> &surfAtoms, allSoseds &sosedi,
        int z_min, atomsCoords &tA, vector<bool> &mask, float *rates) {
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = surface[z_min].size();
    const int xC = surface[z_min][0].size();
    bool result = false;
    bool maskON = !mask.empty();

    int i = rand() * surfAtoms.size() / RAND_MAX;

    auto &surfAtom = surfAtoms[i];
    short x = surfAtom.x;
    short y = surfAtom.y;
    unsigned short z = surfAtom.z;
    unsigned short a = surfAtom.type;
    auto &surfaceZY = surface[z][y];
    if ((maskON && ((z == 0
                     && !mask[ (y - 2)*(surfaceZY.size() - 4) + x - 2 ])
                    || (z + tA[a].z >= 0.5)))
        || !maskON) {
        float randN = (float) rand() / RAND_MAX;
        int bonds = surfaceZY[x][a].fNbCount;
        if (!bonds) {
            surfAtoms.erase(surfAtoms.begin() + i);
            --i;
            //continue;
            return false;
        }

        if ((bonds == 1 && randN < P1)
            || (bonds == 2 && randN < P2)
            || (bonds == 3 && randN < P3)) {
            delAtom(surface, surfAtoms, x, y, z, a, i);
            result = true;
            /*
            удалим и соостветсвенный краевой атом
             */
            if (x == 4 || x == 5) {
                delAtom(surface, surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                delAtom(surface, surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                delAtom(surface, surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                delAtom(surface, surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }
            if (z >= surface.size() - 3) {
                addLayer(surface, sosedi, xC, yC, surface.size());
            }
        }
    }

    return result;
}

bool
selAtomCA(surface3D &surface, vector<atomType> &surfAtoms, int z_min,
          atomsCoords &tA, vector<bool> &mask, float* rates) {
    P1 = rates[0];
    P2 = rates[1];
    P3 = rates[2];
    const int yC = surface[z_min].size();
    const int xC = surface[z_min][0].size();
    bool result = false;
    bool maskON = !mask.empty();

    size_t surfAtomsSize = surfAtoms.size();

    for (unsigned int i = 0; i < surfAtomsSize; ++i) {
        const auto &surfAtom = surfAtoms[i];
        short x = surfAtom.x;
        short y = surfAtom.y;
        unsigned short z = surfAtom.z;
        unsigned short a = surfAtom.type;
        const auto &tAz = tA[a].z;
        auto &surfaceZY = surface[z][y];
        if ((maskON && (((z + tAz < 0.5) && !mask[ (y - 2)*(surfaceZY.size() - 4) + x - 2 ])
                        || (z + tAz >= 0.5)))
            || !maskON) {
            float randN = (float) rand() / RAND_MAX;
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
                if (z == surface.size() - 3)
                    result = true;
            }
        }
    }

    size_t i = 0;
    auto surfAtomsEnd = surfAtoms.end();
    for (auto surfAtomIter = surfAtoms.begin(); surfAtomIter != surfAtomsEnd;
         ++surfAtomIter) {
        if (surfAtomIter->toDel) {
            short x = surfAtomIter->x;
            short y = surfAtomIter->y;
            unsigned short z = surfAtomIter->z;
            unsigned short a = surfAtomIter->type;

            delAtom(surface, surfAtoms, x, y, z, a, i);

            if (x == 4 || x == 5) {
                delAtom(surface, surfAtoms, 5 - x, y, z, a, -1);
            }
            if (x == xC - 6 || x == xC - 5) {
                delAtom(surface, surfAtoms, 2 * (xC - 1) - 5 - x, y, z, a, -1);
            }
            if (y == 4 || y == 5) {
                delAtom(surface, surfAtoms, x, 5 - y, z, a, -1);
            }
            if (y == yC - 6 || y == yC - 5) {
                delAtom(surface, surfAtoms, x, 2 * yC - 7 - y, z, a, -1);
            }

            surfAtomIter = surfAtoms.begin() + --i;
            surfAtomsEnd = surfAtoms.end();
        }
        ++i;
    }
    return result;
}

void
delAtom(surface3D &surface, vector<atomType> &surfAtoms, int x, int y, int z,
        int type, int i) {
    if (i != -1) {
        swap(surfAtoms[i], surfAtoms.back());
        surfAtoms.pop_back();
        swap(surfAtoms[i], surfAtoms.back());
        //        surfAtoms.erase(surfAtoms.begin() + i);
    }
    recallNeighbours(surface, surfAtoms, x, y, z, type);
    surface[z][y][x][type].fNbCount = 0;
    surface[z][y][x][type].deleted = true;
    surface[z][y][x][type].neighbours.clear();
}

double
VectorQuad(const coords3D &V) {
    //Возведение вектора в квадрат

    static auto xPrev = V.x;
    static auto yPrev = V.y;
    static auto zPrev = V.z;

    static auto resultPrev = pow(xPrev, 2) + pow(yPrev, 2) + pow(zPrev, 2);
    if (xPrev == V.x && yPrev == V.y && zPrev == V.z)
        return resultPrev;

    xPrev = V.x;
    yPrev = V.y;
    zPrev = V.z;
    resultPrev = pow(xPrev, 2) + pow(yPrev, 2) + pow(zPrev, 2);

    return resultPrev;
}

coords3D operator +(const coords3D& v1, const coords3D& v2) {
    coords3D temp = {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
    return temp;
}

coords3D operator *(const int& n, const coords3D& v) {
    coords3D temp = {n * v.x, n * v.y, n * v.z};
    return temp;
}

void
optimizeSurface(surface3D &surface, int z_min) {
    if (z_min > 1)
        surface[z_min - 2].clear();
    surface.shrink_to_fit();
}

bool operator==(const atomType &a1, const atomType &a2) {
    return ((a1.x == a2.x) && (a1.y == a2.y) && (a1.z == a2.z) && (a1.type == a2.type));
}
