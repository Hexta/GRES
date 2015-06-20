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

#define GL_GLEXT_PROTOTYPES

#define NOMINMAX
#include "Render.h"

#include "functions.h"
#include "SelectAtomMenu.h"
#include "GeometryPrimitives.h"
#include "Surface3D.h"
#include "Cell.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#include <QTime>
#include <QMenu>
#include <QFileDialog>
#include <QWidget>
#include <QMessageBox>
#include <QDebug>
#include <QMouseEvent>
#include <QKeyEvent>

#include <cstdint>

#define BUFFER_OFFSET(i) ((char*)NULL + (i))
#ifdef _WIN32
#define glBindBufferARB pglBindBufferARB
#define glDeleteBuffersARB pglDeleteBuffersARB
#define glBufferDataARB pglBufferDataARB
#define glBufferSubDataARB pglBufferSubDataARB
#endif

namespace {
QTime t;

void sphereTemplate(float R) {
    GLUquadricObj* quadobj;
    quadobj = gluNewQuadric();
    gluSphere(quadobj, R, 10, 10);
    gluDeleteQuadric(quadobj);
    // glutSolidSphere(sR-0.06, 10, 10);
    // glutSolidCube(2*sR-0.04);
    // glutSolidCube(sR-0.04);
}

struct AtomName {
    int name;
    int xC;
    int yC;
    int zC;
    float x;
    float y;
    float z;
    std::uint8_t type;
    int fNbCount;
};

typedef std::vector<AtomName> AtomNames;

struct Bond {
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
};

typedef std::vector<Bond> Bonds;

Coords3D normcrossprod(const Coords3D& in1, const Coords3D& in2) {
    Coords3D out;
    out.x = in1.y * in2.z - in1.z * in2.y;
    out.y = in1.z * in2.x - in1.x * in2.z;
    out.z = in1.x * in2.y - in1.y * in2.x;
    out.normalize();
    return out;
}

} // namespace

struct Render::Private {
    GLfloat drotationX;
    GLfloat drotationY;
    GLfloat drotationZ;
    GLfloat movX;
    GLfloat movY;
    GLfloat sumMovX;
    GLfloat sumMovY;
    GLfloat scale;
    QPoint lastPos;

    SelectAtomMenu* selAtomMenu;

    AtomName selAtomType;

    AtomNames typeAndCoordsOfAtoms;
    AtomNames atNames;

    float sR;

    Coords3D* Vx, * Vy, * Vz;
    GLuint theSphere;
    QColor clearColor;
    Cell cellAtoms;
    Surface3DPtr surfaceXYZ;

    float xs, ys, zs;
    int z_center, z_min;
    int SIZE_X, SIZE_Y;
    bool scribble;
    bool moving;
    bool rotating;
    bool dataChanged;
    bool dataInitialized;
    float scaling;
    Atoms coordsOfAtoms; // список координат атомов
    Bonds bonds;

    QAction* exitAction;

    GRES::VizType visualizationType;
    Cells surfPoints;
    Coords3DList surfVertex;
    Coords3DList surfNormals;
    std::vector<GLuint> buffers;
    int sphereQual;
    int vSize1, vSize2, vSize3;
    std::vector<GLfloat> matrix;
    QGLWidget* parent;

#ifdef _WIN32
    PFNGLBINDBUFFERARBPROC pglBindBufferARB;
    PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB;
    PFNGLBUFFERDATAARBPROC pglBufferDataARB;
    PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB;
#endif

    Private(QGLWidget* parent) :
        drotationX(0),
        drotationY(0),
        drotationZ(0),
        movX(0),
        movY(0),
        sumMovX(0),
        sumMovY(0),
        scale(1.0),
        selAtomMenu(new SelectAtomMenu),
        sR(0.05F),
        surfaceXYZ(new Surface3D),
        xs(1),
        ys(1),
        zs(1),
        z_center(2),
        z_min(0),
        SIZE_X(14),
        SIZE_Y(14),
        scribble(false),
        moving(false),
        rotating(false),
        dataChanged(false),
        dataInitialized(false),
        exitAction(new QAction(QIcon(":/images/exit.png"), tr("E&xit"), parent)),
        visualizationType(GRES::VizType::CELLS_SURFACE),
        sphereQual(3),
        parent(parent)
#ifdef _WIN32
        ,
        pglBindBufferARB(reinterpret_cast<PFNGLBINDBUFFERARBPROC>(wglGetProcAddress("glBindBufferARB"))),
        pglDeleteBuffersARB(reinterpret_cast<PFNGLDELETEBUFFERSARBPROC>(wglGetProcAddress(
                                                                            "glDeleteBuffersARB"))),
        pglBufferDataARB(reinterpret_cast<PFNGLBUFFERDATAARBPROC>(wglGetProcAddress("glBufferDataARB"))),
        pglBufferSubDataARB(reinterpret_cast<PFNGLBUFFERSUBDATAARBPROC>(wglGetProcAddress(
                                                                            "glBufferSubDataARB")))
#endif
    {
        selAtomType.x =
            selAtomType.y =
                selAtomType.z =
                    selAtomType.type =
                        selAtomType.fNbCount = -1;
    }

    void view(Surface3DPtr surface,
        Cell& atTypes,
        float Xsize,
        float Ysize,
        float Zsize,
        int center,
        int min,
        int width,
        int height,
        Coords3D& vX,
        Coords3D& vY,
        Coords3D& vZ,
        GRES::VizType vT) {
        visualizationType = vT;

        xs = Xsize;
        ys = Ysize;
        zs = Zsize;
        surfaceXYZ = surface;
        z_min = min;
        z_center = center;
        cellAtoms = atTypes;
        SIZE_X = width;
        SIZE_Y = height;
        Vx = &vX;
        Vy = &vY;
        Vz = &vZ;
        dataInitialized = true;
        dataChanged = true;

        const float W = 10;
        const float H = 10;

        float h = (SIZE_Y - 4) * ys;
        float w = (SIZE_X - 4) * xs;

        scaling = w > h ? W / w : H / h;

        surfPoints.clear();
        coordsOfAtoms.clear();
        bonds.clear();
        atNames.clear();
        surfVertex.clear();
        surfNormals.clear();
        initMatrix();
        changeVizType(vT);
    }

    void initMatrix() {
        std::vector<GLfloat> matrix1 = {
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0
        };
        matrix = matrix1;
    }

    void changeVizType(GRES::VizType type) {
        visualizationType = type;

        if (!dataInitialized
            || !(dataChanged || coordsOfAtoms.empty()
                 || surfPoints.empty())) {
            parent->updateGL();
            dataChanged = false;
            return;
        }

        atNames.clear();
        bonds.clear();
        surfVertex.clear();
        surfNormals.clear();

        if (!buffers.empty())
            glDeleteBuffersARB(1, &buffers[0]);

        switch (visualizationType) {
            case GRES::VizType::CELLS_SURFACE:
                if (!surfPoints.empty()) {
                    break;
                }

                createSurfacePoints();
                buffers.clear();
                buffers.push_back(1);
                break;
            case GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK: {
                if (!coordsOfAtoms.empty()) {
                    break;
                }

                createAtomsAndBonds();
                GeometryPrimitives::createSphere(0.09 * scaling, 10, 10, vSize1, vSize2,
                    vSize3);
                buffers.clear();
                buffers.push_back(1);
                buffers.push_back(2);
                buffers.push_back(3);
            }
            break;
            case GRES::VizType::ATOMS_AND_BONDS_SURFACE: {
                if (!coordsOfAtoms.empty()) {
                    break;
                }

                createAtomsAndBondes();
                GeometryPrimitives::createSphere(0.09 * scaling, 10, 10, vSize1, vSize2,
                    vSize3);
                buffers.clear();
                buffers.push_back(1);
                buffers.push_back(2);
                buffers.push_back(3);
            }
            break;
            case GRES::VizType::ATOMS_SURFACE_AND_BULK: {
                if (!coordsOfAtoms.empty()) {
                    break;
                }

                createAtomsAndBonds();
                GeometryPrimitives::createSphere(0.2 * scaling, 10, 10, vSize1, vSize2,
                    vSize3);
                buffers.clear();
                buffers.push_back(1);
                buffers.push_back(2);
                buffers.push_back(3);
            }
            break;
            case GRES::VizType::ATOMS_SURFACE: {
                if (!coordsOfAtoms.empty()) {
                    break;
                }

                createAtomsAndBondes();
                GeometryPrimitives::createSphere(0.2 * scaling, 10, 10, vSize1, vSize2,
                    vSize3);
                buffers.clear();
                buffers.push_back(1);
                buffers.push_back(2);
                buffers.push_back(3);
            }
            break;
            default:
                break;
        }

        parent->updateGL();
        dataChanged = false;
    }

    void processSelection(int x, int y) {
#define BUFFER_LENGTH 512
        GLfloat fAspect;
        static GLuint selectBuff[BUFFER_LENGTH];
        GLint hits, viewport[4];
        glSelectBuffer(BUFFER_LENGTH, selectBuff);
        glGetIntegerv(GL_VIEWPORT, viewport);
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glRenderMode(GL_SELECT);
        glLoadIdentity();

        gluPickMatrix(x, viewport[3] - y, 1, 1, viewport);
        fAspect = static_cast<float>(viewport[2]) / static_cast<float>(viewport[3]);
        gluPerspective(40.0f, fAspect, 1.0, 80.0);
        paintGL();
        hits = glRenderMode(GL_RENDER);

        if (hits == 1)
            processAtom(selectBuff);
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        parent->updateGL();
    }

    void paintGL() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        draw();
    }

    void createSurfacePoints() {
        auto& surface = *surfaceXYZ;

        const size_t dX = surface[z_min][0].size() - 3;
        const size_t dY = surface[z_min].size() - 3;

        AtomInfo stubAtom(AtomType(Coords3D(-1.0, -1.0, -1.0)));
        Cells points(dY, Cell(Atoms(dX, stubAtom)));

        for (size_t z = z_min; z < surface.size() - 2; ++z)
            for (size_t y = surface[z].size() - 2; --y >= 2;)
                for (size_t x = surface[z][y].size() - 2; --x >= 2;)
                    for (const auto& atom : surface[z][y][x].getAtoms())
                        if (!atom.deleted) {
                            auto& atom_ = points[y - 2].getAtoms()[x - 2];
                            if (cmp_float(atom_.type.coords.x, -1.0)) {
                                atom_.type.coords = Coords3D(scaling * (x - 2) * xs,
                                    scaling * (y - 2) * ys, -scaling * z * zs);
                            }
                        }

        for (size_t i = 0; i < dY; ++i) {
            auto const& coords = points[i].getAtoms()[dX - 2].type.coords;

            points[i].getAtoms()[dX - 1].type.coords = Coords3D(
                coords.x + scaling * xs,
                coords.y,
                coords.z);
        }

        for (size_t i = 0; i < dX; ++i) {
			auto const& coords = points[dY - 2].getAtoms()[i].type.coords;

            points[dY - 1].getAtoms()[i].type.coords = Coords3D(
                coords.x,
                coords.y + scaling * ys,
                coords.z);
        }
        surfVertex.reserve((SIZE_Y - 4) * (SIZE_X - 4));
        surfNormals.reserve((SIZE_Y - 4) * (SIZE_X - 4));

        for (int i = 0; i < SIZE_Y - 4; i++)
            for (int j = 0; j < SIZE_X - 4; j++) {
                Cell& cell = points[i];
                Cell& cellNext = points[i + 1];

                auto const& coordsCellJNext = cell.getAtoms()[j + 1].type.coords;
                auto const& coordsCellNextJNext = cellNext.getAtoms()[j + 1].type.coords;

                auto const& coordsCellJ = cell.getAtoms()[j].type.coords;
                auto const& coordsCellNextJ = cellNext.getAtoms()[j].type.coords;

                Coords3D const& v1 = coordsCellJ;
                Coords3D const v2 = {coordsCellJNext.x, coordsCellJNext.y, coordsCellJ.z};
                Coords3D const v3 = {coordsCellNextJNext.x, coordsCellNextJNext.y, coordsCellJ.z};
                Coords3D const& v4 = coordsCellNextJ;

                if (!cmp_float(coordsCellJ.z, coordsCellJNext.z)) {
                    Coords3D const v5 = {coordsCellJNext.x, coordsCellJNext.y, coordsCellJ.z};
                    Coords3D const& v6 = coordsCellJNext;
                    Coords3D const& v7 = coordsCellNextJNext;
                    Coords3D const v8 = {coordsCellNextJNext.x, coordsCellNextJNext.y, coordsCellJ.z};

                    surfNormals.insert(surfNormals.end(), 4, normcrossprod(v6 + -1 * v5, v7 + -1 * v5));

                    surfVertex.push_back(v5);
                    surfVertex.push_back(v6);
                    surfVertex.push_back(v7);
                    surfVertex.push_back(v8);
                }

                if (!cmp_float(coordsCellJ.z, coordsCellNextJ.z)) {
                    Coords3D const v9 = {coordsCellNextJNext.x, coordsCellNextJNext.y, coordsCellJ.z};
                    Coords3D const v10 = {coordsCellNextJNext.x, coordsCellNextJNext.y, coordsCellNextJ.z};
                    Coords3D const& v11 = coordsCellNextJ;
                    Coords3D const v12 = {coordsCellNextJ.x, coordsCellNextJ.y, coordsCellJ.z};

                    surfNormals.insert(surfNormals.end(), 4, normcrossprod(v10 + -1 * v9, v11 + -1 * v9));

                    surfVertex.push_back(v9);
                    surfVertex.push_back(v10);
                    surfVertex.push_back(v11);
                    surfVertex.push_back(v12);
                }

                surfVertex.push_back(v1);
                surfVertex.push_back(v2);
                surfVertex.push_back(v3);
                surfVertex.push_back(v4);
                surfNormals.insert(surfNormals.end(), 4, normcrossprod(v2 + -1 * v1, v3 + -1 * v1));
            }

        glBindBufferARB(GL_ARRAY_BUFFER, 1);
        glBufferDataARB(GL_ARRAY_BUFFER, (surfVertex.size() + surfNormals.size()) * 3 * sizeof(float),
            0, GL_STATIC_DRAW);
        glBufferSubDataARB(GL_ARRAY_BUFFER, 0, surfVertex.size() * 3 * sizeof(float),
            &surfVertex[0].x);
        glBufferSubDataARB(GL_ARRAY_BUFFER, surfVertex.size() * 3 * sizeof(float),
            surfNormals.size() * 3 * sizeof(float), &surfNormals[0].x);
        glBindBufferARB(GL_ARRAY_BUFFER, 0);
    }

    void createAtomsAndBonds() {
        auto& cellAts = cellAtoms.getAtoms();
        auto& surface = *surfaceXYZ;

        int name = 0;
        for (size_t z = z_min; z < surface.size() - 2; ++z)
            for (size_t y = surface[z].size() - 2; --y >= 2;)
                for (size_t x = surface[z][y].size() - 2; --x >= 2;) {
                    float x0 = scaling * x * xs;
                    float y0 = scaling * y * ys;
                    float z0 = -scaling * (z - z_min) * zs;
                    const char atomsCount = static_cast<char>(surface[z][y][x].size());
                    for (unsigned char a = atomsCount; --a > 0;) {
                        if (surface[z][y][x].getAtoms()[a].fNbCount) {
                            float xA = x0 + scaling * cellAts[a].type.coords.x;
                            float yA = y0 + scaling * cellAts[a].type.coords.y;
                            float zA = z0 - scaling * cellAts[a].type.coords.z;

                            ++name;
                            AtomName temp = {name, static_cast<int>(x),
                                             static_cast<int>(y),
                                             static_cast<int>(z), xA, yA, zA, a,
                                             static_cast<int>(atomsCount)};
                            atNames.push_back(temp);

                            for (auto& nb : surface[z][y][x].getAtoms()[a].neighbors) {
                                float xNb = scaling * (xs * nb.x + cellAts[nb.type].type.coords.x);
                                float yNb = scaling * (ys * nb.y + cellAts[nb.type].type.coords.y);
                                float zNb = scaling * (-zs * nb.z - cellAts[nb.type].type.coords.z);

                                Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                                bonds.push_back(bondT);
                            }
                        }
                    }
                }
    }

    void createAtomsAndBondes() {
        auto& atoms = cellAtoms.getAtoms();

        int name = 0;

        auto& surface = *surfaceXYZ;
        for (auto const& surfAtom : surface.getSurfaceAtoms()) {
            const int x = surfAtom.x;
            const int y = surfAtom.y;
            const int z = surfAtom.z;
            const unsigned char a = surfAtom.type;

            auto& neighbors = surface[z][y][x].getAtoms()[a].neighbors;

            if (!neighbors.empty()) {
                float x0 = scaling * x * xs;
                float y0 = scaling * y * ys;
                float z0 = -scaling * (z - z_min) * zs;

                float xA = x0 + scaling * atoms[a].type.coords.x;
                float yA = y0 + scaling * atoms[a].type.coords.y;
                float zA = z0 - scaling * atoms[a].type.coords.z;

                ++name;
                AtomName temp = {name, x, y, z, xA, yA, zA, a, static_cast<int>(neighbors.size())};
                atNames.push_back(temp);

                for (auto& nb : neighbors) {
                    float xNb = scaling * (xs * nb.x + atoms[nb.type].type.coords.x);
                    float yNb = scaling * (ys * nb.y + atoms[nb.type].type.coords.y);
                    float zNb = scaling * (-zs * nb.z - atoms[nb.type].type.coords.z);

                    Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                    bonds.push_back(bondT);
                }
            }
        }

        for (auto const& surfAtom : surface.getSurfaceAtoms()) {
            int const x_ = surfAtom.x;
            int const y_ = surfAtom.y;
            int const z_ = surfAtom.z;
            unsigned char a_ = surfAtom.type;

            for (auto& nb : surface[z_][y_][x_].getAtoms()[a_].neighbors) {
                int x = nb.x;
                int y = nb.y;
                int z = nb.z;
                unsigned char a = nb.type;

                auto& surfaceZYXA = surface[z][y][x].getAtoms()[a];

                if (x > 1 && x < static_cast<decltype(x)>(surface[z_][y_].size()) - 2
                    && y > 1 && y < static_cast<decltype(y)>(surface[z_].size()) - 2)
                    if (!surfaceZYXA.deleted) {
                        float x0 = scaling * x * xs;
                        float y0 = scaling * y * ys;
                        float z0 = -scaling * (z - z_min) * zs;

                        float xA = x0 + scaling * atoms[a].type.coords.x;
                        float yA = y0 + scaling * atoms[a].type.coords.y;
                        float zA = z0 - scaling * atoms[a].type.coords.z;

                        ++name;
                        AtomName temp = {name, x, y, z, xA, yA, zA, a,
                                         static_cast<int>(surfaceZYXA.neighbors.size())};
                        atNames.push_back(temp);

                        for (auto& nb_int : surfaceZYXA.neighbors) {
                            float xNb = scaling * (xs * nb_int.x + atoms[nb_int.type].type.coords.x);
                            float yNb = scaling * (ys * nb_int.y + atoms[nb_int.type].type.coords.y);
                            float zNb = scaling * (-zs * nb_int.z - atoms[nb_int.type].type.coords.z);

                            Bond bondT = {xA, yA, zA,
                                          xNb, yNb, zNb};
                            bonds.push_back(bondT);
                        }
                    }
            }
        }
    }

    void draw() {
        if (dataInitialized) {
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();

            if (visualizationType == GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK
                || visualizationType == GRES::VizType::ATOMS_SURFACE_AND_BULK
                || visualizationType == GRES::VizType::ATOMS_AND_BONDS_SURFACE
                || visualizationType == GRES::VizType::ATOMS_SURFACE) {
                /*GLfloat ambientLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
                  GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };*/
                GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
                GLfloat diffuseLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
                GLfloat specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
                GLfloat light_position[] = {10.0, 10.0, 20.0, 0.0f};
                // GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
                // GLfloat mat_shininess[] = { 200.0 };
                // glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
                // glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

                // GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
                // GLfloat mat_shininess[] = { 100.0 };
                // glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
                // glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
                GLfloat mat_specular[] = {0.3f, 0.3f, 0.3f, 1.0f};
                GLfloat mat_shininess[] = {100.0};
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
                glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
                glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
                glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
                glLightfv(GL_LIGHT0, GL_POSITION, light_position);
                glEnable(GL_LIGHTING);
                glEnable(GL_LIGHT0);
                glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
                glEnable(GL_COLOR_MATERIAL);
                glInitNames();
                glPushName(0);

                setGeometry(-20.0 - zs * scaling * (z_center - z_min));

                glTranslatef(-6.0, -6.0, -20.0);
                drawAxis();
                glTranslatef(6.0, 6.0, 20.0);

                glTranslatef(-5.0 - 2.0 * xs * scaling,
                    -5.0 - 2.0 * ys * scaling, -20.0);
                glColor3f(0.98f, 0.625f, 0.12f);
                // sphere( 5,10,10 );

                if (!(visualizationType == GRES::VizType::ATOMS_SURFACE_AND_BULK
                      || visualizationType == GRES::VizType::ATOMS_SURFACE)
                    || rotating || scribble || moving) {
                    glColor3f(0, 0, 1.0);
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glVertexPointer(3, GL_FLOAT, 0, &bonds[0].x1);
                    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2 * bonds.size()));
                    glDisableClientState(GL_VERTEX_ARRAY);
                }

                if (!rotating && !scribble && !moving) {
                    int name = 0;
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glEnableClientState(GL_NORMAL_ARRAY);
                    // t.start();
                    for (auto const& atomInfo : atNames) {
                        glColor3f(105.0f / 255, 68.0f / 255, 9.0f / 255);
                        if (atomInfo.fNbCount == 3)
                            glColor3f(0.08f, 0.925f, 0.02f);

                        if (atomInfo.fNbCount == 2)
                            glColor3f(1.0, 1.0, 0.0);

                        if (atomInfo.fNbCount == 1)
                            glColor3f(1.0, 0.0, 0.0);

                        if (cmp_float(atomInfo.x, selAtomType.x)
                            && cmp_float(atomInfo.y, selAtomType.y)
                            && cmp_float(atomInfo.z, selAtomType.z))
                            glColor3f(1.0, 0.0, 1.0);

                        glPushMatrix();

                        glTranslatef(atomInfo.x, atomInfo.y, atomInfo.z);
                        ++name;
                        glLoadName(name);

                        glBindBufferARB(GL_ARRAY_BUFFER, 1);
                        glVertexPointer(3, GL_FLOAT, 0, 0);
                        glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof(float) * vSize1));

                        glDrawArrays(GL_TRIANGLE_FAN, 0, vSize1);

                        glBindBufferARB(GL_ARRAY_BUFFER, 2);
                        glVertexPointer(3, GL_FLOAT, 0, 0);
                        glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof(float) * vSize2));

                        glDrawArrays(GL_TRIANGLE_FAN, 0, vSize2);

                        glBindBufferARB(GL_ARRAY_BUFFER, 3);
                        glVertexPointer(3, GL_FLOAT, 0, 0);
                        glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof(float) * vSize3));

                        glDrawArrays(GL_QUAD_STRIP, 0, vSize3);

                        glPopMatrix();
                    }

                    // qDebug("Time elapsed: %d ms", t.elapsed());
                    glDisableClientState(GL_VERTEX_ARRAY);
                    glDisableClientState(GL_NORMAL_ARRAY);
                    glBindBufferARB(GL_ARRAY_BUFFER, 0);
                }

                glFlush();
            } else if (visualizationType == GRES::VizType::CELLS_SURFACE) {
                GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
                GLfloat diffuseLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
                GLfloat specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
                GLfloat light_position[] = {5.0f, 5.0f, 5.0f, 0.0f};
                // GLfloat light_position[] = { 10.0, 10.0, 20.0, 0.0f };
                GLfloat mat_specular[] = {0.3f, 0.3f, 0.3f, 1.0f};
                GLfloat mat_shininess[] = {100.0};
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
                glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

                glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
                glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
                glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
                glLightfv(GL_LIGHT0, GL_POSITION, light_position);
                glEnable(GL_LIGHTING);
                glEnable(GL_LIGHT0);

                glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
                glEnable(GL_COLOR_MATERIAL);

                setGeometry(-20.0);

                glTranslatef(-6.0, -6.0, -20.0);
                drawAxis();
                glTranslatef(6.0, 6.0, 20.0);

                glColor3f(0.18f, 0.325f, 0.82f);

                glTranslatef(-5.0, -5.0, -20.0 + zs * z_min * scaling);

                glPointSize(10.0f);
                glDisable(GL_CULL_FACE);

                glBindBufferARB(GL_ARRAY_BUFFER, 1);
                glVertexPointer(3, GL_FLOAT, 0, 0);
                glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof(float) * surfVertex.size()));

                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glDrawArrays(GL_QUADS, 0, static_cast<GLsizei>(surfVertex.size()));
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_NORMAL_ARRAY);
                glBindBufferARB(GL_ARRAY_BUFFER, 0);

                /*
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, &surfVertex[0].x);
                glNormalPointer(GL_FLOAT, 0, &surfNormals[0].x);
                glDrawArrays(GL_QUADS, 0, surfVertex.size());
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_NORMAL_ARRAY);*/
                glFlush();
                glEnable(GL_CULL_FACE);
            }
        }
        // qDebug("Time elapsed: %d ms", t.elapsed());
    }

    void setGeometry(GLfloat zCenter) {
        glTranslatef(0, 0, zCenter);

        if (!cmp_float(drotationX, 0))
            glRotatef(drotationX, 1.0, 0.0, 0.0);

        if (!cmp_float(drotationY, 0))
            glRotatef(drotationY, 0.0, 1.0, 0.0);

        glTranslatef(0, 0, -zCenter);

        drotationX = drotationY = 0.0;
        glMultMatrixf(&matrix[0]);
        glGetFloatv(GL_MODELVIEW_MATRIX, &matrix[0]);

        glLoadIdentity();

        glTranslatef(0, 0, zCenter);

        if (!cmp_float(scale, 1.0))
            glScalef(scale, scale, scale);
        glTranslatef(0, 0, -zCenter);
        glTranslatef(sumMovX, -sumMovY, 0.0);
        glMultMatrixf(&matrix[0]);
    }

    void drawAxis() {
        glLineWidth(7.0);
        glColor3f(1.0, 1.0, 0.0);
        glLineWidth(5.0);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glBegin(GL_LINES);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(4.0, 0.0, 0.0);

        glColor3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, -4.0);

        glColor3f(1.0, 1.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 4.0, 0.0);
        glEnd();
        glDisable(GL_LINE_SMOOTH);
        glDisable(GL_POLYGON_SMOOTH);
        glDisable(GL_BLEND);

        glColor3f(1.0, 1.0, 0.0);
        QString xText = "x (" + QString::number(Vx->x) +
                        "," + QString::number(Vx->y) + "," + QString::number(Vx->z) + ")";
        QString yText = "y (" + QString::number(Vy->x) +
                        "," + QString::number(Vy->y) + "," + QString::number(Vy->z) + ")";
        QString zText = "z (" + QString::number(Vz->x) +
                        "," + QString::number(Vz->y) + "," + QString::number(Vz->z) + ")";
        parent->renderText(4.2, -0.1, 0.0, xText);
        glColor3f(1.0, 0.0, 1.0);
        parent->renderText(0.0, -0.1, -4.3, zText);
        glColor3f(1.0, 1.0, 1.0);
        parent->renderText(0.0, 4.3, 0.0, yText);
    }

    void initializeGL() {
#ifdef _WIN32
        pglBindBufferARB = reinterpret_cast<PFNGLBINDBUFFERARBPROC>(wglGetProcAddress("glBindBufferARB"));
        pglDeleteBuffersARB =
            reinterpret_cast<PFNGLDELETEBUFFERSARBPROC>(wglGetProcAddress("glDeleteBuffersARB"));
        pglBufferDataARB = reinterpret_cast<PFNGLBUFFERDATAARBPROC>(wglGetProcAddress("glBufferDataARB"));
        pglBufferSubDataARB =
            reinterpret_cast<PFNGLBUFFERSUBDATAARBPROC>(wglGetProcAddress("glBufferSubDataARB"));
#endif
    }

    void processAtom(const GLuint* pSelectBuff) {
        int id;
        id = pSelectBuff[3];
        for (auto const& atName : atNames) {
            if (atName.name == id) {
                selAtomType = atName;
                break;
            }
        }
    }

    void processSelectionMenu() {
        selAtomMenu->setWindowFlags(Qt::Tool);
        selAtomMenu->setWindowTitle("Atom`s information");
        selAtomMenu->move(QCursor::pos());
        selAtomMenu->setInfo(selAtomType.xC, selAtomType.yC,
            selAtomType.zC, selAtomType.type, selAtomType.fNbCount);
        selAtomMenu->show();
    }

    void mouseMoveEvent(QMouseEvent* event) {
        GLfloat dx = static_cast<GLfloat>(event->x() - lastPos.x()) / parent->width();
        GLfloat dy = static_cast<GLfloat>(event->y() - lastPos.y()) / parent->height();
        if (event->buttons() & Qt::LeftButton) {
            drotationX = 180 * dy;
            drotationY = 180 * dx;
            rotating = true;
            parent->updateGL();
        } else if (event->buttons() & Qt::RightButton) {
            scribble = true;
            scale += dy * scale;
            parent->updateGL();
        } else if (event->buttons() & Qt::MidButton) {
            movX = 18 * dx;
            movY = 18 * dy;
            sumMovX += movX;
            sumMovY += movY;
            parent->updateGL();
        }

        lastPos = event->pos();
    }
};

Render::Render(QWidget* parent) : QGLWidget(parent),
                                  d(new Private(this)) {
}

void Render::changeVizType(GRES::VizType type) {
    d->changeVizType(type);
}

void Render::initializeGL() {
    d->initializeGL();

    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    d->theSphere = glGenLists(1);
    glNewList(d->theSphere, GL_COMPILE);
    sphereTemplate(d->sR * 10);
    glEndList();
}

void Render::resizeGL(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40.0,
        static_cast<GLfloat>(width) / static_cast<GLfloat>(height),
        1.0, 80.0);
    glMatrixMode(GL_MODELVIEW);
}

void Render::mousePressEvent(QMouseEvent* event) {
    d->lastPos = event->pos();
    if (event->buttons() & Qt::LeftButton) {
        if (event->modifiers() == Qt::ControlModifier)
            d->processSelection(event->x(), event->y());
        else
            d->rotating = true;
    } else if (event->buttons() & Qt::MidButton) {
        d->moving = true;
    }
}

void Render::mouseReleaseEvent(QMouseEvent* event) {
    if ((!d->scribble) && (event->button() == Qt::RightButton)) {
        d->processSelectionMenu();
    }
    if ((d->scribble) && (event->button() == Qt::RightButton)) {
        d->scribble = false;
        updateGL();
    }
    if (d->rotating && event->button() == Qt::LeftButton) {
        d->rotating = false;
        updateGL();
    }
    if (d->moving && event->button() == Qt::MidButton) {
        d->moving = false;
        updateGL();
    }
}

void Render::mouseMoveEvent(QMouseEvent* event) {
    d->mouseMoveEvent(event);
}

void Render::keyPressEvent(QKeyEvent* event) {
    if (event->text() == "p")
        emit(etching());
}

void Render::saveResult() {
    QPixmap outPixmap;
    this->raise();
    outPixmap = QPixmap::grabWindow(this->winId());

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
        "./untitled.png",
        tr("Images (*.png *.jpg *.tiff *.bmp *.xpm)"));
    if (!fileName.isNull()) { // check that user don't canceled file selection
        if (!outPixmap.save(fileName))
            QMessageBox::information(this, tr("GRES"), tr("Cannot save %1.").arg(fileName));
    }
}

QSize Render::minimumSizeHint() const {
    return QSize(400, 400);
}

QSize Render::sizeHint() const {
    return QSize(600, 600);
}

void Render::view(Surface3DPtr surface,
    Cell& atTypes,
    float Xsize,
    float Ysize,
    float Zsize,
    int center,
    int min,
    int width,
    int height,
    Coords3D& vX,
    Coords3D& vY,
    Coords3D& vZ,
    GRES::VizType vT) {
    d->view(surface, atTypes, Xsize, Ysize, Zsize,
        center, min, width, height, vX, vY, vZ, vT);
}

Render::~Render() {
}

void Render::paintGL() {
    d->paintGL();
}
