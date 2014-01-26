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

#define GL_GLEXT_PROTOTYPES

#define NOMINMAX 
#include "Render.h"
#include "selectAtomMenu.h"
#include "geometry.h"

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

#define BUFFER_OFFSET(i) ((char*)NULL + (i))
#ifdef _WIN32
#define glBindBufferARB           d->pglBindBufferARB
#define glDeleteBuffersARB        d->pglDeleteBuffersARB
#define glBufferDataARB           d->pglBufferDataARB
#define glBufferSubDataARB        d->pglBufferSubDataARB
#endif

QTime t;

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

    atomName selAtomType;

    AtomsNames typeAndCoordsOfAtoms;
    AtomsNames atNames;

    float sR;

    Coords3D *Vx, *Vy, *Vz;
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
    Atoms coordsOfAtoms; //список координат атомов
    Bonds bonds;

    QAction *exitAction;

    GRES::VizType visualizationType;
    Cells surfPoints;
    Atoms surfVertex;
    Atoms surfNormals;
    std::vector <GLuint> buffers;
    int sphereQual;
    int vSize1, vSize2, vSize3;
    std::vector<GLfloat> matrix;

#ifdef _WIN32
    PFNGLBINDBUFFERARBPROC pglBindBufferARB;
    PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB;
    PFNGLBUFFERDATAARBPROC pglBufferDataARB;
    PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB;
#endif

    Private(QObject* parent) :
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
        sphereQual(3)
#ifdef _WIN32
      , pglBindBufferARB(reinterpret_cast<PFNGLBINDBUFFERARBPROC>(wglGetProcAddress("glBindBufferARB"))),
        pglDeleteBuffersARB(reinterpret_cast<PFNGLDELETEBUFFERSARBPROC>(wglGetProcAddress("glDeleteBuffersARB"))),
        pglBufferDataARB(reinterpret_cast<PFNGLBUFFERDATAARBPROC>(wglGetProcAddress("glBufferDataARB"))),
        pglBufferSubDataARB(reinterpret_cast<PFNGLBUFFERSUBDATAARBPROC>(wglGetProcAddress("glBufferSubDataARB")))
#endif
    {
        selAtomType.x =
            selAtomType.y =
            selAtomType.z =
            selAtomType.type =
            selAtomType.fNbCount = -1;
    }
};

Render::Render(QWidget *parent) : QGLWidget(parent), d(new Private(this)) {
    createActions();
}

void
Render::changeVizType(GRES::VizType type) {
    d->visualizationType = type;

    if (!d->dataInitialized
        || !(d->dataChanged || d->coordsOfAtoms.empty()
        || d->surfPoints.empty())) {
        updateGL();
        d->dataChanged = false;
        return;
    }

    d->atNames.clear();
    d->bonds.clear();
    d->surfVertex.clear();
    d->surfNormals.clear();

    if (!d->buffers.empty())
        glDeleteBuffersARB(1, &d->buffers[0]);

    switch (d->visualizationType) {
        case GRES::VizType::CELLS_SURFACE:
            if (!d->surfPoints.empty()) {
                break;
            }

            createSurfacePoints();
            d->buffers.clear();
            d->buffers.push_back(1);
            break;
        case GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK:
        {
            if (!d->coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBonds();
            createSphere(0.09 * d->scaling, 10, 10, d->vSize1, d->vSize2,
                d->vSize3);
            d->buffers.clear();
            d->buffers.push_back(1);
            d->buffers.push_back(2);
            d->buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_AND_BONDS_SURFACE:
        {
            if (!d->coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBondes();
            createSphere(0.09 * d->scaling, 10, 10, d->vSize1, d->vSize2,
                d->vSize3);
            d->buffers.clear();
            d->buffers.push_back(1);
            d->buffers.push_back(2);
            d->buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_SURFACE_AND_BULK:
        {
            if (!d->coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBonds();
            createSphere(0.2 * d->scaling, 10, 10, d->vSize1, d->vSize2,
                d->vSize3);
            d->buffers.clear();
            d->buffers.push_back(1);
            d->buffers.push_back(2);
            d->buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_SURFACE:
        {
            if (!d->coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBondes();
            createSphere(0.2 * d->scaling, 10, 10, d->vSize1, d->vSize2,
                d->vSize3);
            d->buffers.clear();
            d->buffers.push_back(1);
            d->buffers.push_back(2);
            d->buffers.push_back(3);
        }
            break;
        default: break;
    }

    updateGL();
    d->dataChanged = false;
}

void
Render::createActions() {
}

void Render::createAtomsAndBonds()
{
    auto& cellAts = d->cellAtoms.atoms;
    auto& surface = *d->surfaceXYZ;

    int name = 0;
    for (size_t z = d->z_min; z < surface.size() - 2; ++z)
        for (size_t y = surface[z].size() - 2; --y >= 2;)
            for (size_t x = surface[z][y].size() - 2; --x >= 2;) {
                float x0 = d->scaling * x * d->xs;
                float y0 = d->scaling * y * d->ys;
                float z0 = -d->scaling * (z - d->z_min) * d->zs;
                const char atomsCount = static_cast<char>(surface[z][y][x].size());
                for (unsigned char a = atomsCount; --a > 0;) {

                    if (surface[z][y][x][a].fNbCount) {
                        float xA = x0 + d->scaling * cellAts[a].x;
                        float yA = y0 + d->scaling * cellAts[a].y;
                        float zA = z0 - d->scaling * cellAts[a].z;

                        ++name;
                        atomName temp = {name, static_cast<int> (x),
                                         static_cast<int> (y),
                                         static_cast<int> (z), xA, yA, zA, a,
                                         static_cast<int> (atomsCount)};
                        d->atNames.push_back(temp);

                        for (auto &nb : surface[z][y][x][a].neighbors) {
                            float xNb = d->scaling * (d->xs * nb.x + cellAts[nb.type].x);
                            float yNb = d->scaling * (d->ys * nb.y + cellAts[nb.type].y);
                            float zNb = d->scaling * (-d->zs * nb.z - cellAts[nb.type].z);

                            Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                            d->bonds.push_back(bondT);
                        }
                    }
                }
            }
}

void Render::createAtomsAndBondes() {
    auto& cellAtoms = d->cellAtoms.atoms;

    int name = 0;
    auto const& scaling = d->scaling;
    auto& surface = *d->surfaceXYZ;
    for (auto const& surfAtom : surface.getSurfaceAtoms()) {
        const int x = surfAtom.x;
        const int y = surfAtom.y;
        const int z = surfAtom.z;
        const unsigned char a = surfAtom.type;

        auto &neighbors = surface[z][y][x][a].neighbors;

        if (!neighbors.empty()) {
            float x0 = scaling * x * d->xs;
            float y0 = scaling * y * d->ys;
            float z0 = -scaling * (z - d->z_min) * d->zs;

            float xA = x0 + scaling * cellAtoms[a].x;
            float yA = y0 + scaling * cellAtoms[a].y;
            float zA = z0 - scaling * cellAtoms[a].z;

            ++name;
            atomName temp = {name, x, y, z, xA, yA, zA, a, static_cast<int> (neighbors.size())};
            d->atNames.push_back(temp);

            for (auto &nb : neighbors) {
                float xNb = scaling * (d->xs * nb.x + cellAtoms[nb.type].x);
                float yNb = scaling * (d->ys * nb.y + cellAtoms[nb.type].y);
                float zNb = scaling * (-d->zs * nb.z - cellAtoms[nb.type].z);

                Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                d->bonds.push_back(bondT);
            }
        }
    }

    for (auto const& surfAtom : surface.getSurfaceAtoms()) {
        int const x_ = surfAtom.x;
        int const y_ = surfAtom.y;
        int const z_ = surfAtom.z;
        unsigned char a_ = surfAtom.type;

        for (auto &nb : surface[z_][y_][x_][a_].neighbors) {
            int x = nb.x;
            int y = nb.y;
            int z = nb.z;
            unsigned char a = nb.type;

            auto &surfaceZYXA = surface[z][y][x][a];

            if (x > 1 && x < static_cast<decltype(x)>(surface[z_][y_].size()) - 2
                && y > 1 && y < static_cast<decltype(y)>(surface[z_].size()) - 2)
            if (!surfaceZYXA.deleted) {

                float x0 = scaling * x * d->xs;
                float y0 = scaling * y * d->ys;
                float z0 = -scaling * (z - d->z_min) * d->zs;

                float xA = x0 + scaling * cellAtoms[a].x;
                float yA = y0 + scaling * cellAtoms[a].y;
                float zA = z0 - scaling * cellAtoms[a].z;

                ++name;
                atomName temp = {name, x, y, z, xA, yA, zA, a,
                    static_cast<int> (surfaceZYXA.neighbors.size())};
                d->atNames.push_back(temp);

                for (auto &nb_int : surfaceZYXA.neighbors) {
                    float xNb = scaling * (d->xs * nb_int.x + cellAtoms[nb_int.type].x);
                    float yNb = scaling * (d->ys * nb_int.y + cellAtoms[nb_int.type].y);
                    float zNb = scaling * (-d->zs * nb_int.z - cellAtoms[nb_int.type].z);

                    Bond bondT = {xA, yA, zA,
                        xNb, yNb, zNb};
                    d->bonds.push_back(bondT);
                }
            }
        }
    }
}

void Render::createSurfacePoints()
{
    auto& surface = *d->surfaceXYZ;
    auto const& z_min = d->z_min;
    const size_t dX = surface[z_min][0].size() - 3;
    const size_t dY = surface[z_min].size() - 3;

    Coords3D emptyP = {-1.0, -1.0, -1.0};
    Cells points(dY, Atoms(dX, emptyP));

    for (size_t z = z_min; z < surface.size() - 2; ++z)
        for (size_t y = surface[z].size() - 2; --y >= 2;)
            for (size_t x = surface[z][y].size() - 2; --x >= 2;)
                for (const auto& atom : surface[z][y][x])
                if (!atom.deleted) {
                    auto& atom_ = points[y - 2].atoms[x - 2];
                    if (cmp_float(atom_.x, -1.0)) {
                        atom_.x = d->scaling * (x - 2) * d->xs;
                        atom_.z = -d->scaling * z * d->zs;
                        atom_.y = d->scaling * (y - 2) * d->ys;
                        }
                }      

    for (size_t i = 0; i < dY; ++i) {
        points[i].atoms[dX - 1].x = points[i].atoms[dX - 2].x +
            d->scaling * d->xs;
        points[i].atoms[dX - 1].y = points[i].atoms[dX - 2].y;
        points[i].atoms[dX - 1].z = points[i].atoms[dX - 2].z;
    }

    for (size_t i = 0; i < dX; ++i) {
        points[dY - 1].atoms[i].x = points[dY - 2].atoms[i].x;
        points[dY - 1].atoms[i].y = points[dY - 2].atoms[i].y +
            d->scaling * d->ys;
        points[dY - 1].atoms[i].z = points[dY - 2].atoms[i].z;
    }
    d->surfVertex.reserve((d->SIZE_Y - 4)*(d->SIZE_X - 4));
    d->surfNormals.reserve((d->SIZE_Y - 4)*(d->SIZE_X - 4));

    for (int i = 0; i < d->SIZE_Y - 4; i++)
        for (int j = 0; j < d->SIZE_X - 4; j++) {
            Cell& cell = points[i];
            Cell& cellNext = points[i + 1];
            Coords3D const& v1 = cell.atoms[j];
            Coords3D const v2 = {cell.atoms[j + 1].x, cell.atoms[j + 1].y, cell.atoms[j].z};
            Coords3D const v3 = {cellNext.atoms[j + 1].x, cellNext.atoms[j + 1].y, cell.atoms[j].z};
            Coords3D const& v4 = cellNext.atoms[j];

            if (!cmp_float(cell.atoms[j].z, cell.atoms[j + 1].z)) {
                Coords3D const v5 = {cell.atoms[j + 1].x, cell.atoms[j + 1].y, cell.atoms[j].z};
                Coords3D const& v6 = cell.atoms[j + 1];
                Coords3D const& v7 = cellNext.atoms[j + 1];
                Coords3D const v8 = {cellNext.atoms[j + 1].x, cellNext.atoms[j + 1].y, cell.atoms[j].z};

                d->surfNormals.insert(d->surfNormals.end(), 4, normcrossprod(v6 + -1 * v5, v7 + -1 * v5));

                d->surfVertex.push_back(v5);
                d->surfVertex.push_back(v6);
                d->surfVertex.push_back(v7);
                d->surfVertex.push_back(v8);
            }

            if (!cmp_float(cell.atoms[j].z, cellNext.atoms[j].z)) {
                Coords3D const v9 = {cellNext.atoms[j + 1].x, cellNext.atoms[j + 1].y, cell.atoms[j].z};
                Coords3D const v10 = {cellNext.atoms[j + 1].x, cellNext.atoms[j + 1].y, cellNext.atoms[j].z};
                Coords3D const& v11 = cellNext.atoms[j];
                Coords3D const v12 = {cellNext.atoms[j].x, cellNext.atoms[j].y, cell.atoms[j].z};

                d->surfNormals.insert(d->surfNormals.end(), 4, normcrossprod(v10 + -1 * v9, v11 + -1 * v9));

                d->surfVertex.push_back(v9);
                d->surfVertex.push_back(v10);
                d->surfVertex.push_back(v11);
                d->surfVertex.push_back(v12);
            }

            d->surfVertex.push_back(v1);
            d->surfVertex.push_back(v2);
            d->surfVertex.push_back(v3);
            d->surfVertex.push_back(v4);
            d->surfNormals.insert(d->surfNormals.end(), 4, normcrossprod(v2 + -1 * v1, v3 + -1 * v1));
        }

        glBindBufferARB(GL_ARRAY_BUFFER, 1);
        glBufferDataARB(GL_ARRAY_BUFFER, (d->surfVertex.size() + d->surfNormals.size()) * 3 * sizeof (float),
            0, GL_STATIC_DRAW);
        glBufferSubDataARB(GL_ARRAY_BUFFER, 0, d->surfVertex.size() * 3 * sizeof (float),
            &d->surfVertex[0].x);
        glBufferSubDataARB(GL_ARRAY_BUFFER, d->surfVertex.size() * 3 * sizeof (float),
            d->surfNormals.size() * 3 * sizeof (float), &d->surfNormals[0].x);
        glBindBufferARB(GL_ARRAY_BUFFER, 0);
}

void
Render::drawAxis() {
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
    QString xText = "x (" + QString::number(d->Vx->x)
        + "," + QString::number(d->Vx->y) + "," + QString::number(d->Vx->z) + ")";
    QString yText = "y (" + QString::number(d->Vy->x)
        + "," + QString::number(d->Vy->y) + "," + QString::number(d->Vy->z) + ")";
    QString zText = "z (" + QString::number(d->Vz->x)
        + "," + QString::number(d->Vz->y) + "," + QString::number(d->Vz->z) + ")";
    this->renderText(4.2, -0.1, 0.0, xText);
    glColor3f(1.0, 0.0, 1.0);
    this->renderText(0.0, -0.1, -4.3, zText);
    glColor3f(1.0, 1.0, 1.0);
    this->renderText(0.0, 4.3, 0.0, yText);
}

void
Render::initializeGL() {
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    d->theSphere = glGenLists(1);
    glNewList(d->theSphere, GL_COMPILE);
    sphereTemplate(d->sR * 10);
    glEndList();
}

void
Render::resizeGL(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40.0,
            static_cast<GLfloat> (width) / static_cast<GLfloat> (height),
            1.0, 80.0);
    glMatrixMode(GL_MODELVIEW);
}

void
Render::initMatrix() {
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
        1.0};
    d->matrix = matrix1;
}

void
Render::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    draw();
}

void
Render::draw() {
    if (d->dataInitialized) {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        if (d->visualizationType == GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK
            || d->visualizationType == GRES::VizType::ATOMS_SURFACE_AND_BULK
            || d->visualizationType == GRES::VizType::ATOMS_AND_BONDS_SURFACE
            || d->visualizationType == GRES::VizType::ATOMS_SURFACE) {/*
			GLfloat ambientLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };*/
            GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
            GLfloat diffuseLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
            GLfloat specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
            GLfloat light_position[] = {10.0, 10.0, 20.0, 0.0f};
            //GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
            //GLfloat mat_shininess[] = { 200.0 };
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

            //GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
            //GLfloat mat_shininess[] = { 100.0 };
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
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

            setGeometry(-20.0 - d->zs * d->scaling * (d->z_center - d->z_min));

            glTranslatef(-6.0, -6.0, -20.0);
            drawAxis();
            glTranslatef(6.0, 6.0, 20.0);

            glTranslatef(-5.0 - 2.0 * d->xs * d->scaling,
                -5.0 - 2.0 * d->ys * d->scaling, -20.0);
            glColor3f(0.98f, 0.625f, 0.12f);
            //sphere( 5,10,10 );

            if (!(d->visualizationType == GRES::VizType::ATOMS_SURFACE_AND_BULK
                || d->visualizationType == GRES::VizType::ATOMS_SURFACE)
                || d->rotating || d->scribble || d->moving) {
                glColor3f(0, 0, 1.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, &d->bonds[0].x1);
                glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(2 * d->bonds.size()));
                glDisableClientState(GL_VERTEX_ARRAY);
            }

            if (!d->rotating && !d->scribble && !d->moving) {
                int name = 0;
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                //t.start();
                for (auto const& atomInfo : d->atNames) {
                    glColor3f(105.0f / 255, 68.0f / 255, 9.0f / 255);
                    if (atomInfo.fNbCount == 3)
                        glColor3f(0.08f, 0.925f, 0.02f);
                    if (atomInfo.fNbCount == 2)
                        glColor3f(1.0, 1.0, 0.0);
                    if (atomInfo.fNbCount == 1)
                        glColor3f(1.0, 0.0, 0.0);
                    if (cmp_float(atomInfo.x, d->selAtomType.x)
                        && cmp_float(atomInfo.y, d->selAtomType.y)
                        && cmp_float(atomInfo.z, d->selAtomType.z))
                        glColor3f(1.0, 0.0, 1.0);
                    glPushMatrix();
                    glTranslatef(atomInfo.x, atomInfo.y, atomInfo.z);
                    ++name;
                    glLoadName(name);

                    glBindBufferARB(GL_ARRAY_BUFFER, 1);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *d->vSize1));

                    glDrawArrays(GL_TRIANGLE_FAN, 0, d->vSize1);

                    glBindBufferARB(GL_ARRAY_BUFFER, 2);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *d->vSize2));

                    glDrawArrays(GL_TRIANGLE_FAN, 0, d->vSize2);

                    glBindBufferARB(GL_ARRAY_BUFFER, 3);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *d->vSize3));

                    glDrawArrays(GL_QUAD_STRIP, 0, d->vSize3);

                    glPopMatrix();
                }
                //	qDebug("Time elapsed: %d ms", t.elapsed());
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_NORMAL_ARRAY);
                glBindBufferARB(GL_ARRAY_BUFFER, 0);
            }
            glFlush();
        }
        else if (d->visualizationType == GRES::VizType::CELLS_SURFACE) {

            GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
            GLfloat diffuseLight[] = {0.7f, 0.7f, 0.7f, 1.0f};
            GLfloat specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
            GLfloat light_position[] = {5.0f, 5.0f, 5.0f, 0.0f};
            //GLfloat light_position[] = { 10.0, 10.0, 20.0, 0.0f };
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

            glTranslatef(-5.0, -5.0, -20.0 + d->zs * d->z_min * d->scaling);

            glPointSize(10.0f);
            glDisable(GL_CULL_FACE);

            glBindBufferARB(GL_ARRAY_BUFFER, 1);
            glVertexPointer(3, GL_FLOAT, 0, 0);
            glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *d->surfVertex.size()));

            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glDrawArrays(GL_QUADS, 0, static_cast<GLsizei>(d->surfVertex.size()));
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

void
Render::mousePressEvent(QMouseEvent *event) {
    d->lastPos = event->pos();
    if (event->buttons() & Qt::LeftButton) {
        if (event->modifiers() == Qt::ControlModifier)
            processSelection(event->x(), event->y());
        else
            d->rotating = true;
    } else if (event->buttons() & Qt::MidButton) {
        d->moving = true;
    }
}

void
Render::mouseReleaseEvent(QMouseEvent *event) {
    if ((!d->scribble) && (event->button() == Qt::RightButton)) {
        processSelectionMenu();
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

void
Render::mouseMoveEvent(QMouseEvent *event) {
    GLfloat dx = static_cast<GLfloat>(event->x() - d->lastPos.x()) / width();
    GLfloat dy = static_cast<GLfloat>(event->y() - d->lastPos.y()) / height();
    if (event->buttons() & Qt::LeftButton) {
        d->drotationX = 180 * dy;
        d->drotationY = 180 * dx;
        d->rotating = true;
        updateGL();
    } else if (event->buttons() & Qt::RightButton) {
        d->scribble = true;
        d->scale += dy * d->scale;
        updateGL();
    } else if (event->buttons() & Qt::MidButton) {

        d->movX = 18 * dx;
        d->movY = 18 * dy;
        d->sumMovX += d->movX;
        d->sumMovY += d->movY;
        updateGL();
    }

    d->lastPos = event->pos();
}

void
Render::keyPressEvent(QKeyEvent * event) {
    if (event->text() == "p")
        emit(etching());
}

void
Render::processSelectionMenu() {
    d->selAtomMenu->setWindowFlags(Qt::Tool);
    d->selAtomMenu->setWindowTitle("Atom`s information");
    d->selAtomMenu->move(QCursor::pos());
    d->selAtomMenu->setInfo(d->selAtomType.xC, d->selAtomType.yC,
        d->selAtomType.zC, d->selAtomType.type, d->selAtomType.fNbCount);
    d->selAtomMenu->show();
}

void
Render::saveResult() {
    QPixmap outPixmap;
    this->raise();
    outPixmap = QPixmap::grabWindow(this->winId());

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
            "./untitled.png",
            tr("Images (*.png *.jpg *.tiff *.bmp *.xpm)"));
    if (!fileName.isNull())//проверяем не отменил ли юзер выбор файла
    {
        if (!outPixmap.save(fileName))
            QMessageBox::information(this, tr("GRES"), tr("Cannot save %1.").arg(fileName));
    }
}

void
Render::setGeometry(GLfloat zCenter) {
    glTranslatef(0, 0, zCenter);
    if (!cmp_float(d->drotationX, 0))
        glRotatef(d->drotationX, 1.0, 0.0, 0.0);
    if (!cmp_float(d->drotationY, 0))
        glRotatef(d->drotationY, 0.0, 1.0, 0.0);
    glTranslatef(0, 0, -zCenter);

    d->drotationX = d->drotationY = 0.0;
    glMultMatrixf(&d->matrix[0]);
    glGetFloatv(GL_MODELVIEW_MATRIX, &d->matrix[0]);

    glLoadIdentity();

    glTranslatef(0, 0, zCenter);

    if (!cmp_float(d->scale, 1.0))
        glScalef(d->scale, d->scale, d->scale);
    glTranslatef(0, 0, -zCenter);
    glTranslatef(d->sumMovX, -d->sumMovY, 0.0);
    glMultMatrixf(&d->matrix[0]);
}

void
Render::sphereTemplate(float R) {

    GLUquadricObj *quadobj;
    quadobj = gluNewQuadric();
    gluSphere(quadobj, R, 10, 10);
    gluDeleteQuadric(quadobj);
    //  	glutSolidSphere(sR-0.06, 10, 10);
    // glutSolidCube(2*sR-0.04);
    // glutSolidCube(sR-0.04);
}

QSize
Render::minimumSizeHint() const {

    return QSize(400, 400);
}

QSize
Render::sizeHint() const {

    return QSize(600, 600);
}

void
Render::processSelection(int x, int y) {
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
    fAspect = static_cast<float> (viewport[2]) / static_cast<float> (viewport[3]);
    gluPerspective(40.0f, fAspect, 1.0, 80.0);
    paintGL();
    hits = glRenderMode(GL_RENDER);

    if (hits == 1)
        processAtom(selectBuff);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    updateGL();
}

void Render::processAtom(const GLuint *pSelectBuff) {
    int id;
    id = pSelectBuff[3];
    for (auto const& atName : d->atNames) {
        if (atName.name == id) {
            d->selAtomType = atName;
            break;
        }
    }
}

void Render::view(Surface3DPtr surface, Cell &atTypes, float Xsize, float Ysize,
    float Zsize, int center, int min, int width, int height, Coords3D &vX,
    Coords3D &vY, Coords3D &vZ, GRES::VizType vT)
{
    d->visualizationType = vT;

    d->xs = Xsize;
    d->ys = Ysize;
    d->zs = Zsize;
    d->surfaceXYZ = surface;
    d->z_min = min;
    d->z_center = center;
    d->cellAtoms = atTypes;
    d->SIZE_X = width;
    d->SIZE_Y = height;
    d->Vx = &vX;
    d->Vy = &vY;
    d->Vz = &vZ;
    d->dataInitialized = true;
    d->dataChanged = true;

    const float W = 10;
    const float H = 10;

    float h = (d->SIZE_Y - 4) * d->ys;
    float w = (d->SIZE_X - 4) * d->xs;

    d->scaling = w > h ? W / w : H / h;

    d->surfPoints.clear();
    d->coordsOfAtoms.clear();
    d->bonds.clear();
    d->atNames.clear();
    d->surfVertex.clear();
    d->surfNormals.clear();
    initMatrix();
    changeVizType(vT);
}

Render::~Render() {
}
