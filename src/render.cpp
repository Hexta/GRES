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
#include "render.h"
#include "selectAtomMenu.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#include <QTime>
#include <QMenu>
#include <QFileDialog>
#include <QMessageBox>
#include "geometry.h"
#include "QDebug"

#define BUFFER_OFFSET(i) ((char*)NULL + (i))
#ifdef _WIN32
#define glBindBufferARB           pglBindBufferARB
#define glDeleteBuffersARB        pglDeleteBuffersARB
#define glBufferDataARB           pglBufferDataARB
#define glBufferSubDataARB        pglBufferSubDataARB
#endif

QTime t;

Render::Render(QWidget *parent) : QGLWidget(parent) {
    scribble =
            rotating =
            moving = false;

    drotationX =
            drotationY =
            drotationZ =
            movX =
            movY =
            sumMovX =
            sumMovY = 0.0;

    sphereQual = 3;
    scale = 1.0;
    sR = 0.05F;

    selAtomType.x =
            selAtomType.y =
            selAtomType.z =
            selAtomType.type =
            selAtomType.fNbCount = -1;
    createActions();

    SIZE_Y =
            SIZE_X = 14;

    xs =
            ys =
            zs = 1;

    z_min = 0;
    z_center = 2;
    dataInitialized = false;
    dataChanged = false;
    buffers.clear();
    visualType = GRES::VizType::CELLS_SURFACE;
}

void
Render::changeVizType(GRES::VizType type) {
    visualType = type;

    if (!dataInitialized
        || !(dataChanged || coordsOfAtoms.empty() || surfPoints.empty())) {
        updateGL();
        dataChanged = false;
        return;
    }

    atNames.clear();
    bonds.clear();
    surfVertex.clear();
    surfNormals.clear();

    if (!buffers.empty())
        glDeleteBuffersARB(1, &buffers[0]);

    switch (visualType) {
        case GRES::VizType::CELLS_SURFACE:
            if (!surfPoints.empty()) {
                break;
            }

            createSurfacePoints(*surfaceXYZ, xs, ys, zs, z_min);
            buffers.clear();
            buffers.push_back(1);
            break;
        case GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK:
        {
            if (!coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBonds(*surfaceXYZ, *cellAtoms, xs, ys, zs, z_min,
                    atNames, bonds);
            createSphere(0.09 * scaling, 10, 10, vSize1, vSize2, vSize3);
            buffers.clear();
            buffers.push_back(1);
            buffers.push_back(2);
            buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_AND_BONDS_SURFACE:
        {
            if (!coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBondes(*surfaceXYZ, *surfAtoms, *cellAtoms,
                    xs, ys, zs, z_min, scaling, atNames,
                    bonds);
            createSphere(0.09 * scaling, 10, 10, vSize1, vSize2, vSize3);
            buffers.clear();
            buffers.push_back(1);
            buffers.push_back(2);
            buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_SURFACE_AND_BULK:
        {
            if (!coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBonds(*surfaceXYZ, *cellAtoms, xs, ys, zs,
                    z_min, atNames, bonds);
            createSphere(0.2 * scaling, 10, 10, vSize1, vSize2, vSize3);
            buffers.clear();
            buffers.push_back(1);
            buffers.push_back(2);
            buffers.push_back(3);
        }
            break;
        case GRES::VizType::ATOMS_SURFACE:
        {
            if (!coordsOfAtoms.empty()) {
                break;
            }

            createAtomsAndBondes(*surfaceXYZ, *surfAtoms, *cellAtoms,
                    xs, ys, zs, z_min, scaling, atNames,
                    bonds);
            createSphere(0.2 * scaling, 10, 10, vSize1, vSize2, vSize3);
            buffers.clear();
            buffers.push_back(1);
            buffers.push_back(2);
            buffers.push_back(3);
        }
            break;
        default: break;
    }

    updateGL();
    dataChanged = false;
}

void
Render::createActions() {
    exitAction = new QAction(QIcon(":/images/exit.png"), tr("E&xit"), this);
    // 	connect (exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));
}

void
Render::createAtomsAndBonds(surface3D &surface, atomsCoords &cellAts, float xs_,
                            float ys_, float zs_, int z_min, AtomsNames &atN,
                            Bonds &outBonds) {

    int name = 0;
    for (size_t z = z_min; z < surface.size() - 2; ++z)
        for (size_t y = surface[z].size() - 2; --y >= 2;)
            for (size_t x = surface[z][y].size() - 2; --x >= 2;) {
                float x0 = scaling * x*xs_;
                float y0 = scaling * y*ys_;
                float z0 = -scaling * (z - z_min) * zs_;
                const char atomsCount = static_cast<char>(surface[z][y][x].size());
                for (unsigned char a = atomsCount; --a > 0;) {

                    if (surface[z][y][x][a].fNbCount) {
                        float xA = x0 + scaling * cellAts[a].x;
                        float yA = y0 + scaling * cellAts[a].y;
                        float zA = z0 - scaling * cellAts[a].z;

                        ++name;
                        atomName temp = {name, static_cast<int> (x),
                                         static_cast<int> (y),
                                         static_cast<int> (z), xA, yA, zA, a,
                                         static_cast<int> (atomsCount)};
                        atN.push_back(temp);

                        for (auto &nb : surface[z][y][x][a].neighbours) {
                            float xNb = scaling * (xs * nb.x + cellAts[nb.type].x);
                            float yNb = scaling * (ys * nb.y + cellAts[nb.type].y);
                            float zNb = scaling * (-zs * nb.z - cellAts[nb.type].z);

                            Bond bondT = {xA, yA, zA, xNb, yNb, zNb};
                            outBonds . push_back(bondT);
                        }
                    }
                }
            }
}

void Render::createSurfacePoints(const surface3D &surface, float Xsize, float Ysize, float Zsize,
                            int z_min) {
    const size_t dX = surface[z_min][0].size() - 3;
    const size_t dY = surface[z_min].size() - 3;

    coords3D emptyP = {-1.0, -1.0, -1.0};
    cells points(dY, atomsCoords(dX, emptyP));

    for (size_t z = z_min; z < surface.size() - 2; ++z)
        for (size_t y = surface[z].size() - 2; --y >= 2;)
            for (size_t x = surface[z][y].size() - 2; --x >= 2;)
                for (const auto& atom : surface[z][y][x])
                    if (!atom.deleted)
                        if (cmp_float(points[y - 2][x - 2].x, -1.0)) {
                            points[y - 2][x - 2].x = scaling * (x - 2) * Xsize;
                            points[y - 2][x - 2].z = -scaling * z*Zsize;
                            points[y - 2][x - 2].y = scaling * (y - 2) * Ysize;
                        }

    for (int i = 0; i < dY; ++i) {
        points[i][dX - 1].x = points[i][dX - 2].x + scaling*Xsize;
        points[i][dX - 1].y = points[i][dX - 2].y;
        points[i][dX - 1].z = points[i][dX - 2].z;
    }

    for (int i = 0; i < dX; ++i) {
        points[dY - 1][i].x = points[dY - 2][i].x;
        points[dY - 1][i].y = points[dY - 2][i].y + scaling*Ysize;
        points[dY - 1][i].z = points[dY - 2][i].z;
    }
    surfVertex.reserve((SIZE_Y - 4)*(SIZE_X - 4));
    surfNormals.reserve((SIZE_Y - 4)*(SIZE_X - 4));

    for (int i = 0; i < SIZE_Y - 4; i++)
        for (int j = 0; j < SIZE_X - 4; j++) {
            coords3D v1 = {points[i][j].x, points[i][j].y, points[i][j].z};
            coords3D v2 = {points[i][j + 1].x, points[i][j + 1].y, points[i][j].z};
            coords3D v3 = {points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i][j].z};
            coords3D v4 = {points[i + 1][j].x, points[i + 1][j].y, points[i][j].z};

            if (!cmp_float(points[i][j].z, points[i][j + 1].z)) {
                coords3D v5 = {points[i][j + 1].x, points[i][j + 1].y, points[i][j].z};
                coords3D v6 = {points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z};
                coords3D v7 = {points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i][j + 1].z};
                coords3D v8 = {points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i][j].z};

                surfNormals.insert(surfNormals.end(), 4, normcrossprod(v6 + -1 * v5, v7 + -1 * v5));

                surfVertex.push_back(v5);
                surfVertex.push_back(v6);
                surfVertex.push_back(v7);
                surfVertex.push_back(v8);
            }

            if (!cmp_float(points[i][j].z, points[i + 1][j].z)) {
                coords3D v9 = {points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i][j].z};
                coords3D v10 = {points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j].z};
                coords3D v11 = {points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z};
                coords3D v12 = {points[i + 1][j].x, points[i + 1][j].y, points[i][j].z};

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
    glBufferDataARB(GL_ARRAY_BUFFER, (surfVertex.size() + surfNormals.size())*3 * sizeof (float),
            0, GL_STATIC_DRAW);
    glBufferSubDataARB(GL_ARRAY_BUFFER, 0, surfVertex.size()*3 * sizeof (float),
            &surfVertex[0].x);
    glBufferSubDataARB(GL_ARRAY_BUFFER, surfVertex.size()*3 * sizeof (float),
            surfNormals.size()*3 * sizeof (float), &surfNormals[0].x);
    glBindBufferARB(GL_ARRAY_BUFFER, 0);
}

void
Render::drawAxis() {
    /*
    Рисует оси координат
     */
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
    QString xText = "x (" + QString::number(Vx->x)
            + "," + QString::number(Vx->y) + "," + QString::number(Vx->z) + ")";
    QString yText = "y (" + QString::number(Vy->x)
            + "," + QString::number(Vy->y) + "," + QString::number(Vy->z) + ")";
    QString zText = "z (" + QString::number(Vz->x)
            + "," + QString::number(Vz->y) + "," + QString::number(Vz->z) + ")";
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

    theSphere = glGenLists(1);
    glNewList(theSphere, GL_COMPILE);
    sphereTemplate(sR * 10);
    glEndList();
#ifdef _WIN32
    pglBindBufferARB = (PFNGLBINDBUFFERARBPROC) wglGetProcAddress("glBindBufferARB");
    pglDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) wglGetProcAddress("glDeleteBuffersARB");
    pglBufferDataARB = (PFNGLBUFFERDATAARBPROC) wglGetProcAddress("glBufferDataARB");
    pglBufferSubDataARB = (PFNGLBUFFERSUBDATAARBPROC) wglGetProcAddress("glBufferSubDataARB");
#endif
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
Render::initMatrix(vector<GLfloat>* m) {
    vector<GLfloat> matrix1 = {1.0,
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
    *m = matrix1;
}

void
Render::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    draw();
}

void
Render::draw() {
    if (dataInitialized) {
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        if (visualType == GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK
            || visualType == GRES::VizType::ATOMS_SURFACE_AND_BULK
            || visualType == GRES::VizType::ATOMS_AND_BONDS_SURFACE
            || visualType == GRES::VizType::ATOMS_SURFACE) {/*
			GLfloat ambientLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
			GLfloat diffuseLight[] = { 0.8, 0.8, 0.8, 1.0 };*/
            GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
            GLfloat diffuseLight[] = {0.7, 0.7, 0.7, 1.0};
            GLfloat specularLight[] = {1.0, 1.0, 1.0, 1.0};
            GLfloat light_position[] = {10.0, 10.0, 20.0, 0.0f};
            //GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
            //GLfloat mat_shininess[] = { 200.0 };
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

            //GLfloat mat_specular[] = { 0.8, 0.8, 0.8, 1.0 };
            //GLfloat mat_shininess[] = { 100.0 };
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
            //glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
            GLfloat mat_specular[] = {0.3, 0.3, 0.3, 1.0};
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

            glTranslatef(-5.0 - 2.0 * xs*scaling, -5.0 - 2.0 * ys*scaling, -20.0);
            glColor3f(0.98, 0.625, 0.12);
            //sphere( 5,10,10 );

            if (!(visualType == GRES::VizType::ATOMS_SURFACE_AND_BULK || visualType == GRES::VizType::ATOMS_SURFACE)
                || rotating || scribble || moving) {
                glColor3f(0, 0, 1.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, &bonds[0].x1);
                glDrawArrays(GL_LINES, 0, static_cast<GLsizei> (2 * bonds.size()));
                glDisableClientState(GL_VERTEX_ARRAY);
            }

            if (!rotating && !scribble && !moving) {
                int name = 0;
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                //t.start();
                for (unsigned int i = 0; i < atNames.size(); ++i) {
                    glColor3f(105.0 / 255, 68.0 / 255, 9.0 / 255);
                    if (atNames[i].fNbCount == 3)
                        glColor3f(0.08, 0.925, 0.02);
                    if (atNames[i].fNbCount == 2)
                        glColor3f(1.0, 1.0, 0.0);
                    if (atNames[i].fNbCount == 1)
                        glColor3f(1.0, 0.0, 0.0);
                    if (cmp_float(atNames[i].x, selAtomType.x)
                        && cmp_float(atNames[i].y, selAtomType.y)
                        && cmp_float(atNames[i].z, selAtomType.z))
                        glColor3f(1.0, 0.0, 1.0);
                    glPushMatrix();
                    glTranslatef(atNames[i].x, atNames[i].y, atNames[i].z);
                    ++name;
                    glLoadName(name);

                    glBindBufferARB(GL_ARRAY_BUFFER, 1);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *vSize1));

                    glDrawArrays(GL_TRIANGLE_FAN, 0, vSize1);

                    glBindBufferARB(GL_ARRAY_BUFFER, 2);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *vSize2));

                    glDrawArrays(GL_TRIANGLE_FAN, 0, vSize2);

                    glBindBufferARB(GL_ARRAY_BUFFER, 3);
                    glVertexPointer(3, GL_FLOAT, 0, 0);
                    glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *vSize3));

                    glDrawArrays(GL_QUAD_STRIP, 0, vSize3);

                    glPopMatrix();
                }
                //	qDebug("Time elapsed: %d ms", t.elapsed());
                glDisableClientState(GL_VERTEX_ARRAY);
                glDisableClientState(GL_NORMAL_ARRAY);
                glBindBufferARB(GL_ARRAY_BUFFER, 0);
            }
            glFlush();
        }/*-------------------------------*/
        else if (visualType == GRES::VizType::CELLS_SURFACE) {

            GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
            GLfloat diffuseLight[] = {0.7, 0.7, 0.7, 1.0};
            GLfloat specularLight[] = {1.0, 1.0, 1.0, 1.0};
            GLfloat light_position[] = {5.0, 5.0, 5.0, 0.0f};
            //GLfloat light_position[] = { 10.0, 10.0, 20.0, 0.0f };
            GLfloat mat_specular[] = {0.3, 0.3, 0.3, 1.0};
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

            glColor3f(0.18, 0.325, 0.82);

            glTranslatef(-5.0, -5.0, -20.0 + zs * z_min * scaling);

            glPointSize(10.0f);
            glDisable(GL_CULL_FACE);

            glBindBufferARB(GL_ARRAY_BUFFER, 1);
            glVertexPointer(3, GL_FLOAT, 0, 0);
            glNormalPointer(GL_FLOAT, 0, BUFFER_OFFSET(3 * sizeof (float) *surfVertex.size()));

            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glDrawArrays(GL_QUADS, 0, static_cast<GLsizei> (surfVertex.size()));
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
    lastPos = event->pos();
    if (event->buttons() & Qt::LeftButton) {
        if (event->modifiers() == Qt::ControlModifier)
            processSelection(event->x(), event->y());
        else
            rotating = true;
    } else if (event->buttons() & Qt::MidButton) {

        moving = true;
    }
}

void
Render::mouseReleaseEvent(QMouseEvent *event) {
    if ((!scribble) && (event->button() == Qt::RightButton)) {
        processSelectionMenu();
    }
    if ((scribble) && (event->button() == Qt::RightButton)) {
        scribble = false;
        updateGL();
    }
    if (rotating && event->button() == Qt::LeftButton) {
        rotating = false;
        updateGL();
    }
    if (moving && event->button() == Qt::MidButton) {

        moving = false;
        updateGL();
    }
}

void
Render::mouseMoveEvent(QMouseEvent *event) {
    GLfloat dx = static_cast<GLfloat>(event->x() - lastPos.x()) / width();
    GLfloat dy = static_cast<GLfloat>(event->y() - lastPos.y()) / height();
    if (event->buttons() & Qt::LeftButton) {
        drotationX = 180 * dy;
        drotationY = 180 * dx;
        rotating = true;
        updateGL();
    } else if (event->buttons() & Qt::RightButton) {
        scribble = true;
        scale += dy*scale;
        updateGL();
    } else if (event->buttons() & Qt::MidButton) {

        movX = 18 * dx;
        movY = 18 * dy;
        sumMovX += movX;
        sumMovY += movY;
        updateGL();
    }

    lastPos = event->pos();
}

void
Render::keyPressEvent(QKeyEvent * event) {
    if (event->text() == "p")
        emit(etching());
}

void
Render::processSelectionMenu() {
    selAtomMenu = new SelectAtomMenu;
    selAtomMenu->setWindowFlags(Qt::Tool);
    selAtomMenu->setWindowTitle("Atom`s information");
    selAtomMenu->move(QCursor::pos());
    selAtomMenu->setInfo(selAtomType.xC, selAtomType.yC, selAtomType.zC,
            selAtomType.type, selAtomType.fNbCount);
    selAtomMenu->show();
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
    for (auto &atName : atNames) {
        if (atName.name == id) {
            selAtomType.xC = atName.xC;
            selAtomType.yC = atName.yC;
            selAtomType.zC = atName.zC;
            selAtomType.x = atName.x;
            selAtomType.y = atName.y;
            selAtomType.z = atName.z;
            selAtomType.type = atName.type;
            selAtomType.fNbCount = static_cast<char> (atName.fNbCount);
            break;
        }
    }
}

void
Render::view(surface3D &surface, vector<atomType> &surfAt, atomsCoords &atTypes,
             float Xsize, float Ysize, float Zsize, int center, int min,
             int width, int height, coords3D &vX, coords3D &vY, coords3D &vZ,
             GRES::VizType vT) {
    visualType = vT;
    surfAtoms = &surfAt;

    xs = Xsize;
    ys = Ysize;
    zs = Zsize;
    surfaceXYZ = &surface;
    z_min = min;
    z_center = center;
    cellAtoms = &atTypes;
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
    initMatrix(&matrix);
    changeVizType(vT);
}
