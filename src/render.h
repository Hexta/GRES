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

#include <QWidget>
#include <QGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QColorDialog>
#include <QVector>
#include <vector>
#include <list>
#include <cmath>
#include "functions.h"
#include "consts.h"
#include "dataTypes.h"
#ifdef _WIN32
#include <windows.h>
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif
using namespace std;

class SelectAtomMenu;

class Render : public QGLWidget {
    Q_OBJECT
public:
    Render(QWidget *parent = 0);
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
public slots:
    void view(surface3D &surface, vector<atomType>&, atomsCoords&, float Xsize,
            float Ysize, float Zsize, int z_min, int z_center, int width,
            int height, coords3D &Vx, coords3D &Vy, coords3D &Vz,
            GRES::VizType vizualType);
    void changeVizType(GRES::VizType type);
    void saveResult();
signals:
    void etching(void);
protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent * event);
private:
    void draw();
    GLfloat drotationX;
    GLfloat drotationY;
    GLfloat drotationZ;
    GLfloat movX;
    GLfloat sumMovX;
    GLfloat sumMovY;
    GLfloat movY;
    GLfloat scale;
    QPoint lastPos;
    SelectAtomMenu* selAtomMenu;

    struct SelAtomType {
        short xC;
        short yC;
        short zC;
        float x;
        float y;
        float z;
        char type;
        char fNbCount;
    } selAtomType;

    struct TypeAndCoordsOfAtom {
        float x;
        float y;
        float z;
        char fNbCount;
    };
    vector<TypeAndCoordsOfAtom> typeAndCoordsOfAtoms;
    AtomsNames atNames;
    float sR;
    void sphereTemplate(float);
    void cylinderTemplate(float);
    coords3D *Vx, *Vy, *Vz;
    GLuint theSphere;
    QColor clearColor;
    atomsCoords *cellAtoms;
    surface3D *surfaceXYZ;
    float xs, ys, zs;
    int z_center, z_min;
    int SIZE_X, SIZE_Y;
    bool scribble;
    bool moving;
    bool rotating;
    bool dataChanged;
    float scaling;
    cellAtom coordsOfAtoms; //список координат атомов
    Bonds bonds;
    QAction *exitAction;
    void createActions();
    void processSelection(int x, int y);
    void processSelectionMenu();
    void processAtom(GLuint *pSelectBuff);
    bool dataInitialized;
    GRES::VizType visualType; //тип визуализации
    cells surfPoints;
    atomsCoords surfVertex;
    atomsCoords surfNormals;
    vector <GLuint> buffers;
    int sphereQual;
    int vSize1, vSize2, vSize3;
    vector<GLfloat> matrix;
    void createAtomsAndBonds(surface3D &surface, atomsCoords&, float Xsize,
            float Ysize, float Zsize, int z_min, AtomsNames&, Bonds&);
    void createSurfacePoints(surface3D &surface, float Xsize, float Ysize,
            float Zsize, int z_min, atomsCoords &surfV, atomsCoords &surfN);
    void initMatrix(vector<GLfloat>*);
    void drawAxis();
    void setGeometry(GLfloat zCenter = 0);
    vector<atomType> *surfAtoms;

#ifdef _WIN32
    PFNGLBINDBUFFERARBPROC pglBindBufferARB;
    PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB;
    PFNGLBUFFERDATAARBPROC pglBufferDataARB;
    PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB;
#endif
};