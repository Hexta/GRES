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

#pragma once

#include "functions.h"
#include "consts.h"
#include "geometry.h"
#include "Surface3D.h"

#include <QtOpenGL/QGLWidget>

#ifdef _WIN32
#include <windows.h>
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <memory>

class QWidget;
class QMouseEvent;
class QKeyEvent;
class Cell;

class Render : public QGLWidget {
    Q_OBJECT

public:
    Render(QWidget* parent = 0);
    ~Render();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

public slots:
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
        GRES::VizType vT);
    void changeVizType(GRES::VizType type);
    void saveResult();

signals:
    void etching(void);

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);

private:
    struct Private;
    std::unique_ptr<Private> d;
};
