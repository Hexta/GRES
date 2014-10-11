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

#define _STDCALL_SUPPORTED
#define GLUT_NO_LIB_PRAGMA

#define NOMINMAX
#include "MainWindow.h"

#include <QApplication>

#include <QGLFormat>

int main(int argc, char* argv[]) {
    Q_INIT_RESOURCE(resources);
    QApplication app(argc, argv);
#ifdef _WIN32
    QGLFormat fmt;
    fmt.setAlpha(true);
    fmt.setDoubleBuffer(false);
    QGLFormat::setDefaultFormat(fmt);
#endif
    MainWindow* content;
    content = new MainWindow(0, argc, argv);
    content->setWindowTitle(QObject::tr("GRES"));
    content->show();
    return app.exec();
}
