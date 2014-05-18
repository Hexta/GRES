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

#include "functions.h"
#include "consts.h"

#include <QtGui>
#include <QMainWindow>
#include <memory>

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0, int argc = 0, char *const argv[] = NULL);
    virtual ~MainWindow();

public slots:
    void newDocument();
    void getSettings(int, int, int, int, int, int);
    void changeVizType(QAction* type);
    void etch(int, int, float *rates);

private:
    void createActions();
    void createMenus();
    void createToolBars();

private slots:
    void showEtchMenu();
    void showMenuMask();
    void setMask(std::vector<bool> inMask);

private:
    struct Private;
    std::unique_ptr<Private> d;
};
