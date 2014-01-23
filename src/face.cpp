/******************************************************************************
 * Copyright (c) 2009-2013 Artur Molchanov <artur.molchanov@gmail.com>     *
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 ******************************************************************************/

#include "face.h"
#include "Render.h"
#include "settings.h"
#include "consts.h"
#include "etchingMenu.h"
#include "maskMenu.h"

#include <QDockWidget>
#include <QActionGroup>
#include <QMenuBar>
#include <QToolBar>

MainW::MainW(QWidget *parent, int, char * const *) : QMainWindow(parent) {
    surfaceXYZ.reserve(5000);
    h = 1;
    k = 0;
    l = 0;
    SIZE_X = 54;
    SIZE_Y = 54;
    SIZE_Z = 4;
    z_min = 0;

    createActions();
    createMenus();
    createToolBars();
    result = new Render;
    connect(saveAct, SIGNAL(triggered()), result, SLOT(saveResult()));
    setCentralWidget(result);
    result->setFocusPolicy(Qt::ClickFocus);

    settings = new Settings;
    QDockWidget *settingsDock = new QDockWidget(tr("Settings"), this);
    settingsDock->setAllowedAreas(Qt::LeftDockWidgetArea);
    settingsDock->setWidget(settings);
    addDockWidget(Qt::LeftDockWidgetArea, settingsDock);
    result->setFocus();
    setSettings();
    vizualType = GRES::VizType::CELLS_SURFACE;

    etchMenu = new EtchingMenu;

    connect(settings,
            SIGNAL(settingsChanged(int, int, int, int, int, int)),
            this,
            SLOT(getSettings(int, int, int, int, int, int)));

    //     connect (result, SIGNAL (etching()), this, SLOT (etch()));
    connect(etchMenu, SIGNAL(startEtching(int, int, float*)), this,
            SLOT(etch(int, int, float*)));
}

void
MainW::changeVizType(QAction* type) {
    if (type == viewAsAts_SurfaceAndBulkAct)
        vizualType = GRES::VizType::ATOMS_SURFACE_AND_BULK;

    else if (type == viewAsAts_SurfaceAct)
        vizualType = GRES::VizType::ATOMS_SURFACE;

    else if (type == viewAsAtsAndBonds_SurfaceAndBulkAct)
        vizualType = GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK;

    else if (type == viewAsAtsAndBonds_SurfaceAct)
        vizualType = GRES::VizType::ATOMS_AND_BONDS_SURFACE;

    else if (type == viewAsCellsSurface)
        vizualType = GRES::VizType::CELLS_SURFACE;

    result->changeVizType(vizualType);
}

void
MainW::createActions() {
    newAct = new QAction(QIcon(":/images/new.png"), tr("N&ew"), this);
    connect(newAct, SIGNAL(triggered()), this, SLOT(newDocument()));

    saveAct = new QAction(QIcon(":/images/filesave.png"), tr("S&ave"), this);
    saveAct->setEnabled(false);

    viewAsAtsAndBonds_SurfaceAndBulkAct = new QAction(tr("Atoms and bonds"), this);
    viewAsAtsAndBonds_SurfaceAndBulkAct ->setCheckable(true);

    viewAsAtsAndBonds_SurfaceAct = new QAction(tr("Atoms and bonds"), this);
    viewAsAtsAndBonds_SurfaceAct ->setCheckable(true);

    viewAsAts_SurfaceAndBulkAct = new QAction(tr("Atoms"), this);
    viewAsAts_SurfaceAndBulkAct ->setCheckable(true);

    viewAsAts_SurfaceAct = new QAction(tr("Atoms"), this);
    viewAsAts_SurfaceAct ->setCheckable(true);

    viewAsCellsSurface = new QAction(tr("Surface"), this);
    viewAsCellsSurface ->setCheckable(true);

    viewGroup = new QActionGroup(this);
    viewGroup -> addAction(viewAsAts_SurfaceAndBulkAct);
    viewGroup -> addAction(viewAsAts_SurfaceAct);
    viewGroup -> addAction(viewAsAtsAndBonds_SurfaceAct);
    viewGroup -> addAction(viewAsAtsAndBonds_SurfaceAndBulkAct);
    viewGroup -> addAction(viewAsCellsSurface);
    viewAsCellsSurface->setChecked(true);


    etchingAction = new QAction(tr("Etch"), this);
    etchingAction->setEnabled(false);
    connect(etchingAction, SIGNAL(triggered()), this, SLOT(showEtchMenu()));

    connect(viewGroup, SIGNAL(triggered(QAction*)), this, SLOT(changeVizType(QAction*)));

    maskAction = new QAction(tr("Mask"), this);
    maskAction->setEnabled(false);
    connect(maskAction, SIGNAL(triggered()), this, SLOT(showMenuMask()));

    exitAction = new QAction(QIcon(":/images/exit.png"), tr("E&xit"), this);
    connect(exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));

    aboutQtAction = new QAction(QIcon(":/images/qt.png"), tr("About &Qt"), this);
    connect(aboutQtAction, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void
MainW::createMenus() {

    fileMenu = menuBar()->addMenu(tr("File"));
    fileMenu->addAction(newAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(exitAction);

    viewSurfaceMenu = new QMenu(tr("Surface"));
    viewBulkMenu = new QMenu(tr("Surface and bulk"));
    viewMenu = menuBar()->addMenu(tr("View"));
    viewMenu->addMenu(viewSurfaceMenu);
    viewMenu->addMenu(viewBulkMenu);

    viewBulkMenu->addAction(viewAsAts_SurfaceAndBulkAct);
    viewBulkMenu->addAction(viewAsAtsAndBonds_SurfaceAndBulkAct);

    viewSurfaceMenu->addAction(viewAsCellsSurface);
    viewSurfaceMenu->addAction(viewAsAtsAndBonds_SurfaceAct);
    viewSurfaceMenu->addAction(viewAsAts_SurfaceAct);

    helpMenu = menuBar()->addMenu(tr("Help"));
    helpMenu->addAction(aboutQtAction);
}

void
MainW::createToolBars() {
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(newAct);
    fileToolBar->addAction(saveAct);

    etchingToolBar = addToolBar(tr("Etching"));
    etchingToolBar ->addAction(maskAction);
    etchingToolBar ->addAction(etchingAction);
}

void
MainW::drawResult() {
    result->view(surfaceXYZ, surfAtoms, cell, Xsize, Ysize, Zsize,
                 z_center, z_min, SIZE_X, SIZE_Y, Vx, Vy, Vz, vizualType);
}

void
MainW::etch(int simType_, int IterCount, float *rates) {
    GRES::SimType sT = static_cast<GRES::SimType> (simType_);
    QTime t;
    //    t.start();
    int count = IterCount;
    for (int n = 0; n < count; ++n) {
        if (sT == GRES::SimType::KMC)
            perfect = surfaceXYZ.selAtom(surfAtoms, sosedi, z_min, cell, mask,
            rates);
        else if (sT == GRES::SimType::CA)
            perfect = surfaceXYZ.selAtomCA(surfAtoms, z_min, cell, mask, rates);

        if (perfect) {
            surfaceXYZ.addLayer(sosedi, SIZE_X, SIZE_Y, SIZE_Z);
            ++SIZE_Z;
            perfect = true;
        }

        z_min = surfaceXYZ.findZmin(z_min);
        surfaceXYZ.optimize(z_min);
    }
    //    qDebug() << "etch: " << t.elapsed();
    z_center = (SIZE_Z - 2 + z_min) / 2;
    drawResult();
}

void
MainW::getSettings(int hGet, int kGet, int lGet, int xS, int yS, int zS) {
    h = hGet;
    k = kGet;
    l = lGet;
    SIZE_X = xS + 4;
    SIZE_Y = yS + 4;
    SIZE_Z = zS + 2;
}

void
MainW::newDocument() {
    mask.clear();
    etchingAction->setEnabled(true);
    maskAction->setEnabled(true);
    saveAct->setEnabled(true);
    perfect = true;
    QTime t;
    z_min = 0;
    int xMax = SIZE_X;
    int yMax = SIZE_Y;
    int zMax = SIZE_Z;
    z_center = z_min + (zMax - 2 - z_min) / 2;

    cell = findCell(h, k, l, Xsize, Ysize, Zsize, Vx, Vy, Vz);
    const unsigned int NUMBER_OF_ATOMS_IN_CELL
            = static_cast<unsigned int> (cell.size());
    sosedi = cell.findSoseds(Xsize, Ysize, Zsize);

    surfaceXYZ.clear();
    surfAtoms.clear();
    surfaceXYZ.reserve(SIZE_Z);
    for (int z = 0; z < SIZE_Z; ++z) {
        Surface2D surfaceXY;
        surfaceXY.reserve(SIZE_Y);
        for (int y = 0; y < SIZE_Y; ++y) {
            Surface1D surfaceX;
            surfaceX.reserve(SIZE_X);
            for (int x = 0; x < SIZE_X; ++x) {
                vector<AtomInfo> cell;
                for (unsigned char a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a) {
                    Neighbors neighbs;
                    char numberNeighbs = 0; //Число первых соседей
                    for (int nb = 0; nb < 4; ++nb) {
                        auto &sosediANb = sosedi[a][nb];
                        if (x + sosediANb.x >= 0 && y + sosediANb.y >= 0
                            && z + sosediANb.z >= 0 &&
                            x + sosediANb.x < xMax && y + sosediANb.y < yMax
                            && z + sosediANb.z < zMax + 1) {
                            ++numberNeighbs;
                            AtomType neighb = {x + sosediANb.x, y + sosediANb.y,
                                               z + sosediANb.z, sosediANb.type, false};
                            neighbs.push_back(neighb);
                        }
                    }
                    if ((x > 1 && x < SIZE_X - 2) && (y > 1 && y < SIZE_Y - 2) && z < SIZE_Z - 2)
                        if (numberNeighbs && numberNeighbs < 4) {
                            AtomType aT = {x, y, z, a, false};
                            surfAtoms.push_back(aT);
                        }
                    AtomInfo atom = {neighbs, numberNeighbs, !numberNeighbs};
                    cell.push_back(atom);
                }
                surfaceX.push_back(cell);
            }
            surfaceXY.push_back(surfaceX);
        }
        surfaceXYZ.push_back(surfaceXY);
    }
    drawResult();
}

void
MainW::showMenuMask() {
    maskMenu = new MaskMenu(0, SIZE_X - 4, SIZE_Y - 4);
    connect(maskMenu, SIGNAL(maskChanged(std::vector<bool>)), this,
            SLOT(setMask(std::vector<bool>)));
    maskMenu->show();
}

void
MainW::setMask(const vector<bool> inMask) {
    mask = inMask;
}

void
MainW::setSettings() {
    settings->set(h, k, l, SIZE_X - 4, SIZE_Y - 4, SIZE_Z - 2);
}

void
MainW::showEtchMenu() {
    etchMenu->show();
}
