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

#include "MainWindow.h"
#include "Render.h"
#include "Settings.h"
#include "consts.h"
#include "etchingMenu.h"
#include "maskMenu.h"

#include <QDockWidget>
#include <QActionGroup>
#include <QMenuBar>
#include <QToolBar>

struct MainWindow::Private {
    QAction *aboutQtAction;
    QAction *exitAction;
    QAction *newAct;
    QAction *saveAct;
    QAction *viewAsAtsAndBonds_SurfaceAndBulkAct;
    QAction *viewAsAtsAndBonds_SurfaceAct;

    QAction *viewAsAts_SurfaceAndBulkAct;
    QAction *viewAsAts_SurfaceAct;
    QAction *viewAsCellsSurface;
    QAction *etchingAction;
    QAction *maskAction;

    QMenu *viewSurfaceMenu;
    QMenu *viewBulkMenu;
    QActionGroup *viewGroup;

    MaskMenu *maskMenu;
    EtchingMenu* etchMenu;

    Render* result;
    Settings *settings;

    int h, k, l;
    int SIZE_X, SIZE_Y, SIZE_Z;
    int z_center, z_min;
    Cell cell;
    Surface3DPtr surfaceXYZ;
    float Xsize, Ysize, Zsize;
    AtomTypes surfAtoms;

    Coords3D Vx, Vy, Vz;
    AllNeighbors neighbours;
    Cells surfacePoints;
    GRES::VizType vizualType; //тип визуализации

    std::vector<bool> mask;

    Private(QObject* parent) :
        aboutQtAction(new QAction(QIcon(":/images/qt.png"), tr("About &Qt"), parent)),
        exitAction(new QAction(QIcon(":/images/exit.png"), tr("E&xit"), parent)),
        newAct(new QAction(QIcon(":/images/new.png"), tr("N&ew"), parent)),
        saveAct(new QAction(QIcon(":/images/filesave.png"), tr("S&ave"), parent)),
        viewAsAtsAndBonds_SurfaceAndBulkAct(new QAction(tr("Atoms and bonds"), parent)),
        viewAsAtsAndBonds_SurfaceAct(new QAction(tr("Atoms and bonds"), parent)),
        viewAsAts_SurfaceAndBulkAct(new QAction(tr("Atoms"), parent)),
        viewAsAts_SurfaceAct(new QAction(tr("Atoms"), parent)),
        viewAsCellsSurface(new QAction(tr("Surface"), parent)),
        etchingAction(new QAction(tr("Etch"), parent)),
        maskAction(new QAction(tr("Mask"), parent)),

        viewSurfaceMenu(new QMenu(tr("Surface"))),
        viewBulkMenu(new QMenu(tr("Surface and bulk"))),
        viewGroup(new QActionGroup(parent)),
        etchMenu(new EtchingMenu),
        result(new Render),
        settings(new Settings),

        h(1),
        k(0),
        l(0),
        SIZE_X(54),
        SIZE_Y(54),
        SIZE_Z(4),
        z_min(0),
        surfaceXYZ(new Surface3D),
        vizualType(GRES::VizType::CELLS_SURFACE) {

        saveAct->setEnabled(false);

        viewAsAtsAndBonds_SurfaceAndBulkAct->setCheckable(true);
        viewAsAtsAndBonds_SurfaceAct->setCheckable(true);
        viewAsAts_SurfaceAndBulkAct->setCheckable(true);
        viewAsAts_SurfaceAct->setCheckable(true);
        viewAsCellsSurface->setCheckable(true);

        viewGroup->addAction(viewAsAts_SurfaceAndBulkAct);
        viewGroup->addAction(viewAsAts_SurfaceAct);
        viewGroup->addAction(viewAsAtsAndBonds_SurfaceAct);
        viewGroup->addAction(viewAsAtsAndBonds_SurfaceAndBulkAct);
        viewGroup->addAction(viewAsCellsSurface);
        viewAsCellsSurface->setChecked(true);

        etchingAction->setEnabled(false);
        maskAction->setEnabled(false);

        surfaceXYZ->reserve(5000);
    }
};


MainWindow::MainWindow(QWidget *parent, int, char * const *) : QMainWindow(parent),
    d(new Private(this)) {

    createActions();
    createMenus();
    createToolBars();
    
    connect(d->saveAct, SIGNAL(triggered()), d->result, SLOT(saveResult()));
    setCentralWidget(d->result);
    d->result->setFocusPolicy(Qt::ClickFocus);

    QDockWidget *settingsDock = new QDockWidget(tr("Settings"), this);
    settingsDock->setAllowedAreas(Qt::LeftDockWidgetArea);
    settingsDock->setWidget(d->settings);
    addDockWidget(Qt::LeftDockWidgetArea, settingsDock);
    d->result->setFocus();
    setSettings();
    
    connect(d->settings,
            SIGNAL(settingsChanged(int, int, int, int, int, int)),
            this,
            SLOT(getSettings(int, int, int, int, int, int)));

    //     connect (result, SIGNAL (etching()), this, SLOT (etch()));
    connect(d->etchMenu, SIGNAL(startEtching(int, int, float*)), this,
            SLOT(etch(int, int, float*)));
}

void
MainWindow::changeVizType(QAction* type) {
    if (type == d->viewAsAts_SurfaceAndBulkAct)
        d->vizualType = GRES::VizType::ATOMS_SURFACE_AND_BULK;

    else if (type == d->viewAsAts_SurfaceAct)
        d->vizualType = GRES::VizType::ATOMS_SURFACE;

    else if (type == d->viewAsAtsAndBonds_SurfaceAndBulkAct)
        d->vizualType = GRES::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK;

    else if (type == d->viewAsAtsAndBonds_SurfaceAct)
        d->vizualType = GRES::VizType::ATOMS_AND_BONDS_SURFACE;

    else if (type == d->viewAsCellsSurface)
        d->vizualType = GRES::VizType::CELLS_SURFACE;

    d->result->changeVizType(d->vizualType);
}

void
MainWindow::createActions() {
    connect(d->newAct, SIGNAL(triggered()), this, SLOT(newDocument()));
    connect(d->etchingAction, SIGNAL(triggered()), this, SLOT(showEtchMenu()));
    connect(d->viewGroup, SIGNAL(triggered(QAction*)), this, SLOT(changeVizType(QAction*)));
    connect(d->maskAction, SIGNAL(triggered()), this, SLOT(showMenuMask()));
    connect(d->exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));
    connect(d->aboutQtAction, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void
MainWindow::createMenus() {
    QMenu* fileMenu = menuBar()->addMenu(tr("File"));
    fileMenu->addAction(d->newAct);
    fileMenu->addAction(d->saveAct);
    fileMenu->addAction(d->exitAction);

    QMenu* viewMenu = menuBar()->addMenu(tr("View"));
    viewMenu->addMenu(d->viewSurfaceMenu);
    viewMenu->addMenu(d->viewBulkMenu);

    d->viewBulkMenu->addAction(d->viewAsAts_SurfaceAndBulkAct);
    d->viewBulkMenu->addAction(d->viewAsAtsAndBonds_SurfaceAndBulkAct);

    d->viewSurfaceMenu->addAction(d->viewAsCellsSurface);
    d->viewSurfaceMenu->addAction(d->viewAsAtsAndBonds_SurfaceAct);
    d->viewSurfaceMenu->addAction(d->viewAsAts_SurfaceAct);

    QMenu* helpMenu = menuBar()->addMenu(tr("Help"));
    helpMenu->addAction(d->aboutQtAction);
}

void
MainWindow::createToolBars() {
    QToolBar* fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(d->newAct);
    fileToolBar->addAction(d->saveAct);

    QToolBar* etchingToolBar = addToolBar(tr("Etching"));
    etchingToolBar->addAction(d->maskAction);
    etchingToolBar->addAction(d->etchingAction);
}

void
MainWindow::drawResult() {
    d->result->view(d->surfaceXYZ, d->surfAtoms, d->cell, d->Xsize, d->Ysize,
        d->Zsize, d->z_center, d->z_min, d->SIZE_X, d->SIZE_Y, d->Vx, d->Vy,
        d->Vz, d->vizualType);
}

void
MainWindow::etch(int simType_, int IterCount, float *rates) {
    GRES::SimType sT = static_cast<GRES::SimType> (simType_);
    QTime t;
    //    t.start();
    bool perfect = false;
    for (int n = 0; n < IterCount; ++n) {
        if (sT == GRES::SimType::KMC)
            perfect = d->surfaceXYZ->selAtom(d->surfAtoms, d->neighbours,
            d->z_min, d->cell, d->mask, rates);
        else if (sT == GRES::SimType::CA)
            perfect = d->surfaceXYZ->selAtomCA(d->surfAtoms, d->z_min,
            d->cell, d->mask, rates);

        if (perfect) {
            d->surfaceXYZ->addLayer(d->neighbours, d->SIZE_X, d->SIZE_Y, d->SIZE_Z);
            ++d->SIZE_Z;
            perfect = true;
        }

        d->z_min = d->surfaceXYZ->findZmin(d->z_min);
        d->surfaceXYZ->optimize(d->z_min);
    }
    //    qDebug() << "etch: " << t.elapsed();
    d->z_center = (d->SIZE_Z - 2 + d->z_min) / 2;
    drawResult();
}

void
MainWindow::getSettings(int hGet, int kGet, int lGet, int xS, int yS, int zS) {
    d->h = hGet;
    d->k = kGet;
    d->l = lGet;
    d->SIZE_X = xS + 4;
    d->SIZE_Y = yS + 4;
    d->SIZE_Z = zS + 2;
}

void
MainWindow::newDocument() {
    d->mask.clear();
    d->etchingAction->setEnabled(true);
    d->maskAction->setEnabled(true);
    d->saveAct->setEnabled(true);

    QTime t;
    d->z_min = 0;
    int xMax = d->SIZE_X;
    int yMax = d->SIZE_Y;
    int zMax = d->SIZE_Z;
    d->z_center = d->z_min + (zMax - 2 - d->z_min) / 2;

    d->cell = findCell(d->h, d->k, d->l, d->Xsize, d->Ysize, d->Zsize, d->Vx,
        d->Vy, d->Vz);
    const unsigned int NUMBER_OF_ATOMS_IN_CELL
        = static_cast<unsigned int> (d->cell.size());
    d->neighbours = d->cell.findSoseds(d->Xsize, d->Ysize, d->Zsize);

    d->surfaceXYZ->clear();
    d->surfAtoms.clear();
    d->surfaceXYZ->reserve(d->SIZE_Z);
    for (int z = 0; z < d->SIZE_Z; ++z) {
        Surface2D surfaceXY;
        surfaceXY.reserve(d->SIZE_Y);
        for (int y = 0; y < d->SIZE_Y; ++y) {
            Surface1D surfaceX;
            surfaceX.reserve(d->SIZE_X);
            for (int x = 0; x < d->SIZE_X; ++x) {
                vector<AtomInfo> cell;
                for (unsigned char a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a) {
                    Neighbors neighbs;
                    char numberNeighbs = 0; //Число первых соседей
                    for (int nb = 0; nb < 4; ++nb) {
                        auto &sosediANb = d->neighbours[a][nb];
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
                    if ((x > 1 && x < d->SIZE_X - 2)
                        && (y > 1 && y < d->SIZE_Y - 2)
                        && z < d->SIZE_Z - 2)
                        if (numberNeighbs && numberNeighbs < 4) {
                            AtomType aT = {x, y, z, a, false};
                            d->surfAtoms.push_back(aT);
                        }
                    AtomInfo atom = {neighbs, numberNeighbs, !numberNeighbs};
                    cell.push_back(atom);
                }
                surfaceX.push_back(cell);
            }
            surfaceXY.push_back(surfaceX);
        }
        d->surfaceXYZ->push_back(surfaceXY);
    }
    drawResult();
}

void
MainWindow::showMenuMask() {
    d->maskMenu = new MaskMenu(0, d->SIZE_X - 4, d->SIZE_Y - 4);
    connect(d->maskMenu, SIGNAL(maskChanged(std::vector<bool>)), this,
            SLOT(setMask(std::vector<bool>)));
    d->maskMenu->show();
}

void
MainWindow::setMask(const vector<bool> inMask) {
    d->mask = inMask;
}

void
MainWindow::setSettings() {
    d->settings->set(d->h, d->k, d->l, d->SIZE_X - 4, d->SIZE_Y - 4,
        d->SIZE_Z - 2);
}

void
MainWindow::showEtchMenu() {
    d->etchMenu->show();
}

MainWindow::~MainWindow() {
}
