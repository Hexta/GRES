/******************************************************************************
 * Copyright (c) 2009-2016 Artur Molchanov <artur.molchanov@gmail.com>     *
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
#include "EtchingMenu.h"
#include "MaskMenu.h"

#include "Surface3D.h"
#include "Cell.h"

#include "SimulationType.h"

#include <QDockWidget>
#include <QActionGroup>
#include <QMenuBar>
#include <QToolBar>

struct MainWindow::Private {
    QAction* aboutQtAction;
    QAction* exitAction;
    QAction* newAct;
    QAction* saveAct;
    QAction* viewAsAtsAndBonds_SurfaceAndBulkAct;
    QAction* viewAsAtsAndBonds_SurfaceAct;

    QAction* viewAsAts_SurfaceAndBulkAct;
    QAction* viewAsAts_SurfaceAct;
    QAction* viewAsCellsSurface;
    QAction* etchingAction;
    QAction* maskAction;

    QMenu* viewSurfaceMenu;
    QMenu* viewBulkMenu;
    QActionGroup* viewGroup;

    MaskMenu* maskMenu;
    EtchingMenu* etchMenu;

    Render* result;
    Settings* settings;

    int h, k, l;
    int SIZE_X, SIZE_Y, SIZE_Z;
    int z_center, z_min;
    Cell cell;
    Surface3DPtr surfaceXYZ;
    float Xsize, Ysize, Zsize;

    Coords3D Vx, Vy, Vz;
    AllNeighbors neighbors;
    Cells surfacePoints;
    Render::VizType vizualizationType;

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
        surfaceXYZ(nullptr),
        vizualizationType(Render::VizType::CELLS_SURFACE) {
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
    }

    void newDocument() {
        mask.clear();
        etchingAction->setEnabled(true);
        maskAction->setEnabled(true);
        saveAct->setEnabled(true);

        QTime t;
        z_min = 0;
        int const xMax = SIZE_X;
        int const yMax = SIZE_Y;
        int const zMax = SIZE_Z;
        z_center = z_min + (zMax - 2 - z_min) / 2;

        cell = Cell(h, k, l);
        Xsize = cell.getXSize();
        Ysize = cell.getYSize();
        Zsize = cell.getZSize();

        Vx = cell.getVx();
        Vy = cell.getVy();
        Vz = cell.getVz();

        cell.optimize();

        auto const numberOfAtomInCell =
            static_cast<unsigned int>(cell.size());
        neighbors = cell.findNeighbors(Xsize, Ysize, Zsize);

        surfaceXYZ.reset(new Surface3D(cell));
        surfaceXYZ->clear();
        surfaceXYZ->reserve(SIZE_Z);
        for (int z = 0; z < SIZE_Z; ++z) {
            Surface2D surfaceXY;
            surfaceXY.reserve(SIZE_Y);
            for (int y = 0; y < SIZE_Y; ++y) {
                Surface1D surfaceX;
                surfaceX.reserve(SIZE_X);
                for (int x = 0; x < SIZE_X; ++x) {
                    Cell cell;
                    for (unsigned char a = 0; a < numberOfAtomInCell; ++a) {
                        Neighbors neighbs;
                        char numberNeighbs = 0; // number of the first neighbors
                        for (int nb = 0; nb < 4; ++nb) {
                            auto& neighborsANb = neighbors[a][nb];
                            if (x + neighborsANb.x >= 0
                                && y + neighborsANb.y >= 0
                                && z + neighborsANb.z >= 0
                                && x + neighborsANb.x < xMax
                                && y + neighborsANb.y < yMax
                                && z + neighborsANb.z < zMax + 1) {
                                ++numberNeighbs;
                                AtomType neighb = {x + neighborsANb.x, y + neighborsANb.y,
                                                   z + neighborsANb.z, neighborsANb.type,
                                                   false};
                                neighbs.push_back(neighb);
                            }
                        }

                        AtomInfo atom;
                        atom.neighbors = neighbs;
                        atom.firstNeighborsCount = numberNeighbs;
                        atom.deleted = numberNeighbs == 0;
                        // TODO: atom.type =
                        cell.addAtom(atom);
                    }
                    surfaceX.push_back(cell);
                }
                surfaceXY.push_back(surfaceX);
            }
            surfaceXYZ->push_back(surfaceXY);
        }
        surfaceXYZ->rebuildSurfaceAtoms();
        drawResult();
    }

    void drawResult() {
        result->view(surfaceXYZ, cell, Xsize, Ysize,
            Zsize, z_center, z_min, SIZE_X, SIZE_Y, Vx, Vy,
            Vz, vizualizationType);
    }

    void etch(int simType_, int IterCount, float* rates) {
        GRES::SimType sT = static_cast<GRES::SimType>(simType_);
        QTime t;
        // t.start();
        bool perfect = false;
        auto surface = surfaceXYZ;

        for (int n = 0; n < IterCount; ++n) {
            if (sT == GRES::SimType::KMC)
                perfect = surface->deleteRandomAtomKmc(neighbors, z_min, mask, rates);
            else if (sT == GRES::SimType::CA)
                perfect = surface->deleteRandomAtomCa(z_min, mask, rates);

            if (perfect) {
                surface->addLayer(neighbors, SIZE_X, SIZE_Y, SIZE_Z);
                ++SIZE_Z;
                perfect = true;
            }

            z_min = surface->findZmin(z_min);
            surface->optimize(z_min);
        }
        // qDebug() << "etch: " << t.elapsed();
        z_center = (SIZE_Z - 2 + z_min) / 2;
        drawResult();
    }

    void setSettings() {
        settings->set(h, k, l, SIZE_X - 4, SIZE_Y - 4, SIZE_Z - 2);
    }

    void getSettings(int hGet, int kGet, int lGet, int xS, int yS, int zS) {
        h = hGet;
        k = kGet;
        l = lGet;
        SIZE_X = xS + 4;
        SIZE_Y = yS + 4;
        SIZE_Z = zS + 2;
    }
};

MainWindow::MainWindow(QWidget* parent, int, char* const*) :
    QMainWindow(parent),
    d(new Private(this)) {
    createActions();
    createMenus();
    createToolBars();

    connect(d->saveAct, SIGNAL(triggered()), d->result, SLOT(saveResult()));
    setCentralWidget(d->result);
    d->result->setFocusPolicy(Qt::ClickFocus);

    QDockWidget* settingsDock = new QDockWidget(tr("Settings"), this);
    settingsDock->setAllowedAreas(Qt::LeftDockWidgetArea);
    settingsDock->setWidget(d->settings);
    addDockWidget(Qt::LeftDockWidgetArea, settingsDock);
    d->result->setFocus();
    d->setSettings();

    connect(d->settings,
        SIGNAL(settingsChanged(int, int, int, int, int, int)),
        this,
        SLOT(getSettings(int, int, int, int, int, int)));

    // connect (result, SIGNAL (etching()), this, SLOT (etch()));
    connect(d->etchMenu, SIGNAL(startEtching(int, int, float*)), this,
        SLOT(etch(int, int, float*)));
}

void MainWindow::changeVizType(QAction* type) {
    if (type == d->viewAsAts_SurfaceAndBulkAct)
        d->vizualizationType = Render::VizType::ATOMS_SURFACE_AND_BULK;

    else if (type == d->viewAsAts_SurfaceAct)
        d->vizualizationType = Render::VizType::ATOMS_SURFACE;

    else if (type == d->viewAsAtsAndBonds_SurfaceAndBulkAct)
        d->vizualizationType = Render::VizType::ATOMS_AND_BONDS_SURFACE_AND_BULK;

    else if (type == d->viewAsAtsAndBonds_SurfaceAct)
        d->vizualizationType = Render::VizType::ATOMS_AND_BONDS_SURFACE;

    else if (type == d->viewAsCellsSurface)
        d->vizualizationType = Render::VizType::CELLS_SURFACE;

    d->result->changeVizType(d->vizualizationType);
}

void MainWindow::createActions() {
    connect(d->newAct, SIGNAL(triggered()), this, SLOT(newDocument()));
    connect(d->etchingAction, SIGNAL(triggered()), this, SLOT(showEtchMenu()));
    connect(d->viewGroup, SIGNAL(triggered(QAction*)), this, SLOT(changeVizType(QAction*)));
    connect(d->maskAction, SIGNAL(triggered()), this, SLOT(showMenuMask()));
    connect(d->exitAction, SIGNAL(triggered()), qApp, SLOT(quit()));
    connect(d->aboutQtAction, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void MainWindow::createMenus() {
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

void MainWindow::createToolBars() {
    QToolBar* fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(d->newAct);
    fileToolBar->addAction(d->saveAct);

    QToolBar* etchingToolBar = addToolBar(tr("Etching"));
    etchingToolBar->addAction(d->maskAction);
    etchingToolBar->addAction(d->etchingAction);
}

void MainWindow::etch(int simType_, int IterCount, float* rates) {
    d->etch(simType_, IterCount, rates);
}

void MainWindow::getSettings(int h, int k, int l, int xS, int yS, int zS) {
    d->getSettings(h, k, l, xS, yS, zS);
}

void MainWindow::newDocument() {
    d->newDocument();
}

void MainWindow::showMenuMask() {
    d->maskMenu = new MaskMenu(0, d->SIZE_X - 4, d->SIZE_Y - 4);
    connect(d->maskMenu, SIGNAL(maskChanged(std::vector<bool>)), this,
        SLOT(setMask(std::vector<bool>)));
    d->maskMenu->show();
}

void MainWindow::setMask(const std::vector<bool> inMask) {
    d->mask = inMask;
}

void MainWindow::showEtchMenu() {
    d->etchMenu->show();
}

MainWindow::~MainWindow() {
}
