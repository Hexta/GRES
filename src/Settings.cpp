/******************************************************************************
 * Copyright (c) 2009-2016 Artur Molchanov <artur.molchanov@gmail.com>        *
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

#include "Settings.h"

#include <QWidget>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QGridLayout>
#include <QComboBox>
#include <QSpinBox>
#include <QGroupBox>
#include <QHBoxLayout>

Settings::Settings(QWidget* parent) : QWidget(parent) {
    QGroupBox* orientationBox = new QGroupBox(tr("Surface Orientation"));
    QLabel* hLabel = new QLabel(tr("H"));
    QLabel* kLabel = new QLabel(tr("K"));
    QLabel* lLabel = new QLabel(tr("L"));
    hSpinBox = new QSpinBox;
    kSpinBox = new QSpinBox;
    lSpinBox = new QSpinBox;

    QHBoxLayout* hbox = new QHBoxLayout;
    hbox->addWidget(hLabel);
    hbox->addWidget(hSpinBox);
    hbox->addWidget(kLabel);
    hbox->addWidget(kSpinBox);
    hbox->addWidget(lLabel);
    hbox->addWidget(lSpinBox);
    orientationBox->setLayout(hbox);
    QGridLayout* grid = new QGridLayout;

    QGroupBox* sizeBox = new QGroupBox(tr("Surface Size"));
    QLabel* xSizeLabel = new QLabel(tr("X"));
    QLabel* ySizeLabel = new QLabel(tr("Y"));
    QLabel* zSizeLabel = new QLabel(tr("Z"));

    xSpinBox = new QSpinBox;
    ySpinBox = new QSpinBox;
    zSpinBox = new QSpinBox;

    xSpinBox->setMaximum(999);
    ySpinBox->setMaximum(999);

    QHBoxLayout* sizeHbox = new QHBoxLayout;
    sizeHbox->addWidget(xSizeLabel);
    sizeHbox->addWidget(xSpinBox);

    sizeHbox->addWidget(ySizeLabel);
    sizeHbox->addWidget(ySpinBox);

    sizeHbox->addWidget(zSizeLabel);
    sizeHbox->addWidget(zSpinBox);

    sizeBox->setLayout(sizeHbox);

    grid->addWidget(orientationBox, 0, 0, Qt::AlignTop);
    grid->addWidget(sizeBox, 1, 0, Qt::AlignTop);
    QSpacerItem* space = new QSpacerItem(10, 10, QSizePolicy::Expanding,
        QSizePolicy::Expanding);
    grid->addItem(space, 2, 0);
    setLayout(grid);

    connect(hSpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
    connect(kSpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
    connect(lSpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
    connect(xSpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
    connect(ySpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
    connect(zSpinBox,
        SIGNAL(valueChanged(int)), this, SLOT(get()));
}

void Settings::get() {
    int h = hSpinBox->value();
    int k = kSpinBox->value();
    int l = lSpinBox->value();
    int xS = xSpinBox->value();
    int yS = ySpinBox->value();
    int zS = zSpinBox->value();
    emit settingsChanged(h, k, l, xS, yS, zS);
}

void Settings::set(int h, int k, int l, int xS, int yS, int zS) {
    hSpinBox->setValue(h);
    kSpinBox->setValue(k);
    lSpinBox->setValue(l);

    xSpinBox->setValue(xS);
    ySpinBox->setValue(yS);
    zSpinBox->setValue(zS);
}
