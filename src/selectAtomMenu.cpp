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

#include "selectAtomMenu.h"
#include <QLabel>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLineEdit>
#include <QGroupBox>

SelectAtomMenu::SelectAtomMenu(QWidget *parent) : QWidget(parent) {
    QLabel *firstNeighborsCountLabel = new QLabel(tr("first"));

    firstNeighborsCount = new QLabel();
    firstNeighborsCount -> setTextInteractionFlags(Qt::TextSelectableByMouse);

    atomPosition = new QLabel();

    QGroupBox *atomPositionBox = new QGroupBox(tr("Atom`s position"));

    QGridLayout *atomPositionGrid = new QGridLayout;

    atomPositionGrid ->addWidget(atomPosition, 0, 0);

    QSpacerItem *positionSpace = new QSpacerItem(1, 1, QSizePolicy::Expanding,
                                                 QSizePolicy::Minimum);
    atomPositionGrid->addItem(positionSpace, 0, 1);
    atomPositionBox->setLayout(atomPositionGrid);

    QGridLayout *grid = new QGridLayout;
    grid->addWidget(atomPositionBox, 0, 0);

    QGroupBox *neighborsBox = new QGroupBox(tr("Atom`s neighbors"));
    QGridLayout *neighborGrid = new QGridLayout;
    neighborGrid -> addWidget(firstNeighborsCountLabel, 0, 0);
    neighborGrid -> addWidget(firstNeighborsCount, 0, 1);
    QSpacerItem *neighbSpace = new QSpacerItem(1, 1, QSizePolicy::Expanding,
                                               QSizePolicy::Minimum);
    neighborGrid -> addItem(neighbSpace, 0, 2);
    neighborsBox->setLayout(neighborGrid);

    grid->addWidget(neighborsBox, 1, 0);

    QSpacerItem *space = new QSpacerItem(10, 10, QSizePolicy::Expanding,
                                         QSizePolicy::Expanding);
    grid->addItem(space, 2, 0);

    setLayout(grid);
}

void
SelectAtomMenu::setInfo(int x, int y, int z, int type, int fNbCount) {
    QString position;
    position = "(" + QString::number(x) + "," + QString::number(y)
            + "," + QString::number(z) +
            "," + QString::number(type) + ")";
    atomPosition -> setText(position);
    firstNeighborsCount->setText(QString::number(fNbCount));
}

QSize
SelectAtomMenu::sizeHint() const {
    return QSize(100, 300);
}
