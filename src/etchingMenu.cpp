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

#include <QLabel>
#include <QLineEdit>
#include <QGridLayout>
#include <QPushButton>
#include <QtDebug>
#include <QGroupBox>
#include <QComboBox>
#include "etchingMenu.h"

EtchingMenu::EtchingMenu(QWidget *parent) : QWidget(parent) {
    const float P1 = 1.0;
    const float P2 = 0.030;
    const float P3 = 0.000010;
    QLabel *iterationsLabel = new QLabel(tr("Number of iterations"));
    iterationsLine = new QLineEdit;
    iterationsLine->setText(QString::number(2000));
    QGridLayout *iterationsGrid = new QGridLayout;

    etchButton = new QPushButton(tr("Start etching"));

    QGroupBox *probabilityBox = new QGroupBox(tr("Removal rates"));
    QGridLayout *ratesGrid = new QGridLayout;
    QLabel *fNbCountLabel = new QLabel(tr("Count of first neighbours"));

    QLabel *f1NbLabel = new QLabel(tr("1"));
    QLabel *f2NbLabel = new QLabel(tr("2"));
    QLabel *f3NbLabel = new QLabel(tr("3"));

    QLabel *ratesLabel = new QLabel(tr("Removal rates"));
    f1NbLine = new QLineEdit;
    f2NbLine = new QLineEdit;
    f3NbLine = new QLineEdit;
    f1NbLine->setText(QString::number(P1));
    f2NbLine->setText(QString::number(P2));
    f3NbLine->setText(QString::number(P3));

    ratesGrid -> addWidget(fNbCountLabel, 0, 0);
    ratesGrid -> addWidget(ratesLabel, 0, 1);
    ratesGrid -> addWidget(f1NbLabel, 1, 0);
    ratesGrid -> addWidget(f2NbLabel, 2, 0);
    ratesGrid -> addWidget(f3NbLabel, 3, 0);
    ratesGrid -> addWidget(f1NbLine, 1, 1);
    ratesGrid -> addWidget(f2NbLine, 2, 1);
    ratesGrid -> addWidget(f3NbLine, 3, 1);
    probabilityBox->setLayout(ratesGrid);

    iterationsGrid -> addWidget(probabilityBox, 0, 0, 1, 2);


    QLabel *simulationTypeLabel = new QLabel(tr("Type of modeling"));
    simTypeLine = new QLineEdit;
    simTypeComboBox = new QComboBox;
    iterationsGrid -> addWidget(simulationTypeLabel, 1, 0);
    iterationsGrid -> addWidget(simTypeComboBox, 1, 1);
    simTypeComboBox->insertItems(0, QStringList()
                                 << tr("CA")
                                 << tr("KMC"));

    iterationsGrid -> addWidget(iterationsLabel, 2, 0);
    iterationsGrid -> addWidget(iterationsLine, 2, 1);
    iterationsGrid -> addWidget(etchButton, 3, 1);

    QSpacerItem *iterationSpace = new QSpacerItem(1, 1, QSizePolicy::Expanding,
                                                  QSizePolicy::Minimum);
    iterationsGrid->addItem(iterationSpace, 4, 0);

    setLayout(iterationsGrid);
    iterationsLine->setFocus(Qt::MouseFocusReason);

    connect(etchButton, SIGNAL(clicked()), this, SLOT(startEtching()));
}

void EtchingMenu::startEtching() {
    int simType = simTypeComboBox->currentIndex();
    int iterCount = iterationsLine->text().toInt();
    removRates[0] = f1NbLine->text().toFloat();
    removRates[1] = f2NbLine->text().toFloat();
    removRates[2] = f3NbLine->text().toFloat();
    emit startEtching(simType, iterCount, &removRates[0]);
    this -> close();
}
