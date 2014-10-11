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

#include "maskMenu.h"

#include <QtDebug>
#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QImage>
#include <QPixmap>

struct MaskMenu::Private {
    Private() :
        maskPreview(new QLabel),
        loadMaskButton(new QPushButton(tr("Load Mask"))),
        setMaskButton(new QPushButton(tr("Set Mask"))),
        maskImage(new QImage) {
    }

    QLabel* maskPreview;
    QPushButton* loadMaskButton;
    QPushButton* setMaskButton;
    QImage* maskImage;
    int xS, yS;
};

MaskMenu::MaskMenu(QWidget* parent, int xS_, int yS_) :
    QWidget(parent),
    d(new Private) {
    d->maskPreview->setMinimumSize(200, 200);

    connect(d->loadMaskButton, SIGNAL(clicked()), this, SLOT(loadMask()));
    connect(d->setMaskButton, SIGNAL(clicked()), this, SLOT(setMask()));

    QGridLayout* grid = new QGridLayout;
    grid->addWidget(d->maskPreview, 0, 0);
    grid->addWidget(d->loadMaskButton, 0, 1);
    grid->addWidget(d->setMaskButton, 1, 1);
    setLayout(grid);

    d->xS = xS_;
    d->yS = yS_;
}

void MaskMenu::loadMask() {
    QString selFilter = "";
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open file with mask"),
        "./masks/", tr("Images (*.bmp *.gif *.png *.tiff *.jpg *.xpm)"),
        &selFilter);

    // check that user canceled file selection
    if (!fileName.isNull()) {
        d->maskImage = new QImage(fileName);
        if (d->maskImage->isNull()) {
            QMessageBox::information(this, tr("GRES"),
                tr("Cannot load %1.").arg(fileName));
            return;
        }

        *d->maskImage = d->maskImage->convertToFormat(QImage::Format_Mono);
        QPixmap maskPixmap = QPixmap::fromImage(*d->maskImage);
        maskPixmap = maskPixmap.scaled(200, 200);
        d->maskPreview->setPixmap(maskPixmap);
    }
}

void MaskMenu::setMask() {
    std::vector<bool> mask;
    auto const& xSize = d->xS;
    auto const& ySize = d->yS;

    mask.reserve(ySize * xSize);
    *d->maskImage = d->maskImage->scaled(xSize, ySize);
    for (int i = 0; i < ySize; ++i)
        for (int j = 0; j < xSize; ++j)
            mask.push_back(!qGray(d->maskImage->pixel(j, i)));
    emit maskChanged(mask);
    this->close();
}

MaskMenu::~MaskMenu() {
}
