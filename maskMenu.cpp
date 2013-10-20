#include "maskMenu.h"
#include <QtDebug>
#include <QLabel>
#include <QGridLayout>
#include <QPushButton>
#include <QFileDialog>
#include <QMessageBox>
#include <QImage>
#include <QPixmap>

MaskMenu::MaskMenu(QWidget *parent, int xS_, int yS_) : QWidget(parent) {
    maskPreview = new QLabel();
    maskPreview->setMinimumSize(200, 200);

    loadMaskButton = new QPushButton(tr("Load Mask"));
    connect(loadMaskButton, SIGNAL(clicked()), this, SLOT(loadMask()));

    setMaskButton = new QPushButton(tr("Set Mask"));
    connect(setMaskButton, SIGNAL(clicked()), this, SLOT(setMask()));

    QGridLayout *grid = new QGridLayout;
    grid ->addWidget(maskPreview, 0, 0);
    grid ->addWidget(loadMaskButton, 0, 1);
    grid ->addWidget(setMaskButton, 1, 1);
    setLayout(grid);
    xS = xS_;
    yS = yS_;
}

void
MaskMenu::loadMask() {
    QString selFilter = "";
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open file with mask"),
                                                    "./masks/", tr("Images (*.bmp *.gif *.png *.tiff *.jpg *.xpm)"),
                                                    &selFilter);
    if (!fileName.isNull())//проверяем не отменил ли юзер выбор файла
    {
        maskImage = new QImage(fileName);
        if (maskImage->isNull()) {
            QMessageBox::information(this, tr("GRES"),
                                     tr("Cannot load %1.").arg(fileName));
            return;
        }
        *maskImage = maskImage->convertToFormat(QImage::Format_Mono);
        QPixmap maskPixmap = QPixmap::fromImage(*maskImage);
        maskPixmap = maskPixmap.scaled(200, 200);
        maskPreview->setPixmap(maskPixmap);
    }
}

void
MaskMenu::setMask() {
    vector<bool> mask;
    *maskImage = maskImage->scaled(xS, yS);
    for (int i = 0; i < yS; ++i)
        for (int j = 0; j < xS; ++j)
            mask.push_back(!qGray(maskImage->pixel(j, i)));
    emit maskChanged(&mask);
    this->close();
}
