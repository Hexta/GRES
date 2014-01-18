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

#include <QWidget>

#include <vector>

using std::vector;

class QLabel;
class QPushButton;
class QImage;

class MaskMenu : public QWidget {
    Q_OBJECT
public:
    MaskMenu(QWidget *parent = 0, int xS = 10, int yS = 10);
signals:
    void maskChanged(std::vector<bool>);
private:

    QLabel* maskPreview;
    QPushButton *loadMaskButton;
    QPushButton *setMaskButton;
    QImage* maskImage;
    int xS, yS;
private slots:
    void loadMask();
    void setMask();
};
