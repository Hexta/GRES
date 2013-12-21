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

class QGroupBox;
class QSpinBox;

class Settings : public QWidget {
    Q_OBJECT
public:
    Settings(QWidget *parent = 0);
public slots:
    void set(int h, int k, int l, int xS, int yS, int zS);
signals:
    void settingsChanged(int, int, int, int, int, int);
private:
    QSpinBox* hSpinBox;
    QSpinBox* kSpinBox;
    QSpinBox* lSpinBox;
    QSpinBox* xSpinBox;
    QSpinBox* ySpinBox;
    QSpinBox* zSpinBox;

private slots:
    void get();

};

