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

#ifndef SELECT_ATOM_MENU
#    define SELECT_ATOM_MENU

#    include <QWidget>

class QLineEdit;
class QLabel;

class SelectAtomMenu : public QWidget
{
    Q_OBJECT
public:
    SelectAtomMenu ( QWidget *parent = 0 );
    QSize sizeHint( ) const;
public slots:
    void setInfo( int x, int y, int z, int type, int fNbCount );
private:
    QLabel *atomPosition;
    QLabel *firstNeighboursCount;
} ;
#endif