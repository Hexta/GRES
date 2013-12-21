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

#ifndef ETCHING_MENU
#define ETCHING_MENU
#include <QWidget>

class QLineEdit;
class QPushButton;
class QComboBox;
class EtchingMenu : public QWidget
{
	Q_OBJECT
	public:
		EtchingMenu (QWidget *parent=0);
	signals:
		void startEtching(int, int iterations, float* rates);
		
	private:
		QComboBox *simTypeComboBox;
		QLineEdit *iterationsLine;
		QLineEdit *f1NbLine;
		QLineEdit *f2NbLine;
		QLineEdit *f3NbLine;
		QLineEdit *f4NbLine;
		QLineEdit *simTypeLine;
		QPushButton *etchButton;
		float removRates[3];
	private slots:
		void startEtching();
};
#endif
