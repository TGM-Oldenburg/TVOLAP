/*-----------------------------------------------------------------------------*\
| Header of QtMainwindow.cpp.                                                   |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| GPL Release: May 2017, License see end of file                                |
\*-----------------------------------------------------------------------------*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    QCustomPlot *inPlot, *irPlot, *outPlot;
    QVector<double> xInPlot, xOutPlot, xIrPlot, yIrPlot, yInPlot, yOutPlot;
};

#endif // MAINWINDOW_H

/*------------------------------License---------------------------------------*\
| Copyright (c) 2012-2017 Hagen Jaeger                                         |
|                                                                              |
| This program is free software: you can redistribute it and/or modify         |
| it under the terms of the GNU General Public License as published by         |
| the Free Software Foundation, either version 3 of the License, or            |
| (at your option) any later version.                                          |
|                                                                              |
| This program is distributed in the hope that it will be useful,              |
| but WITHOUT ANY WARRANTY; without even the implied warranty of               |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                |
| GNU General Public License for more details.                                 |
|                                                                              |
| You should have received a copy of the GNU General Public License            |
| along with this program. If not, see <http://www.gnu.org/licenses/>.         |
\*----------------------------------------------------------------------------*/
