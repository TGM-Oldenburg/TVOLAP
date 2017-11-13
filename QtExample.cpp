/*-----------------------------------------------------------------------------*\
| Qt example for visualizing the result of the TVOLAP processing routine.        |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| GPL Release: May 2017, License see end of file                                |
\*-----------------------------------------------------------------------------*/

#include "QtMainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow mainWin;
    mainWin.setMaximumSize(1200,200);
    mainWin.resize(1200,200);

    mainWin.show();

    return a.exec();
}

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
