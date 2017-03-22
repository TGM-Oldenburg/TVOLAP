#include "mainwindow.h"
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
