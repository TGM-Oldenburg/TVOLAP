#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include "WOLAP.h"

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
