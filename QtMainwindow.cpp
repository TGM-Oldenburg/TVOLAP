/*-----------------------------------------------------------------------------*\
| QtMainwindow which calculates and visualizes the result of the WOLAP          |
| processing routine.                                                           |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| GPL Release: May 2017, License see end of file                                |
\*-----------------------------------------------------------------------------*/

#include "QtMainwindow.h"
#include <stdint.h>
#include <vector>
#include <WOLAP.h>
#include <cmath>

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    //----------------------------TestSignalGeneration----------------------------------
    uint32_t sampFreq = 48000;
    double timeDurAudio = 2.0;
    double timeDurIR = 0.2;
    double tmpSum = 0.0001, alpha, minVal, maxVal, tmpVal;
    uint32_t squareWaveFreq = 3;
    uint32_t oscillatorFreq = 30;
    uint32_t numChansAudio = 2;
    uint32_t lenIR = (double) sampFreq*timeDurIR;
    uint32_t numChansIR = 2;
    uint32_t numSampsAudioPerChan = (double) sampFreq*timeDurAudio;
    uint32_t numSampsAudio = numSampsAudioPerChan*numChansAudio;
    uint32_t numSampsIR = numChansIR*lenIR;
    uint32_t blockLen = 512;
    uint32_t numBlocks = numSampsAudioPerChan/blockLen;

    std::vector<double> testSignal(numSampsAudio, 0.0);
    std::vector<double> testIR(numSampsIR, 1.0);

    //test input is a square wave with iFreq Hz
    for (unsigned int i=0; i<numSampsAudio; i+=numChansAudio)
    {
        tmpVal = 0.5*(((sin(2.0*M_PI*squareWaveFreq*((double)i/numChansAudio/sampFreq)))>0) ? (1):(-1));
        for (unsigned int k=0; k<numChansAudio; k++)
            testSignal.at(i+k) = tmpVal;
    }

    //test impulse response is a damped oscillator with same frequency as test signal
    alpha = 1.0 - 1.0/(timeDurIR*0.1*sampFreq);
    for (unsigned int i=numChansIR; i<numSampsIR; i++)
        testIR.at(i) = alpha*testIR[i-numChansIR];

    for (unsigned int i=0; i<numSampsIR/2; i++) //test input is a delayed delta impulse
        testIR.at(i) = 0.0;

    for (unsigned int i=numSampsIR/2; i<numSampsIR; i+=numChansIR) //test input is a delayed delta impulse
    {
        tmpVal = -sin(2.0*M_PI*oscillatorFreq*((double)(i-numSampsIR/2)/numChansIR/sampFreq));
        for (unsigned int k=0; k<numChansIR; k++)
            testIR.at(i+k) *= tmpVal;
    }

    for (unsigned int i=0; i<numSampsIR; i++)
        tmpSum += (testIR[i]>0) ? (testIR[i]) : (-testIR[i]);

    for (unsigned int i=0; i<numSampsIR; i++)
        testIR[i]/=tmpSum;

    //---------------------------------WOLAP init----------------------------------------

    WOLAP wolapInst(testIR, lenIR, numChansIR, blockLen, numChansAudio);

    //--------------------------------Plotting stuff-------------------------------------

    QBrush plotBrush(QColor(150,150,150,150));
    uint32_t plotChannel = 0;

    inPlot = new QCustomPlot(this);
    inPlot->move(0,0);
    inPlot->resize(QSize(400,200));
    inPlot->addGraph();
    inPlot->graph(0)->setPen(QColor(230,50,50));
    inPlot->graph(0)->setBrush(plotBrush);
    inPlot->xAxis->setRange(1,numBlocks*blockLen);
    minVal = *std::min_element(testSignal.begin(),testSignal.end());
    maxVal = *std::max_element(testSignal.begin(),testSignal.end());
    inPlot->yAxis->setRange(1.2*minVal,1.2*maxVal);
    inPlot->setBackground(this->palette().background().color());
    xInPlot.resize(numSampsAudioPerChan);
    yInPlot.resize(numSampsAudioPerChan);
    for (unsigned int cnt=0; cnt<numBlocks*blockLen; cnt++)
    {
        xInPlot[cnt] = cnt;
        yInPlot[cnt] = testSignal[cnt*numChansAudio+plotChannel];
    }
    inPlot->graph(0)->setData(xInPlot, yInPlot, true);

    for (unsigned int i=0; i<numBlocks; i++)
        wolapInst.process(&testSignal[i*numChansAudio*blockLen]);

    for (unsigned int i=numBlocks*blockLen*numChansAudio; i<numSampsAudio; i++)
        testSignal[i] = 0.0;

    irPlot = new QCustomPlot(this);
    irPlot->move(400,0);
    irPlot->resize(QSize(400,200));
    irPlot->addGraph();
    irPlot->graph(0)->setPen(QColor(230,50,50));
    irPlot->graph(0)->setBrush(plotBrush);
    irPlot->xAxis->setRange(1,lenIR);
    minVal = *std::min_element(testIR.begin(),testIR.end());
    maxVal = *std::max_element(testIR.begin(),testIR.end());
    irPlot->yAxis->setRange(1.2*minVal, 1.2*maxVal);
    irPlot->setBackground(this->palette().background().color());
    xIrPlot.resize(lenIR);
    yIrPlot.resize(lenIR);
    for (unsigned int cnt=0; cnt<lenIR; cnt++)
    {
        xIrPlot[cnt] = cnt;
        yIrPlot[cnt] = testIR[cnt*numChansIR+plotChannel];
    }
    irPlot->graph(0)->setData(xIrPlot, yIrPlot, true);

    outPlot = new QCustomPlot(this);
    outPlot->move(800,0);
    outPlot->resize(QSize(400,200));
    outPlot->addGraph();
    outPlot->graph(0)->setPen(QColor(230,50,50));
    outPlot->graph(0)->setBrush(plotBrush);
    outPlot->xAxis->setRange(1,numBlocks*blockLen);
    minVal = *std::min_element(testSignal.begin(),testSignal.end());
    maxVal = *std::max_element(testSignal.begin(),testSignal.end());
    outPlot->yAxis->setRange(1.2*minVal,1.2*maxVal);
    outPlot->setBackground(this->palette().background().color());
    xOutPlot.resize(numSampsAudioPerChan);
    yOutPlot.resize(numSampsAudioPerChan);
    for (unsigned int cnt=0; cnt<numBlocks*blockLen; cnt++)
    {
        xOutPlot[cnt] = cnt;
        yOutPlot[cnt] = testSignal[cnt*numChansAudio+plotChannel];
    }
    outPlot->graph(0)->setData(xOutPlot, yOutPlot, true);
}

MainWindow::~MainWindow()
{
    delete inPlot;
    delete irPlot;
    delete outPlot;
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
