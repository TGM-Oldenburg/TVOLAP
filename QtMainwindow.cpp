/*-----------------------------------------------------------------------------*\
| QtMainwindow which calculates and visualizes the result of the TVOLAP          |
| processing routine.                                                           |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| GPL Release: May 2017, License see end of file                                |
\*-----------------------------------------------------------------------------*/

#include <stdint.h>
#include <vector>
#include <cmath>
#include "QtMainwindow.h"
#include "TVOLAP.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    //----------------------------TestSignalGeneration----------------------------------
    uint32_t sampFreq = 48000;
    double timeDurAudio = 3.0;
    double timeDurIR = 0.25;
    double tmpSum = 0.0001, alpha, minVal, maxVal, tmpVal;
    uint32_t squareWaveFreq = 3;
    uint32_t oscillatorFreq = 10;
    uint32_t freqIncrement = 10;
    uint32_t numChansAudio = 3;
    uint32_t numSampsIRPerChan = (double) sampFreq*timeDurIR;
    uint32_t numIR = 4;
    uint32_t numChansIR = 3;
    uint32_t numSampsAudioPerChan = (double) sampFreq*timeDurAudio;
    uint32_t blockLen = 512;
    uint32_t numBlocks = numSampsAudioPerChan/blockLen;
    uint32_t numBlocksPerIR = numBlocks/numIR;

    std::vector<double> testSignal(numChansAudio*numSampsAudioPerChan, 0.0);
    std::vector< std::vector< std::vector<double> > > testIR;
    std::vector<double> interleavedIR(numIR*numSampsIRPerChan*numChansIR, 0.0);

    testIR.resize(numIR);
    for (unsigned int i = 0; i < numIR; i++)
    {
    	testIR.at(i).resize(numChansIR);
    	for (unsigned int j = 0; j < numChansIR; j++)
    		testIR.at(i).at(j).resize(numSampsIRPerChan, 0.0);
    }

    //test input is a square wave with iFreq Hz
    for (unsigned int i=0; i<numSampsAudioPerChan*numChansAudio; i+=numChansAudio)
    {
        tmpVal = 0.7071*(((sin(2.0*M_PI*squareWaveFreq*((double)i/numChansAudio/sampFreq)))>0) ? (1):(-1));
        for (unsigned int j=0; j<numChansAudio; j++)
            testSignal.at(i+j) = tmpVal;
    }

    //test impulse response is a damped oscillator with same frequency as test signal
    for (unsigned int i=0; i<numIR; i++)
    {
		for (unsigned int j=0; j<numChansIR; j++)
			testIR.at(i).at(j).at(numSampsIRPerChan/2) = 1.0;

		alpha = 1.0 - 1.0/(timeDurIR*0.1*sampFreq);
		for (unsigned int k=numSampsIRPerChan/2+1; k<numSampsIRPerChan; k++)
		{
			for (unsigned int j=0; j<numChansIR; j++)
				testIR.at(i).at(j).at(k) = alpha*testIR.at(i).at(j).at(k-1);
		}
		for (unsigned int k=numSampsIRPerChan/2; k<numSampsIRPerChan; k++)
		{
			tmpVal = -sin(2.0*M_PI*oscillatorFreq*((double)(k-numSampsIRPerChan/2)/sampFreq));
			for (unsigned int j=0; j<numChansIR; j++)
				testIR.at(i).at(j).at(k) *= tmpVal;
		}
		oscillatorFreq += freqIncrement;
    }

    for (unsigned int i=0; i<numIR; i++)
    {
    	for (unsigned int j=0; j<numChansIR; j++)
    	{
			tmpSum = 0.0;
			for (unsigned int k=0; k<numSampsIRPerChan; k++)
			{
				tmpVal = testIR.at(i).at(j).at(k);
				tmpSum += (tmpVal>0) ? (tmpVal) : (-tmpVal);
			}
			for (unsigned int k=0; k<numSampsIRPerChan; k++)
			{
				testIR.at(i).at(j).at(k)/= tmpSum;
				interleavedIR.at(k+(i*numChansIR+j)*numSampsIRPerChan) = testIR.at(i).at(j).at(k);
			}
    	}
    }

    //---------------------------------TVOLAP init----------------------------------------

    TVOLAP TVOLAPInst(interleavedIR, numIR, numSampsIRPerChan, numChansIR, blockLen, numChansAudio);

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

    for (unsigned int i=0, j=0, k=0; i<numBlocks; i++, j++)
    {
    	if (j>=numBlocksPerIR)
    	{
    		j = 0;
    		TVOLAPInst.setIR(++k);
    	}

        TVOLAPInst.process(&testSignal[i*numChansAudio*blockLen]);
    }

    for (unsigned int i=numBlocks*blockLen*numChansAudio; i<numSampsAudioPerChan*numChansAudio; i++)
        testSignal[i] = 0.0;

    irPlot = new QCustomPlot(this);
    irPlot->move(400,0);
    irPlot->resize(QSize(400,200));
    irPlot->addGraph();
    irPlot->graph(0)->setPen(QColor(230,50,50));
    irPlot->graph(0)->setBrush(plotBrush);
    irPlot->xAxis->setRange(1,numSampsIRPerChan*numChansIR*numIR);
    minVal = *std::min_element(interleavedIR.begin(),interleavedIR.end());
    maxVal = *std::max_element(interleavedIR.begin(),interleavedIR.end());
    irPlot->yAxis->setRange(1.2*minVal, 1.2*maxVal);
    irPlot->setBackground(this->palette().background().color());
    xIrPlot.resize(numSampsIRPerChan*numChansIR*numIR);
    yIrPlot.resize(numSampsIRPerChan*numChansIR*numIR);
    for (unsigned int cnt=0; cnt<numSampsIRPerChan*numChansIR*numIR; cnt++)
    {
        xIrPlot[cnt] = cnt;
        yIrPlot[cnt] = interleavedIR[cnt];
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
