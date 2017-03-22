#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    //----------------------------TestSignalGeneration----------------------------------
    uint32_t iFs = 48000;
    double fTimeDurAudio = 2.0;
    double fTimeDurIR = 0.2;
    double fTmpSum=0.0001, fAlpha, fMin, fMax, fTmp;
    uint32_t iFreqAudio = 3;
    uint32_t iFreqIR = 30;
    uint32_t iNumChansAudio = 2;
    uint32_t iLenIR = (double) iFs*fTimeDurIR;
    uint32_t iNumChansIR = 2;
    uint32_t iNumSampsAudioPerChan = (double) iFs*fTimeDurAudio;
    uint32_t iNumSampsAudio = iNumSampsAudioPerChan*iNumChansAudio;
    uint32_t iNumSampsIR = iNumChansIR*iLenIR;
    uint32_t iBlockLen = 512;
    uint32_t iNumBlocks = iNumSampsAudioPerChan/iBlockLen;

    std::vector<double> vTestAudioIn(iNumSampsAudio, 0.0);
    std::vector<double> vTestIR(iNumSampsIR, 1.0);

    //test input is a square wave with iFreq Hz
    for (unsigned int i=0; i<iNumSampsAudio; i+=iNumChansAudio)
    {
        fTmp = 0.5*(((sin(2.0*M_PI*iFreqAudio*((double)i/iNumChansAudio/iFs)))>0) ? (1):(-1));
        for (unsigned int k=0; k<iNumChansAudio; k++)
            vTestAudioIn.at(i+k) = fTmp;
    }

    //test impulse response is a damped oscillator with same frequency as test signal
    fAlpha = 1.0 - 1.0/(fTimeDurIR*0.1*iFs);
    for (unsigned int i=iNumChansIR; i<iNumSampsIR; i++)
        vTestIR.at(i) = fAlpha*vTestIR[i-iNumChansIR];

    for (unsigned int i=0; i<iNumSampsIR/2; i++) //test input is a delayed delta impulse
        vTestIR.at(i) = 0.0;

    for (unsigned int i=iNumSampsIR/2; i<iNumSampsIR; i+=iNumChansIR) //test input is a delayed delta impulse
    {
        fTmp = -sin(2.0*M_PI*iFreqIR*((double)(i-iNumSampsIR/2)/iNumChansIR/iFs));
        for (unsigned int k=0; k<iNumChansIR; k++)
            vTestIR.at(i+k) *= fTmp;
    }

    for (unsigned int i=0; i<iNumSampsIR; i++)
        fTmpSum += (vTestIR[i]>0) ? (vTestIR[i]) : (-vTestIR[i]);

    for (unsigned int i=0; i<iNumSampsIR; i++)
        vTestIR[i]/=fTmpSum;

    //---------------------------------WOLAP call----------------------------------------

    WOLAP wolapInst(vTestIR.data(), iLenIR, iNumChansIR, iBlockLen, iNumChansAudio);

    //---------------------------------GUI stuff-----------------------------------------

    QBrush plotBrush(QColor(150,150,150,150));
    uint32_t plotChannel = 0;

    inPlot = new QCustomPlot(this);
    inPlot->move(0,0);
    inPlot->resize(QSize(400,200));
    inPlot->addGraph();
    inPlot->graph(0)->setPen(QColor(230,50,50));
    inPlot->graph(0)->setBrush(plotBrush);
    inPlot->xAxis->setRange(1,iNumBlocks*iBlockLen);
    fMin = *std::min_element(vTestAudioIn.begin(),vTestAudioIn.end());
    fMax = *std::max_element(vTestAudioIn.begin(),vTestAudioIn.end());
    inPlot->yAxis->setRange(1.2*fMin,1.2*fMax);
    inPlot->setBackground(this->palette().background().color());
    xInPlot.resize(iNumSampsAudioPerChan);
    yInPlot.resize(iNumSampsAudioPerChan);
    for (unsigned int cnt=0; cnt<iNumBlocks*iBlockLen; cnt++)
    {
        xInPlot[cnt] = cnt;
        yInPlot[cnt] = vTestAudioIn[cnt*iNumChansAudio+plotChannel];
    }
    inPlot->graph(0)->setData(xInPlot, yInPlot, true);

    for (unsigned int i=0; i<iNumBlocks; i++)
        wolapInst.process(&vTestAudioIn[i*iNumChansAudio*iBlockLen]);

    for (unsigned int i=iNumBlocks*iBlockLen*iNumChansAudio; i<iNumSampsAudio; i++)
        vTestAudioIn[i] = 0.0;

    irPlot = new QCustomPlot(this);
    irPlot->move(400,0);
    irPlot->resize(QSize(400,200));
    irPlot->addGraph();
    irPlot->graph(0)->setPen(QColor(230,50,50));
    irPlot->graph(0)->setBrush(plotBrush);
    irPlot->xAxis->setRange(1,iLenIR);
    fMin = *std::min_element(vTestIR.begin(),vTestIR.end());
    fMax = *std::max_element(vTestIR.begin(),vTestIR.end());
    irPlot->yAxis->setRange(1.2*fMin, 1.2*fMax);
    irPlot->setBackground(this->palette().background().color());
    xIrPlot.resize(iLenIR);
    yIrPlot.resize(iLenIR);
    for (unsigned int cnt=0; cnt<iLenIR; cnt++)
    {
        xIrPlot[cnt] = cnt;
        yIrPlot[cnt] = vTestIR[cnt*iNumChansIR+plotChannel];
    }
    irPlot->graph(0)->setData(xIrPlot, yIrPlot, true);

    outPlot = new QCustomPlot(this);
    outPlot->move(800,0);
    outPlot->resize(QSize(400,200));
    outPlot->addGraph();
    outPlot->graph(0)->setPen(QColor(230,50,50));
    outPlot->graph(0)->setBrush(plotBrush);
    outPlot->xAxis->setRange(1,iNumBlocks*iBlockLen);
    fMin = *std::min_element(vTestAudioIn.begin(),vTestAudioIn.end());
    fMax = *std::max_element(vTestAudioIn.begin(),vTestAudioIn.end());
    outPlot->yAxis->setRange(1.2*fMin,1.2*fMax);
    outPlot->setBackground(this->palette().background().color());
    xOutPlot.resize(iNumSampsAudioPerChan);
    yOutPlot.resize(iNumSampsAudioPerChan);
    for (unsigned int cnt=0; cnt<iNumBlocks*iBlockLen; cnt++)
    {
        xOutPlot[cnt] = cnt;
        yOutPlot[cnt] = vTestAudioIn[cnt*iNumChansAudio+plotChannel];
    }
    outPlot->graph(0)->setData(xOutPlot, yOutPlot, true);
}

MainWindow::~MainWindow()
{
    delete inPlot;
    delete irPlot;
    delete outPlot;
}
