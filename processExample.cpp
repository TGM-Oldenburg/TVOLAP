/*-----------------------------------------------------------------------------*\
| Example code for the WOLAP processing routine. Generates a test signal        |
| (square wave) and convolutes it with an delayed, damped oscillator. The       |
| result is stored as .wav file, so you can have easily have a look on it using |
| any waveform plotting tool (MATLAB, Audacity, a.s.o.)                         |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#include <stdint.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "WOLAP.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

int main()
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
    uint32_t tmp32;
    uint16_t tmp16;

    std::vector<double> vWOLAPIn(iNumSampsAudio, 0.0);
    std::vector<double> vDirectConvIn(iNumSampsAudio, 0.0);
    std::vector<double> vTestIR(iNumSampsIR, 1.0);
    std::vector<double> vDirectConvOut(iNumSampsAudio+iNumSampsIR-1, 0.0);
    std::vector<int16_t> vTmp16;

    //test input is a square wave with iFreq Hz
    for (uint32_t i=0; i<iNumSampsAudio; i+=iNumChansAudio)
    {
        fTmp = 0.5*(((sin(2.0*M_PI*iFreqAudio*((double)i/iNumChansAudio/iFs)))>0) ? (1):(-1));
        for (uint32_t k=0; k<iNumChansAudio; k++)
            vWOLAPIn.at(i+k) = vDirectConvIn.at(i+k) = fTmp;
    }

    //test impulse response is a damped oscillator with same frequency as test signal
    fAlpha = 1.0 - 1.0/(fTimeDurIR*0.1*iFs);
    for (uint32_t i=iNumChansIR; i<iNumSampsIR; i++)
        vTestIR.at(i) = fAlpha*vTestIR[i-iNumChansIR];

    for (uint32_t i=0; i<iNumSampsIR/2; i++) //test input is a delayed delta impulse
        vTestIR.at(i) = 0.0;

    for (uint32_t i=iNumSampsIR/2; i<iNumSampsIR; i+=iNumChansIR) //test input is a delayed delta impulse
    {
        fTmp = -sin(2.0*M_PI*iFreqIR*((double)(i-iNumSampsIR/2)/iNumChansIR/iFs));
        for (uint32_t k=0; k<iNumChansIR; k++)
            vTestIR.at(i+k) *= fTmp;
    }

    for (uint32_t i=0; i<iNumSampsIR; i++)
        fTmpSum += (vTestIR[i]>0) ? (vTestIR[i]) : (-vTestIR[i]);

    for (uint32_t i=0; i<iNumSampsIR; i++)
        vTestIR[i]/=fTmpSum;

    //---------------------------------WOLAP init----------------------------------------

    WOLAP wolapInst(vTestIR, iLenIR, iNumChansIR, iBlockLen, iNumChansAudio);

    //-------------------------WOLAP Response calculation--------------------------------

    for (uint32_t i=0; i<iNumBlocks; i++)
        wolapInst.process(&vWOLAPIn[i*iNumChansAudio*iBlockLen]);

    //------------------------------WAV file writing-------------------------------------

    std::ofstream outFile("resultWOLAP.wav", std::ios::out | std::ios::binary);

    if (outFile.is_open())
    {
        outFile.write("RIFF", 4);
        tmp32 = vWOLAPIn.size()*sizeof(int16_t)+36;
        outFile.write((char*) &tmp32, 4);
        outFile.write("WAVE", 4);

        outFile.write("fmt ", 4);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 4);
        tmp16 = 1;
        outFile.write((char*) &tmp16, 2);
        tmp16 = iNumChansAudio;
        outFile.write((char*) &tmp16, 2);
        outFile.write((char*) &iFs, 4);
        tmp32 = iFs*sizeof(int16_t)*iNumChansAudio;
        outFile.write((char*) &tmp32, 4); //bytes per second
        tmp32 = sizeof(int16_t)*iNumChansAudio;
        outFile.write((char*) &tmp32, 2);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 2);
        outFile.write("data", 4);
        tmp32 = vWOLAPIn.size()*sizeof(int16_t)-44;
        outFile.write((char*) &tmp32, 4);

        vTmp16.resize(vWOLAPIn.size());
        for (uint32_t i=0; i< vWOLAPIn.size(); i++)
            vTmp16[i] = vWOLAPIn[i]*0x7FFF;

        outFile.write((char*) vTmp16.data(), vTmp16.size()*sizeof(int16_t)-iNumChansAudio*iBlockLen);

        outFile.close();

        return 0;
    }
    else
        return -1;
}

/*------------------------------License---------------------------------------*\
| Copyright (c) 2012-2017 Hagen Jaeger                                         |
|                                                                              |
| This program is free software: you can redistribute it and/or modify         |
| it under the terms of the GNU Lesser General Public License as published by  |
| the Free Software Foundation, either version 3 of the License, or            |
| (at your option) any later version.                                          |
|                                                                              |
| This program is distributed in the hope that it will be useful,              |
| but WITHOUT ANY WARRANTY; without even the implied warranty of               |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                |
| GNU Lesser General Public License for more details.                          |
|                                                                              |
| You should have received a copy of the GNU Lesser General Public License     |
| along with this program. If not, see <http://www.gnu.org/licenses/>.         |
\*----------------------------------------------------------------------------*/
