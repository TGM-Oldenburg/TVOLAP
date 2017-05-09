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
    uint32_t sampFreq = 48000;
    double timeDurAudio = 2.0;
    double timeDurIR = 0.2;
    double tmpSum = 0.0001, alpha, tmpMin, tmpMax, tmpVar;
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
    uint32_t tmp32;
    uint16_t tmp16;

    std::vector<double> inputSignal(numSampsAudio, 0.0);
    std::vector<double> testIR(numSampsIR, 1.0);
    std::vector<int16_t> tmp16Sig;

    //test input is a square wave with iFreq Hz
    for (uint32_t i=0; i<numSampsAudio; i+=numChansAudio)
    {
        tmpVar = 0.5*(((sin(2.0*M_PI*squareWaveFreq*((double)i/numChansAudio/sampFreq)))>0) ? (1):(-1));
        for (uint32_t k=0; k<numChansAudio; k++)
            inputSignal.at(i+k) = tmpVar;
    }

    //test impulse response is a damped oscillator with same frequency as test signal
    alpha = 1.0 - 1.0/(timeDurIR*0.1*sampFreq);
    for (uint32_t i=numChansIR; i<numSampsIR; i++)
        testIR.at(i) = alpha*testIR[i-numChansIR];

    for (uint32_t i=0; i<numSampsIR/2; i++) //test input is a delayed delta impulse
        testIR.at(i) = 0.0;

    for (uint32_t i=numSampsIR/2; i<numSampsIR; i+=numChansIR) //test input is a delayed delta impulse
    {
        tmpVar = -sin(2.0*M_PI*oscillatorFreq*((double)(i-numSampsIR/2)/numChansIR/sampFreq));
        for (uint32_t k=0; k<numChansIR; k++)
            testIR.at(i+k) *= tmpVar;
    }

    for (uint32_t i=0; i<numSampsIR; i++)
        tmpSum += (testIR[i]>0) ? (testIR[i]) : (-testIR[i]);

    for (uint32_t i=0; i<numSampsIR; i++)
        testIR[i]/=tmpSum;

    //---------------------------------WOLAP init----------------------------------------

    WOLAP inputSignalst(testIR, lenIR, numChansIR, blockLen, numChansAudio);

    //-------------------------WOLAP Response calculation--------------------------------

    for (uint32_t i=0; i<numBlocks; i++)
        inputSignalst.process(&inputSignal[i*numChansAudio*blockLen]);

    //------------------------------WAV file writing-------------------------------------

    std::ofstream outFile("resultWOLAP.wav", std::ios::out | std::ios::binary);

    if (outFile.is_open())
    {
        outFile.write("RIFF", 4);
        tmp32 = inputSignal.size()*sizeof(int16_t)+36;
        outFile.write((char*) &tmp32, 4);
        outFile.write("WAVE", 4);

        outFile.write("fmt ", 4);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 4);
        tmp16 = 1;
        outFile.write((char*) &tmp16, 2);
        tmp16 = numChansAudio;
        outFile.write((char*) &tmp16, 2);
        outFile.write((char*) &sampFreq, 4);
        tmp32 = sampFreq*sizeof(int16_t)*numChansAudio;
        outFile.write((char*) &tmp32, 4); //bytes per second
        tmp32 = sizeof(int16_t)*numChansAudio;
        outFile.write((char*) &tmp32, 2);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 2);
        outFile.write("data", 4);
        tmp32 = inputSignal.size()*sizeof(int16_t)-44;
        outFile.write((char*) &tmp32, 4);

        tmp16Sig.resize(inputSignal.size());
        for (uint32_t i=0; i< inputSignal.size(); i++)
            tmp16Sig[i] = inputSignal[i]*0x7FFF;

        outFile.write((char*) tmp16Sig.data(), tmp16Sig.size()*sizeof(int16_t)-numChansAudio*blockLen);

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
