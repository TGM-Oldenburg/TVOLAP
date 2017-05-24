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

    //---------------------------------WOLAP init----------------------------------------

    WOLAP wolapInst(interleavedIR, numIR, numSampsIRPerChan, numChansIR, blockLen, numChansAudio);

    //-------------------------WOLAP Response calculation--------------------------------

    for (unsigned int i=0, j=0, k=0; i<numBlocks; i++, j++)
    {
    	if (j>=numBlocksPerIR)
    	{
    		j = 0;
    		wolapInst.setIR(++k);
    	}

        wolapInst.process(&testSignal[i*numChansAudio*blockLen]);
    }

    //------------------------------WAV file writing-------------------------------------

    uint16_t tmp16;
    uint32_t tmp32;
    std::vector<int16_t> tmp16Sig(testSignal.size(), 0);
    std::ofstream outFile("resultWOLAP.wav", std::ios::out | std::ios::binary);

    if (outFile.is_open())
    {
        outFile.write("RIFF", 4);
        tmp32 = testSignal.size()*sizeof(int16_t)+36;
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
        tmp32 = testSignal.size()*sizeof(int16_t)-44;
        outFile.write((char*) &tmp32, 4);

        for (uint32_t i=0; i< testSignal.size(); i++)
            tmp16Sig[i] = testSignal[i]*0x7FFF;

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
