/*-----------------------------------------------------------------------------*\
| Example code for the TVOLAP processing routine. Generates a test signal        |
| (square wave) and convolutes it with an delayed, damped oscillator. The       |
| result is stored as .wav file, so you can have easily have a look on it using |
| any waveform plotting tool (MATLAB, Audacity, a.s.o.)                         |
|                                                                               |
| Author: (c) Hagen Jaeger                           April 2016 - March 2017    |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "TVOLAP.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

int main()
{
    //----------------------------TestSignalGeneration----------------------------------
    uint32_t sampFreq = 48000;
    double timeDurAudio = 5.0;
    double timeDurIR = 0.2;
    uint32_t numChans = 4;
    uint32_t numIR = 20;
    uint32_t blockLen = 512;

    uint32_t numSampsAudioPerChan = uint32_t(pow(2, round(log2(sampFreq*timeDurAudio))));
    uint32_t numSampsIRPerChan = uint32_t(pow(2, round(log2(sampFreq*timeDurIR))));
    uint32_t numBlocks = numSampsAudioPerChan/blockLen;

    std::vector<double> testSignal(numSampsAudioPerChan, 0.0);
    std::vector<double> testSignalInterleaved(numChans*numSampsAudioPerChan, 0.0);
    std::vector< std::vector< std::vector<double> > > testIR;
    std::vector<double> interleavedIR(numIR*numSampsIRPerChan*numChans, 0.0);
    std::vector<complex_float64> tmpFreq;
    double mem;
    srand(time(NULL));

    testIR.resize(numIR);
    for (uint32_t i = 0; i < numIR; i++)
    {
    	testIR.at(i).resize(numChans);
    	for (uint32_t j = 0; j < numChans; j++)
    		testIR.at(i).at(j).resize(numSampsIRPerChan, 0.0);
    }

    //test input is white noise
    tmpFreq.resize(numSampsAudioPerChan/2+1);
    for (uint32_t i=0; i<numSampsAudioPerChan; i++)
    {
    	testSignal.at(i) = double(rand()%0x10000)/double(0x8000)-1.0;

    	for (uint32_t j=0; j<10; j++)
    		testSignal.at(i) += double(rand()%0x10000)/double(0x8000)-1.0;

    	testSignal.at(i) *= 0.1;
    }

	rfft_double(testSignal.data(), tmpFreq.data(), numSampsAudioPerChan);

    for (uint32_t i=0; i<numSampsAudioPerChan/2+1; i++)
    	tmpFreq.at(i) = complex_mulr(tmpFreq.at(i), 0.33*sqrt(numSampsAudioPerChan/double(i+1)));

	irfft_double(tmpFreq.data(), testSignal.data(), numSampsAudioPerChan);

    for (uint32_t i=0; i<numSampsAudioPerChan-1; i++)
    {
    	for (uint32_t j=0; j<numChans; j++)
    		testSignalInterleaved.at(i*numChans+j) = testSignal.at(i);
    }

    //test impulse response are logarithmically spaced bandpass filters (read from file)
    tmpFreq.resize(numSampsIRPerChan/2+1);
    for (uint32_t i=0; i<numIR; i++)
    {
		for (uint32_t j=0; j<numChans; j++)
		{
			testIR.at(i).at(j).at(0) = 1.0;
			rfft_double(testIR.at(i).at(j).data(), tmpFreq.data(), numSampsIRPerChan);

			uint32_t k = (i+j)%numChans;
			double loBnd = pow(10, log10(0.99/0.01)/numChans*k+log10(numSampsIRPerChan/2*0.01));
			double upBnd = pow(10, log10(0.99/0.01)/numChans*(k+1)+log10(numSampsIRPerChan/2*0.01));

			for (uint32_t l=0; l<uint32_t(loBnd); l++)
				tmpFreq.at(l) = {0.0,0.0};


			for (uint32_t l=uint32_t(upBnd); l<numSampsIRPerChan/2+1; l++)
				tmpFreq.at(l) = {0.0,0.0};

			irfft_double(tmpFreq.data(), testIR.at(i).at(j).data(), numSampsIRPerChan);

			for (uint32_t l=0; l<numSampsIRPerChan/2; l++)
			{
				mem = testIR.at(i).at(j).at(l);
				testIR.at(i).at(j).at(l) = testIR.at(i).at(j).at(l+numSampsIRPerChan/2);
				testIR.at(i).at(j).at(l+numSampsIRPerChan/2) = mem;
			}
		}
    }

    for (uint32_t i=0; i<numIR; i++)
    {
    	for (uint32_t j=0; j<numChans; j++)
    	{
			for (uint32_t k=0; k<numSampsIRPerChan; k++)
				interleavedIR.at(k+(i*numChans+j)*numSampsIRPerChan) = testIR.at(i).at(j).at(k);
    	}
    }

    //---------------------------------TVOLAP init----------------------------------------

    TVOLAP TVOLAPInst(interleavedIR, numIR, numSampsIRPerChan, numChans, blockLen, numChans);

    //-------------------------TVOLAP Response calculation--------------------------------

    for (uint32_t i=0, j=0, k=0; i<numBlocks; i++, j++)
    {
    	if (j>=50)
    	{
    		j = 0;
    		TVOLAPInst.setIR(++k);
    	}

        TVOLAPInst.process(&testSignalInterleaved[i*numChans*blockLen]);
    }

    //------------------------------WAV file writing-------------------------------------

    uint16_t tmp16;
    uint32_t tmp32;
    std::vector<int16_t> tmp16Sig(testSignalInterleaved.size(), 0);
    std::ofstream outFile("tstOut.wav", std::ios::out | std::ios::binary);

    if (outFile.is_open())
    {
        outFile.write("RIFF", 4);
        tmp32 = testSignalInterleaved.size()*sizeof(int16_t)+36;
        outFile.write((char*) &tmp32, 4);
        outFile.write("WAVE", 4);

        outFile.write("fmt ", 4);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 4);
        tmp16 = 1;
        outFile.write((char*) &tmp16, 2);
        tmp16 = numChans;
        outFile.write((char*) &tmp16, 2);
        outFile.write((char*) &sampFreq, 4);
        tmp32 = sampFreq*sizeof(int16_t)*numChans;
        outFile.write((char*) &tmp32, 4); //bytes per second
        tmp32 = sizeof(int16_t)*numChans;
        outFile.write((char*) &tmp32, 2);
        tmp32 = 16;
        outFile.write((char*) &tmp32, 2);
        outFile.write("data", 4);
        tmp32 = testSignalInterleaved.size()*sizeof(int16_t)-44;
        outFile.write((char*) &tmp32, 4);

        for (uint32_t i=0; i< testSignalInterleaved.size(); i++)
        	tmp16Sig[i] = testSignalInterleaved[i]*0x7FFF;

        outFile.write((char*) tmp16Sig.data(), tmp16Sig.size()*sizeof(int16_t)-numChans*blockLen);

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
