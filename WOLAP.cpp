/*-----------------------------------------------------------------------------*\
| Real time transfer function processing with special time-variant 				|
| implementation of the partitioned convolution in frequency domain (to prevent |
| from perceptive noticeable audio artifacts).                                  |
|                                                                               |
| Author: (c) Hagen Jaeger, Uwe Simmer               April 2016 - March 2017    |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#include <stdexcept>
#include "WOLAP.h"

#define M_PI 3.14159265358979323846

WOLAP::WOLAP(std::vector<double> &interleavedIR, uint32_t numIR, uint32_t lenIR, uint32_t numChansIR,
             uint32_t blockLen, uint32_t numChansAudio)
{
    std::vector< std::vector< std::vector<double> > > tmpIR;
    std::vector<double> tmpPartIR;
    uint32_t intLenIR, irCnt = 0, chanCnt=0, convCnt=0, partCnt=0, memCnt=0, sampleCnt=0, cntIR=0;

    if (interleavedIR.size() != numIR*lenIR*numChansIR)
    	throw std::runtime_error("Size of interleaved impulse response is wrong."
    			"Must match number of IRs * length of one IR * number of IR channels.");

    this->blockLen = blockLen;
    this->numChansAudio = numChansAudio;
    this->numChansIR = numChansIR;
    this->numIR = numIR;
    this->processLen = 2*blockLen;
    this->nfft = 2*processLen;
    intLenIR = (((lenIR-1)/processLen+1)*processLen);
    this->numParts = intLenIR/processLen;
    this->overlapFact = 2;
    this->numMems = numParts*overlapFact;
    this->freqSaveCnt = 0;
    this->convSaveCnt = 0;
    this->actIR = 0;

    winVec.resize(processLen);
    for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
        winVec.at(sampleCnt) = 0.5-0.5*cos(2*M_PI*((double)sampleCnt/processLen));

    inBlockWin.resize(nfft);
    for (sampleCnt=0; sampleCnt<nfft; sampleCnt++)
        inBlockWin.at(sampleCnt) = 0.0;

    outBlock.resize(nfft);
    for (sampleCnt=0; sampleCnt<nfft; sampleCnt++)
        outBlock.at(sampleCnt) = 0.0;

    inSpectrumSum.resize(processLen+1);

    tmpIR.resize(numIR);
	for (irCnt=0; irCnt<numIR; irCnt++)
	{
		tmpIR.at(irCnt).resize(numChansIR);
		for (chanCnt=0; chanCnt<numChansIR; chanCnt++)
		{
			tmpIR.at(irCnt).at(chanCnt).resize(intLenIR);
			for (sampleCnt=0, cntIR=(irCnt*numChansIR+chanCnt)*lenIR; sampleCnt<lenIR; sampleCnt++, cntIR++)
				tmpIR.at(irCnt).at(chanCnt).at(sampleCnt) = interleavedIR.at(cntIR);

			for (sampleCnt=lenIR; sampleCnt<intLenIR; sampleCnt++)
				tmpIR.at(irCnt).at(chanCnt).at(sampleCnt) = 0.0;
		}
	}

    tmpPartIR.resize(nfft);
    for (sampleCnt=0; sampleCnt<nfft; sampleCnt++)
        tmpPartIR.at(sampleCnt) = 0.0;

    inBlock.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
    {
        inBlock.at(chanCnt).resize(processLen);
        for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
            inBlock.at(chanCnt).at(sampleCnt) = 0.0;
    }

    outBlockMat.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
    {
        outBlockMat.at(chanCnt).resize(processLen);
        for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
            outBlockMat.at(chanCnt).at(sampleCnt) = 0.0;
    }

    outBlockMatMem.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
    {
        outBlockMatMem.at(chanCnt).resize(blockLen);
        for (sampleCnt=0; sampleCnt<blockLen; sampleCnt++)
            outBlockMatMem.at(chanCnt).at(sampleCnt) = 0.0;
    }

    convMem.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
    {
        convMem.at(chanCnt).resize(overlapFact);
        for(convCnt=0; convCnt<overlapFact; convCnt++)
        {
            convMem.at(chanCnt).at(convCnt).resize(processLen);
            for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
                convMem.at(chanCnt).at(convCnt).at(sampleCnt) = 0.0;
        }
    }

    inSpectrum.resize(numChansIR);
    for (chanCnt=0; chanCnt<numChansIR; chanCnt++)
    {
        inSpectrum.at(chanCnt).resize(numMems);
        for (memCnt=0; memCnt<numMems; memCnt++)
        {
            inSpectrum.at(chanCnt).at(memCnt).resize(processLen+1);
        }
    }

    filterSpectrum.resize(numIR);
	for (irCnt=0; irCnt<numIR; irCnt++)
	{
		filterSpectrum.at(irCnt).resize(numChansIR);
		for (chanCnt=0; chanCnt<numChansIR; chanCnt++)
		{
			filterSpectrum.at(irCnt).at(chanCnt).resize(numParts);
			for (partCnt=0; partCnt<numParts; partCnt++)
			{
				filterSpectrum.at(irCnt).at(chanCnt).at(partCnt).resize(processLen+1);

				for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
					tmpPartIR.at(sampleCnt) = tmpIR.at(irCnt).at(chanCnt).at(partCnt*processLen+sampleCnt);

				rfft_double(tmpPartIR.data(), filterSpectrum.at(irCnt).at(chanCnt).at(partCnt).data(), nfft);
			}
		}
	}
}

void WOLAP::process(double *inBlockInterleaved)
{
    uint32_t chanCntAudio, iChanPosAudio, chanCntIR, partCnt, sampleCnt;
    int32_t freqReadCnt;

    for (chanCntAudio=0, chanCntIR=0; chanCntAudio<numChansAudio && chanCntIR<numChansIR; chanCntAudio++, chanCntIR++)
    {
        for (sampleCnt=0, iChanPosAudio=chanCntAudio; sampleCnt<blockLen; sampleCnt++, iChanPosAudio+=numChansAudio)
        {
            inBlock[chanCntAudio][sampleCnt] = inBlock[chanCntAudio][sampleCnt+blockLen];
            inBlock[chanCntAudio][sampleCnt+blockLen] = inBlockInterleaved[iChanPosAudio];
        }

        for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
            inBlockWin[sampleCnt] = inBlock[chanCntAudio][sampleCnt]*winVec[sampleCnt];

        rfft_double(inBlockWin.data(), inSpectrum[chanCntAudio][freqSaveCnt].data(), nfft);

        for (sampleCnt=0; sampleCnt<processLen+1; sampleCnt++)
            inSpectrumSum[sampleCnt].re = inSpectrumSum[sampleCnt].im = 0.0;

        freqReadCnt = freqSaveCnt;
        for (partCnt=0; partCnt<numParts; partCnt++)
        {
            for (sampleCnt=0; sampleCnt<processLen+1; sampleCnt++)
            {
                inSpectrumSum[sampleCnt] = complex_add(inSpectrumSum[sampleCnt],
                        complex_mul(inSpectrum[chanCntAudio][freqReadCnt][sampleCnt], filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt]));
            }

            freqReadCnt-=overlapFact;
            if (freqReadCnt<0)
                freqReadCnt+=numMems;
        }

        irfft_double(inSpectrumSum.data(), outBlock.data(), nfft);

        for (sampleCnt=0; sampleCnt<processLen; sampleCnt++)
        {
            outBlockMat[chanCntAudio][sampleCnt] = outBlock[sampleCnt]+convMem[chanCntAudio][convSaveCnt][sampleCnt];
            convMem[chanCntAudio][convSaveCnt][sampleCnt] = outBlock[sampleCnt+processLen];
        }

        for (sampleCnt=0, iChanPosAudio=chanCntAudio; sampleCnt<blockLen; sampleCnt++, iChanPosAudio+=numChansAudio)
        {
            inBlockInterleaved[iChanPosAudio] = outBlockMat[chanCntAudio][sampleCnt]+outBlockMatMem[chanCntAudio][sampleCnt];
            outBlockMatMem[chanCntAudio][sampleCnt] = outBlockMat[chanCntAudio][sampleCnt+blockLen];
        }
    }

    freqSaveCnt++;
    convSaveCnt++;
    if(freqSaveCnt>=numMems)
            freqSaveCnt=0;
    if(convSaveCnt>=overlapFact)
            convSaveCnt=0;
}

/*------------------------------License---------------------------------------*\
| Copyright (c) 2012-2017 Hagen Jaeger, Uwe Simmer                             |
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
