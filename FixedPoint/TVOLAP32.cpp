/*-----------------------------------------------------------------------------*\
| Real time head related transfer function processing with special time-variant |
| implementation of the partitioned convolution in frequency domain (to prevent |
| from perceptive noticeable audio artifacts).                                  |
|                                                                               |
| Author: (c) Hagen Jaeger, Uwe Simmer               April 2016 - December 2017 |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#include <stdexcept>
#include "TVOLAP32.h"

#define M_PI 3.14159265358979323846

TVOLAP32::TVOLAP32(std::vector<int32_t> &interleavedIR, uint32_t numIR, uint32_t lenIR, uint32_t numChansIR,
                   uint32_t blockLen, uint32_t numChansAudio)
{
    std::vector< std::vector< std::vector<int32_t> > > tmpIR;
    std::vector<int32_t> tmpPartIR;
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
    for (sampleCnt = 0; sampleCnt < processLen; sampleCnt++)
        winVec.at(sampleCnt) = (int32_t)((0.5 - 0.5*cos(2. * M_PI * sampleCnt / processLen)) * INT32_MAX);

    inBlockWin.resize(nfft);

    outBlock.resize(nfft);

    inSpectrum32.resize(processLen+1);

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
                tmpIR.at(irCnt).at(chanCnt).at(sampleCnt) = 0;
        }
    }

    tmpPartIR.resize(nfft);

    inBlock.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
        inBlock.at(chanCnt).resize(processLen);

    outBlockMat.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
        outBlockMat.at(chanCnt).resize(processLen);

    outBlockMatMem.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
        outBlockMatMem.at(chanCnt).resize(blockLen);

    convMem.resize(numChansAudio);
    for(chanCnt=0; chanCnt<numChansAudio; chanCnt++)
    {
        convMem.at(chanCnt).resize(overlapFact);
        for(convCnt=0; convCnt<overlapFact; convCnt++)
            convMem.at(chanCnt).at(convCnt).resize(processLen);
    }

    inSpectrum.resize(numChansIR);
    for (chanCnt=0; chanCnt<numChansIR; chanCnt++)
    {
        inSpectrum.at(chanCnt).resize(numMems);
        for (memCnt=0; memCnt<numMems; memCnt++)
            inSpectrum.at(chanCnt).at(memCnt).resize(processLen+1);
    }

    // base 2 logarithm
    log2nfft = 0;
    for (uint32_t i = 1; i<nfft; i *= 2)
        log2nfft++;

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

#if (BLOCK_FLOATING_POINT)
                // compute transfer function of partition
                rfft32(tmpPartIR.data(), inSpectrum32.data(), 1, nfft);

                // normalize transfer function and round to 16 bit
                normalize(inSpectrum32, filterSpectrum[irCnt][chanCnt][partCnt], nfft / 2 + 1);

                // block exponent ist stored in the first imaginary part
                filterSpectrum[irCnt][chanCnt][partCnt][0].im += (log2nfft + 1);
#else
                // compute transfer function of partition
                rfft32(tmpPartIR.data(), filterSpectrum[irCnt][chanCnt][partCnt].data(), 1, nfft);
#endif
            }
        }
    }
}

void TVOLAP32::process(int32_t *inBlockInterleaved)
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
            inBlockWin[sampleCnt] = ((int64_t)inBlock[chanCntAudio][sampleCnt] * winVec[sampleCnt]) >> 31;

#if (BLOCK_FLOATING_POINT)
        // fast fourier transform
        rfft32(inBlockWin.data(), inSpectrum32.data(), 1, nfft);

        // normalize spectrum and round to 16 bit
        normalize(inSpectrum32, inSpectrum[chanCntAudio][freqSaveCnt], nfft / 2 + 1);
#else
        // fast fourier transform
        rfft32(inBlockWin.data(), inSpectrum[chanCntAudio][freqSaveCnt].data(), 1, nfft);
#endif

         // set sum to zero
        for (sampleCnt=0; sampleCnt<processLen+1; sampleCnt++)
            inSpectrum32[sampleCnt].re = inSpectrum32[sampleCnt].im = 0;

#if (BLOCK_FLOATING_POINT)
        freqReadCnt = freqSaveCnt;
        for (partCnt=0; partCnt<numParts; partCnt++)
        {
            int block_exp = inSpectrum[chanCntAudio][freqReadCnt][0].im + filterSpectrum[actIR][chanCntAudio][partCnt][0].im;

            if (block_exp <= 0)
            {
                block_exp = -block_exp;
                if (block_exp > 31)
                    block_exp = 31;

                inSpectrum32[0].re += (((int32_t)inSpectrum[chanCntAudio][freqReadCnt][0].re * filterSpectrum[actIR][chanCntAudio][partCnt][0].re) >> block_exp);

                for (sampleCnt = 1; sampleCnt < processLen + 1; sampleCnt++)
                {
                    inSpectrum32[sampleCnt].re += ((((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re)
                                                  - ((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im))
                                                    >> block_exp);

                    inSpectrum32[sampleCnt].im += ((((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im)
                                                  + ((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re))
                                                    >> block_exp);
                }
            }
            else
            {
                if (block_exp > 31)
                    block_exp = 31;

                inSpectrum32[0].re += (((int32_t)inSpectrum[chanCntAudio][freqReadCnt][0].re * filterSpectrum[actIR][chanCntAudio][partCnt][0].re) << block_exp);

                for (sampleCnt = 1; sampleCnt < processLen + 1; sampleCnt++)
                {
                    inSpectrum32[sampleCnt].re += ((((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re)
                                                  - ((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im))
                                                    << block_exp);

                    inSpectrum32[sampleCnt].im += ((((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im)
                                                  + ((int32_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re))
                                                    << block_exp);
                }
            }

            freqReadCnt -= overlapFact;
            if (freqReadCnt < 0)
                freqReadCnt += numMems;
        }
#else
        freqReadCnt = freqSaveCnt;
        for (partCnt = 0; partCnt<numParts; partCnt++)
        {
            int block_exp = 31 - log2nfft;

            for (sampleCnt = 0; sampleCnt < processLen + 1; sampleCnt++)
            {
                inSpectrum32[sampleCnt].re += ((((int64_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re)
                                              - ((int64_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im))
                                              >> block_exp);

                inSpectrum32[sampleCnt].im += ((((int64_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].re * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].im)
                                              + ((int64_t)inSpectrum[chanCntAudio][freqReadCnt][sampleCnt].im * filterSpectrum[actIR][chanCntAudio][partCnt][sampleCnt].re))
                                              >> block_exp);
            }

            freqReadCnt -= overlapFact;
            if (freqReadCnt < 0)
                freqReadCnt += numMems;
        }
#endif

        // inverse fast fourier transform
        irfft32(inSpectrum32.data(), outBlock.data(), 0, nfft);

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
    if(freqSaveCnt >= numMems)
        freqSaveCnt = 0;
    if(convSaveCnt >= overlapFact)
        convSaveCnt = 0;
}

#if (BLOCK_FLOATING_POINT)

void TVOLAP32::normalize(std::vector<complex32>& spectrum32, std::vector<complex16>& spectrum16, int numSpect)
{
    int i, n;
    int block_exp;
    int32_t max_value;
    int32_t abs_value;
    int32_t temp32;

    // maximum value of the spectrum
    max_value = 0;
    for (n = 0; n < numSpect; n++)
    {
        abs_value = labs(spectrum32[n].re);

        if (max_value < abs_value)
            max_value = abs_value;

        abs_value = labs(spectrum32[n].im);

        if (max_value < abs_value)
            max_value = abs_value;
    }

    // shift left until max_value >= 0.25
    for (i = 0; i < 32; i++)
    {
        if (max_value >= 0x20000000)
            break;

        max_value = max_value << 1;
    }

    // set block exponent
    block_exp = i;

    // normalize spectrum and round to 16 bit
    for (n = 0; n < numSpect; n++)
    {
        temp32 = (spectrum32[n].re << block_exp) + 0x8000;
        spectrum16[n].re = (int16_t)(temp32 >> 16);
        temp32 = (spectrum32[n].im << block_exp) + 0x8000;
        spectrum16[n].im = (int16_t)(temp32 >> 16);
    }

    // save block exponent in the first imaginary part of the spectrum
    spectrum16[0].im = -block_exp;
}

#endif

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
