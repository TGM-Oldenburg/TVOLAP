/*----------------------------------------------------------------------------*\
| 32-bit Fixed-Point Fast Convolution using Partitioned Overlap Save (or Add)  |
| Public domain. License: see end of file.                                     |
|                                                                              |
| Authors:  Sascha Bilert, Jan Tinneberg, Hagen Jaeger,                        |
|           Christian Busse, Robert Liebchen, Uwe Simmer                       |
|                                                                              |
| History:  2016-10-18  initial version                                        |
|           2017-06-12  added selection of impulse response and cross-fade     |
|           2017-11-26  changed arrays to std::vector< >                       |
\*----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include "part_conv32.h"

#define M_PI 3.14159265358979323846

PartConv32::PartConv32(int blockLen, int impResLen, int maxNumImpRes,
                       int numChannels, int processChannel)
{
    this->blockLen = blockLen;
    this->impResLen = impResLen;
    this->maxNumImpRes = maxNumImpRes;
    this->numChannels = numChannels;
    this->processChannel = processChannel;

    numPart = (impResLen + blockLen - 1) / blockLen;    // number of partitions
    fftLen = 2 * blockLen;
    idxImpResOld = 0;
    crossfade_active = false;
    h.resize(fftLen);
    sumSpec.resize(fftLen / 2 + 1);
    input.resize(2 * blockLen);
    output.resize(2 * blockLen);
    outold.resize(2 * blockLen);
    outputMem.resize(blockLen);
    outoldMem.resize(blockLen);

    xSpec.resize(numPart);
    for (int i = 0; i < numPart; i++)
        xSpec[i].resize(fftLen / 2 + 1);

    hSpec.resize(maxNumImpRes);
    for (int i = 0; i < maxNumImpRes; i++)
    {
        hSpec[i].resize(numPart);
        for (int j = 0; j < numPart; j++)
            hSpec[i][j].resize(fftLen / 2 + 1);
    }

    partitionIndex.resize(numPart);
    for (int i = 0; i < numPart; i++)
        partitionIndex[i] = i;

    fadeIn.resize(blockLen);
    for (int i = 0; i < blockLen; i++)
        fadeIn[i]  = (int32_t)((0.5 - 0.5*cos(M_PI * i / blockLen)) * INT32_MAX);

    fadeOut.resize(blockLen);
    for (int i = 0; i < blockLen; i++)
        fadeOut[i] = (int32_t)((0.5 + 0.5*cos(M_PI * i / blockLen)) * INT32_MAX);

    // base 2 logarithm
    log2nfft = 0;
    for (int i = 1; i<fftLen; i *= 2)
        log2nfft++;
}

#if (BFP)

void PartConv32::normalize(std::vector<complex32>& spectrum32, std::vector<complex16>& spectrum16, int numSpect)
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
        spectrum16[n].re = (int16_t) (temp32 >> 16);
        temp32 = (spectrum32[n].im << block_exp) + 0x8000;
        spectrum16[n].im = (int16_t) (temp32 >> 16);
    }

    // save block exponent in the first imaginary part of the spectrum
    spectrum16[0].im = -block_exp;
}

#endif

void PartConv32::frequencyDomainConvolution(int idx)
{
    // set sum to zero
    for (int n = 0; n < fftLen / 2 + 1; n++)
        sumSpec[n].re = sumSpec[n].im = 0;

#if (BFP)
    // fast convolution in the frequency domain
    for (int p = 0; p < numPart; p++)
    {
        int pidx = partitionIndex[p];
        int block_exp = xSpec[pidx][0].im + hSpec[idx][p][0].im;

        if (block_exp <= 0)
        {
            block_exp = -block_exp;
            if (block_exp > 31)
                block_exp = 31;

            sumSpec[0].re += (((int32_t)xSpec[pidx][0].re * hSpec[idx][p][0].re) >> block_exp);

            for (int n = 1; n < fftLen / 2 + 1; n++)
            {
                sumSpec[n].re += ((((int32_t)xSpec[pidx][n].re * hSpec[idx][p][n].re)
                                 - ((int32_t)xSpec[pidx][n].im * hSpec[idx][p][n].im)) >> block_exp);

                sumSpec[n].im += ((((int32_t)xSpec[pidx][n].re * hSpec[idx][p][n].im)
                                 + ((int32_t)xSpec[pidx][n].im * hSpec[idx][p][n].re)) >> block_exp);
            }
        }
        else
        {
            if (block_exp > 31)
                block_exp = 31;

            sumSpec[0].re += (((int32_t)xSpec[pidx][0].re * hSpec[idx][p][0].re) << block_exp);

            for (int n = 1; n < fftLen / 2 + 1; n++)
            {
                sumSpec[n].re += ((((int32_t)xSpec[pidx][n].re * hSpec[idx][p][n].re)
                                 - ((int32_t)xSpec[pidx][n].im * hSpec[idx][p][n].im)) << block_exp);

                sumSpec[n].im += ((((int32_t)xSpec[pidx][n].re * hSpec[idx][p][n].im)
                                 + ((int32_t)xSpec[pidx][n].im * hSpec[idx][p][n].re)) << block_exp);
            }
        }
    }
#else
    // fast convolution in the frequency domain
    for (int p = 0; p < numPart; p++)
    {
        int pidx = partitionIndex[p];
        int block_exp = 31 - log2nfft;

        for (int n = 0; n < fftLen / 2 + 1; n++)
        {
            sumSpec[n].re += ((((int64_t)xSpec[pidx][n].re * hSpec[idx][p][n].re)
                             - ((int64_t)xSpec[pidx][n].im * hSpec[idx][p][n].im)) >> block_exp);

            sumSpec[n].im += ((((int64_t)xSpec[pidx][n].re * hSpec[idx][p][n].im)
                             + ((int64_t)xSpec[pidx][n].im * hSpec[idx][p][n].re)) >> block_exp);
        }
    }
#endif
}

void PartConv32::processBlock(int32_t *IOblock, int idxImpRes)
{
    if (idxImpRes < 0 || idxImpRes >= maxNumImpRes)
        return;

    // impulse respose changed (OLA only)
    if (!ols && idxImpRes != idxImpResOld && crossfade_active)
    {
        // save old overlap
        for (int i = 0; i < blockLen; i++)
            outoldMem[i] = outputMem[i];

        // frequency domain convolution (old data, new impulse response)
        frequencyDomainConvolution(idxImpRes);

        // inverse fast fourier transform
        irfft32(sumSpec.data(), output.data(), 0, fftLen);

        // new overlap
        for (int i = 0; i < blockLen; i++)
            outputMem[i] = output[i + blockLen];
    }

    // shift register for partition indices of xSpec
    int index = partitionIndex[numPart - 1];
    for (int i = numPart - 1; i > 0; i--)
        partitionIndex[i] = partitionIndex[i - 1];
    partitionIndex[0] = index;

    if (ols)    // overlap save
    {
        // move old samples
        for (int i = 0; i < blockLen; i++)
            input[i] = input[i + blockLen];

        // copy samples of current input channel
        for (int i = 0; i < blockLen; i++)
            input[i + blockLen] = IOblock[i*numChannels + processChannel];
    }
    else        // overlap add
    {
        // copy samples of current input channel
        for (int i = 0; i < blockLen; i++)
            input[i] = IOblock[i*numChannels + processChannel];

        // zero padding
        for (int i = 0; i < blockLen; i++)
            input[i + blockLen] = 0;
    }

#if (BFP)
    // fast fourier transform
    rfft32(input.data(), sumSpec.data(), 1, fftLen);

    // normalize spectrum and round to 16 bit
    normalize(sumSpec, xSpec[partitionIndex[0]], fftLen / 2 + 1);
#else
    // fast fourier transform
    rfft32(input.data(), xSpec[partitionIndex[0]].data(), 1, fftLen);
#endif

    // frequency domain convolution (new data, new impulse response)
    frequencyDomainConvolution(idxImpRes);

    // inverse fast fourier transform
    irfft32(sumSpec.data(), output.data(), 0, fftLen);

    if (!ols)   // overlap add
    {
        for (int i = 0; i < blockLen; i++)
            output[i] = output[i] + outputMem[i];

        for (int i = 0; i < blockLen; i++)
            outputMem[i] = output[i + blockLen];
    }

    // impulse respose changed (OLS or OLA)
    if (idxImpRes != idxImpResOld && crossfade_active)
    {
        // frequency domain convolution (new data, old impulse response)
        frequencyDomainConvolution(idxImpResOld);

        // inverse fast fourier transform
        irfft32(sumSpec.data(), outold.data(), 0, fftLen);

        if (ols)    // overlap save
        {
            // cross-fade
            for (int i = 0; i < blockLen; i++)
                output[i + blockLen] = ((int64_t)output[i + blockLen] * fadeIn[i]
                                      + (int64_t)outold[i + blockLen] * fadeOut[i]) >> 31;
        }
        else        // overlap add
        {
            for (int i = 0; i < blockLen; i++)
                outold[i] = outold[i] + outoldMem[i];

            // cross-fade
            for (int i = 0; i < blockLen; i++)
                output[i] = ((int64_t)output[i] * fadeIn[i]
                           + (int64_t)outold[i] * fadeOut[i]) >> 31;
        }
    }
    idxImpResOld = idxImpRes;

    if (ols)    // overlap save
    {
        for (int i = 0; i < blockLen; i++)
            IOblock[i*numChannels + processChannel] = output[i + blockLen];
    }
    else        // overlap add
    {
        for (int i = 0; i < blockLen; i++)
            IOblock[i*numChannels + processChannel] = output[i];
    }
}

void PartConv32::setImpRes(int32_t *impRes, int idxImpRes)
{
    // partition and transform impulse response
    int n, remainder = impResLen;

    if (idxImpRes < 0 || idxImpRes >= maxNumImpRes)
        return;

    n = blockLen;
    remainder = impResLen;

    for (int p = 0; p < numPart; p++)
    {
        // last partition of impRes may be shorter than blockLen
        if (remainder < blockLen)
            n = remainder;
        remainder = remainder - n;

        // copy impulse response
        for (int i = 0; i < n; i++)
            h[i] = impRes[i + blockLen * p];

        // zero padding
        for (int i = n; i < fftLen; i++)
            h[i] = 0;

#if (BFP)
        // compute transfer function of partition
        rfft32(h.data(), sumSpec.data(), 1, fftLen);

        // normalize transfer function and round to 16 bit
        normalize(sumSpec, hSpec[idxImpRes][p], fftLen / 2 + 1);

        // block exponent ist stored in the first imaginary part
        hSpec[idxImpRes][p][0].im += (log2nfft + 1);
#else
        // compute transfer function of partition
        rfft32(h.data(), hSpec[idxImpRes][p].data(), 1, fftLen);
#endif
    }
}

PartConv32::~PartConv32()
{
}

//--------------------- License ------------------------------------------------

// Copyright (c) 2016-2017 Sascha Bilert, Jan Tinneberg, Hagen Jaeger, 
//                         Christian Busse, Robert Liebchen, Uwe Simmer

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, 
// including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, 
// and to permit persons to whom the Software is furnished to do so, 
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//------------------------------------------------------------------------------
