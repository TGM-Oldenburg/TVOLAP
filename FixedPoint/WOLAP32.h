/*-----------------------------------------------------------------------------*\
| Header of WOLAP32.cpp, for explanation and example see cpp-file.              |
|                                                                               |
| Author: (c) Hagen Jaeger, Uwe Simmer               April 2016 - December 2017 |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#ifndef WOLAP32_H
#define WOLAP32_H

#include <stdint.h>
#include <vector>
#include "fft32.h"
#include "complex16.h"
#include "complex32.h"

#define BLOCK_FLOATING_POINT 0

class WOLAP32
{

public:
    WOLAP32(std::vector<int32_t> &interleavedIR, uint32_t numIR, uint32_t lenIR, uint32_t numChansIR,
            uint32_t blockLen, uint32_t numChansAudio);

    void process(int32_t *inBlockInterleaved);

    inline int setIR(uint32_t actIR)
    {
        if (actIR >= numIR)
            return -1;
        else
            this->actIR = actIR;
        return 0;
    }

private:
    uint32_t blockLen, processLen, nfft, numIR, actIR, numChansAudio, numChansIR, numParts, numMems, overlapFact, freqSaveCnt, convSaveCnt, log2nfft;
    std::vector<int32_t> winVec, inBlockWin, outBlock;
    std::vector<complex32> inSpectrum32;
    std::vector< std::vector<int32_t> > inBlock, outBlockMat, outBlockMatMem;
    std::vector< std::vector< std::vector<int32_t> > > convMem;
#if (BLOCK_FLOATING_POINT)
    std::vector< std::vector< std::vector<complex16> > > inSpectrum;
    std::vector< std::vector< std::vector< std::vector<complex16> > > > filterSpectrum;
    void normalize(std::vector<complex32>& spectrum32, std::vector<complex16>& spectrum16, int numSpect);
#else
    std::vector< std::vector< std::vector<complex32> > > inSpectrum;
    std::vector< std::vector< std::vector< std::vector<complex32> > > > filterSpectrum;
#endif
};

#endif // WOLAP32_H

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
