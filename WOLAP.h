/*-----------------------------------------------------------------------------*\
| Header of WOLAP.cpp, for explanation and example see cpp-file.                |
|                                                                               |
| Author: (c) Hagen Jaeger, Uwe Simmer               April 2016 - March 2017    |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#ifndef WOLAP_H
#define WOLAP_H

#include <stdint.h>
#include <vector>
#include "fft.h"
#include "complex_functions.h"

class WOLAP
{

public:
    WOLAP(std::vector<double> &interleavedIR, uint32_t numIR, uint32_t lenIR, uint32_t numChansIR,
          uint32_t blockLen, uint32_t numChansAudio);

    void process(double *inBlockInterleaved);

    inline int setIR(uint32_t actIR)
    {
    	if (actIR >= numIR)
			return -1;
    	else
    		this->actIR = actIR;
    	return 0;
    }

private:
    uint32_t blockLen, processLen, nfft, numIR, actIR, numChansAudio, numChansIR, numParts, numMems, overlapFact, freqSaveCnt, convSaveCnt;
    std::vector<double> winVec, inBlockWin, outBlock;
    std::vector<complex_float64> inSpectrumSum;
    std::vector< std::vector<double> > inBlock, outBlockMat, outBlockMatMem;
    std::vector< std::vector< std::vector<double> > > convMem;
    std::vector< std::vector< std::vector<complex_float64> > > inSpectrum;
    std::vector< std::vector< std::vector< std::vector<complex_float64> > > > filterSpectrum;
};

#endif // WOLAP_H

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
