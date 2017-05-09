/*-----------------------------------------------------------------------------*\
| Header of WOLAP.cpp, for explanation and example see cpp-file.                |
|                                                                               |
| Author: (c) Hagen Jaeger, Uwe Simmer               April 2016 - March 2017    |
| LGPL Release: May 2017, License see end of file                               |
\*-----------------------------------------------------------------------------*/

#ifndef WOLAP_H
#define WOLAP_H

#include "fft.h"
#include "complex_functions.h"
#include <math.h>
#include <stdint.h>
#include <vector>

class WOLAP
{

public:
    WOLAP(std::vector<double> &vInterleavedIR, uint32_t iLengthIR, uint32_t iNumChansIR,
          uint32_t iBlockLen, uint32_t iNumChansAudio);

    void process(double *vInBlockInterleaved);

private:
    uint32_t iBlockLen, iProcessLen, iNfft, iNumChansAudio, iNumChansIR, iNumParts, iNumMems, iOverlapFact, iFreqSaveCnt, iConvSaveCnt;
    std::vector<double> vWin, vInBlockWin, vOutBlock;
    std::vector<complex_float64> vInSpectrumSum;
    std::vector< std::vector<double> > mInBlock, mOutBlock, mOutBlockMem;
    std::vector< std::vector< std::vector<double> > > mConvMem;
    std::vector< std::vector< std::vector<complex_float64> > > mInSpectrum, mFilterSpectrum;
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
