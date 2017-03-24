/*----------------------------------------------------------------------------*\
Header of WOLAP.cpp

Author: (c) Hagen Jaeger, Uwe Simmer    April 2016 - March 2017
\*----------------------------------------------------------------------------*/

#ifndef WOLAP_H
#define WOLAP_H

#include "complex_float64.h"
#include "fft.h"
#include "stdint.h"
#include <vector>

class WOLAP
{

public:
    WOLAP(double *vInterleavedIR, uint32_t iLengthIR, uint32_t iNumChansIR,
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
