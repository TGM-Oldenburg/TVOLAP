#ifndef _PartConv32_H
#define _PartConv32_H

#include <stdint.h>
#include <vector>

#include "fft32.h"
#include "complex16.h"
#include "complex32.h"

#define BFP 0

class PartConv32
{
public:
    PartConv32(int blockLen, int impResLen, int maxNumImpRes, int numChannels, int processChannel);
    void processBlock(int32_t *IOblock, int idxImpRes);
    void setImpRes(int32_t *impRes, int idxImpRes);
    void setCrossfadeActive(bool active) { crossfade_active = active; };
    void setOLS(bool active) { ols = active; };
    ~PartConv32();

private:
    int blockLen;
    int impResLen;
    int maxNumImpRes;
    int numChannels;
    int processChannel;
    int numPart;
    int fftLen;
    int idxImpResOld;
    int log2nfft;
    bool crossfade_active;
    bool ols;
    std::vector<int32_t> h;
    std::vector<int32_t> input;
    std::vector<int32_t> output;
    std::vector<int32_t> outold;
    std::vector<int32_t> outputMem;
    std::vector<int32_t> outoldMem;
    std::vector<complex32> sumSpec;
#if (BFP)
    std::vector< std::vector<complex16> > xSpec;
    std::vector< std::vector< std::vector<complex16> > > hSpec;
#else
    std::vector< std::vector<complex32> > xSpec;
    std::vector< std::vector< std::vector<complex32> > > hSpec;
#endif
    std::vector<int32_t> partitionIndex;
    std::vector<int32_t> fadeIn;
    std::vector<int32_t> fadeOut;
#if (BFP)
    void normalize(std::vector<complex32>& spectrum32, std::vector<complex16>& spectrum16, int numSpect);
#endif
    void frequencyDomainConvolution(int idx);
};

#endif  // _PartConv32_H
