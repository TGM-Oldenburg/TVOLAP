#include "WOLAP.h"

WOLAP::WOLAP(double *vInterleavedIR, uint32_t iLengthIR, uint32_t iNumChansIR,
             uint32_t iBlockLen, uint32_t iNumChansAudio)
{
    std::vector< std::vector<double> > vTmpIR;
    std::vector<double> vTmpPartIR;
    uint32_t iIntLengthIR, iChanCnt=0, iConvCnt=0, iPartCnt=0, iMemCnt=0, iSampleCnt=0, iCntIR=0;

    this->iBlockLen = iBlockLen;
    this->iNumChansAudio = iNumChansAudio;
    this->iNumChansIR = iNumChansIR;
    this->iProcessLen = 2*iBlockLen;
    this->iNfft = 2*iProcessLen;
    iIntLengthIR = (((iLengthIR-1)/iProcessLen)*iProcessLen)+iProcessLen;
    this->iNumParts = iIntLengthIR/iProcessLen;
    this->iOverlapFact = 2;
    this->iNumMems = iNumParts*iOverlapFact;
    this->iFreqSaveCnt = 0;
    this->iConvSaveCnt = 0;

    vWin.resize(iProcessLen);
    for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
        vWin.at(iSampleCnt) = 0.5-0.5*cos(2*M_PI*((double)iSampleCnt/iProcessLen));

    vInBlockWin.resize(iNfft);
    for (iSampleCnt=0; iSampleCnt<iNfft; iSampleCnt++)
        vInBlockWin.at(iSampleCnt) = 0.0;

    vOutBlock.resize(iNfft);
    for (iSampleCnt=0; iSampleCnt<iNfft; iSampleCnt++)
        vOutBlock.at(iSampleCnt) = 0.0;

    vInSpectrumSum.resize(iProcessLen+1);

    vTmpIR.resize(iNumChansIR);
    for (iChanCnt=0; iChanCnt<iNumChansIR; iChanCnt++)
    {
        vTmpIR.at(iChanCnt).resize(iIntLengthIR);

        for (iSampleCnt=0, iCntIR=0; iSampleCnt<iLengthIR; iSampleCnt++, iCntIR+=iNumChansIR)
            vTmpIR.at(iChanCnt).at(iSampleCnt) = vInterleavedIR[iCntIR];

        for (iSampleCnt=iLengthIR; iSampleCnt<iIntLengthIR; iSampleCnt++)
            vTmpIR.at(iChanCnt).at(iSampleCnt) = 0.0;
    }

    vTmpPartIR.resize(iNfft);
    for (iSampleCnt=0; iSampleCnt<iNfft; iSampleCnt++)
        vTmpPartIR.at(iSampleCnt) = 0.0;

    mInBlock.resize(iNumChansAudio);
    for(iChanCnt=0; iChanCnt<iNumChansAudio; iChanCnt++)
    {
        mInBlock.at(iChanCnt).resize(iProcessLen);
        for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
            mInBlock.at(iChanCnt).at(iSampleCnt) = 0.0;
    }

    mOutBlock.resize(iNumChansAudio);
    for(iChanCnt=0; iChanCnt<iNumChansAudio; iChanCnt++)
    {
        mOutBlock.at(iChanCnt).resize(iProcessLen);
        for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
            mOutBlock.at(iChanCnt).at(iSampleCnt) = 0.0;
    }

    mOutBlockMem.resize(iNumChansAudio);
    for(iChanCnt=0; iChanCnt<iNumChansAudio; iChanCnt++)
    {
        mOutBlockMem.at(iChanCnt).resize(iBlockLen);
        for (iSampleCnt=0; iSampleCnt<iBlockLen; iSampleCnt++)
            mOutBlockMem.at(iChanCnt).at(iSampleCnt) = 0.0;
    }

    mConvMem.resize(iNumChansAudio);
    for(iChanCnt=0; iChanCnt<iNumChansAudio; iChanCnt++)
    {
        mConvMem.at(iChanCnt).resize(iOverlapFact);
        for(iConvCnt=0; iConvCnt<iOverlapFact; iConvCnt++)
        {
            mConvMem.at(iChanCnt).at(iConvCnt).resize(iProcessLen);
            for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
                mConvMem.at(iChanCnt).at(iConvCnt).at(iSampleCnt) = 0.0;
        }
    }

    mInSpectrum.resize(iNumChansIR);
    for (iChanCnt=0; iChanCnt<iNumChansIR; iChanCnt++)
    {
        mInSpectrum.at(iChanCnt).resize(iNumMems);
        for (iMemCnt=0; iMemCnt<iNumMems; iMemCnt++)
        {
            mInSpectrum.at(iChanCnt).at(iMemCnt).resize(iProcessLen+1);
        }
    }

    mFilterSpectrum.resize(iNumChansIR);
    for (iChanCnt=0; iChanCnt<iNumChansIR; iChanCnt++)
        {
        mFilterSpectrum.at(iChanCnt).resize(iNumParts);
        for (iPartCnt=0; iPartCnt<iNumParts; iPartCnt++)
        {
            mFilterSpectrum.at(iChanCnt).at(iPartCnt).resize(iProcessLen+1);

            for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
                vTmpPartIR.at(iSampleCnt) = vTmpIR.at(iChanCnt).at(iPartCnt*iProcessLen+iSampleCnt);

            fft_double(vTmpPartIR.data(), mFilterSpectrum.at(iChanCnt).at(iPartCnt).data(), iNfft);
        }
    }
}

void WOLAP::process(double *vInBlockInterleaved)
{
    uint32_t iChanCntAudio, iChanPosAudio, iChanCntIR, iPartCnt, iSampleCnt;
    int32_t iFreqReadCnt;

    for (iChanCntAudio=0, iChanCntIR=0; iChanCntAudio<iNumChansAudio && iChanCntIR<iNumChansIR; iChanCntAudio++, iChanCntIR++)
    {
        for (iSampleCnt=0, iChanPosAudio=iChanCntAudio; iSampleCnt<iBlockLen; iSampleCnt++, iChanPosAudio+=iNumChansAudio)
        {
            mInBlock[iChanCntAudio][iSampleCnt] = mInBlock[iChanCntAudio][iSampleCnt+iBlockLen];
            mInBlock[iChanCntAudio][iSampleCnt+iBlockLen] = vInBlockInterleaved[iChanPosAudio];
        }

        for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
            vInBlockWin[iSampleCnt] = mInBlock[iChanCntAudio][iSampleCnt]*vWin[iSampleCnt];

        fft_double(vInBlockWin.data(), mInSpectrum[iChanCntAudio][iFreqSaveCnt].data(), iNfft);

        for (iSampleCnt=0; iSampleCnt<iProcessLen+1; iSampleCnt++)
            vInSpectrumSum[iSampleCnt].re = vInSpectrumSum[iSampleCnt].im = 0.0;

        iFreqReadCnt = iFreqSaveCnt;
        for (iPartCnt=0; iPartCnt<iNumParts; iPartCnt++)
        {
            for (iSampleCnt=0; iSampleCnt<iProcessLen+1; iSampleCnt++)
            {
                vInSpectrumSum[iSampleCnt] = complex_add(vInSpectrumSum[iSampleCnt],
                        complex_mul(mInSpectrum[iChanCntAudio][iFreqReadCnt][iSampleCnt], mFilterSpectrum[iChanCntAudio][iPartCnt][iSampleCnt]));
            }

            iFreqReadCnt-=iOverlapFact;
            if (iFreqReadCnt<0)
                iFreqReadCnt+=iNumMems;
        }

        ifft_double(vInSpectrumSum.data(), vOutBlock.data(), iNfft);

        for (iSampleCnt=0; iSampleCnt<iProcessLen; iSampleCnt++)
        {
            mOutBlock[iChanCntAudio][iSampleCnt] = vOutBlock[iSampleCnt]+mConvMem[iChanCntAudio][iConvSaveCnt][iSampleCnt];
            mConvMem[iChanCntAudio][iConvSaveCnt][iSampleCnt] = vOutBlock[iSampleCnt+iProcessLen];
        }

        for (iSampleCnt=0, iChanPosAudio=iChanCntAudio; iSampleCnt<iBlockLen; iSampleCnt++, iChanPosAudio+=iNumChansAudio)
        {
            vInBlockInterleaved[iChanPosAudio] = mOutBlock[iChanCntAudio][iSampleCnt]+mOutBlockMem[iChanCntAudio][iSampleCnt];
            mOutBlockMem[iChanCntAudio][iSampleCnt] = mOutBlock[iChanCntAudio][iSampleCnt+iBlockLen];
        }
    }

    iFreqSaveCnt++;
    iConvSaveCnt++;
    if(iFreqSaveCnt>=iNumMems)
            iFreqSaveCnt=0;
    if(iConvSaveCnt>=iOverlapFact)
            iConvSaveCnt=0;
}
