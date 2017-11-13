# -*- coding: utf-8 -*-
#
# Real time transfer function processing with special time-variant 
# implementation of the partitioned convolution in frequency domain 
# (to prevent from perceptive noticeable audio artifacts).                                 
#                                                                             
# Author: (c) Hagen Jaeger                      April 2016 - November 2017   
# MIT Release: Nov. 2017, License see end of file                              
#--------------------------------------------------------------------------

import numpy as np

class TVOLAP():
    
    def __init__(self, impulseResponse, blockLength):
        
        if not isinstance(impulseResponse, np.ndarray):
            impulseResponse = np.array(impulseResponse)
            
        if not isinstance(blockLength, int):
            blockLength = int(blockLength)
            
        if impulseResponse.ndim != 3:
            raise ValueError('ImpulseResponse response is not a 3dim matrix.')

        if blockLength != int(2**np.round(np.log2(blockLength))):
            raise ValueError('blockLength has to be a power of two (2^n).')
            
        if blockLength < 32:
            raise ValueError('blockLength has to be greather than 32.')
            
        self.blockLen = blockLength
        self.processLen = int(self.blockLen*2)
        self.nfft = int(2*self.processLen)
        self.numIR = np.size(impulseResponse, 0)
        self.numChans = np.size(impulseResponse, 1)
        
        self.numParts = int(np.ceil(np.size(impulseResponse,2)/self.processLen))
        self.numFreqMems = int(2*self.numParts)
        numSampsIR = int(self.numParts*self.processLen)
        impulseResponse = np.concatenate((impulseResponse, np.zeros([np.size(impulseResponse,0), 
                       np.size(impulseResponse,1), numSampsIR-np.size(impulseResponse,2)])),axis = 2)

        self.convMem = np.zeros([2,self.numChans,self.processLen])
        self.outDataMem = np.zeros([self.numChans, self.blockLen])

        self.win = 0.5+0.5*np.cos(np.linspace(-np.pi,np.pi,self.processLen+1))[:-1]
        
        self.freqIR = np.zeros([self.numIR, self.numParts, self.numChans, self.processLen+1]).astype(np.complex)
        for cntIR in np.arange(self.numIR):
            for partCnt in np.arange(self.numParts):
                for chanCnt in np.arange(self.numChans):
                    tmpIdxLo = partCnt*self.processLen
                    tmpIdxHi = tmpIdxLo+self.processLen
                    self.freqIR[cntIR, partCnt, chanCnt, :] = np.fft.rfft(
                            impulseResponse[cntIR, chanCnt, tmpIdxLo:tmpIdxHi], self.nfft)
    
        self.freqMem = np.zeros([self.numFreqMems, self.numChans, self.processLen+1]).astype(np.complex)
        self.freqSumReset = np.zeros([self.numChans, self.processLen+1]).astype(np.complex)
        
        self.freqSaveCnt = 0
        self.modCnt = 0
        self.outDataMem = np.zeros([self.numChans, self.blockLen])
        self.dataPre = np.zeros([self.numChans, self.blockLen])
        self.IR_ID = 0
        
    def process(self, data):
        inDataProcess = np.concatenate((self.dataPre,data), axis=1)*self.win
        self.dataPre = np.copy(data)
        self.freqMem[self.freqSaveCnt,:,:] = np.fft.rfft(inDataProcess, self.nfft)

        freqSum = np.copy(self.freqSumReset)
        
        freqReadCnt = self.freqSaveCnt
        for partCnt in np.arange(self.numParts):
            freqSum += self.freqMem[freqReadCnt,:,:]*self.freqIR[self.IR_ID,partCnt,:,:]
            freqReadCnt -= 2
            if freqReadCnt<0:
                freqReadCnt += self.numFreqMems

        tmpVec = np.fft.irfft(freqSum, self.nfft)
        outDataProcess = tmpVec[:,:self.processLen]+self.convMem[self.modCnt,:,:]
        self.convMem[self.modCnt,:,:] = np.copy(tmpVec[:,self.processLen:])

        data = outDataProcess[:,:self.blockLen]+self.outDataMem
        self.outDataMem = np.copy(outDataProcess[:,self.blockLen:])

        self.modCnt += 1
        self.freqSaveCnt += 1
        
        if self.modCnt>=2:
            self.modCnt = 0 
        if self.freqSaveCnt>=self.numFreqMems:
            self.freqSaveCnt = 0
            
        return data
        
    def setImResp(self, IR_ID):
            if (IR_ID >= 0) and (IR_ID < self.numIR):
                IR_ID = np.round(IR_ID)
                self.IR_ID = IR_ID
            else:
               raise ValueError('requested impulse response ID is out of range.')
        
#--------------------Licence ---------------------------------------------
# Copyright (c) 2012-2017 Hagen Jaeger                           
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the 
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE. 