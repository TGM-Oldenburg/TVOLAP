# -*- coding: utf-8 -*-
#
# This script evaluates TVOLAP against other known Real time procedures
# regarding their differences against each other
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 22.01.2017 (HJ)
# Ver. 0.2 debugged and tested manually 25.01.2017 (HJ)
#--------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from pyOLA import OLA
from pyOLS import OLS
from pyWOLA import WOLA
from pyTVOLAP import TVOLAP
import time
def loadKemarHRIR():
    IR0 = np.fromfile('Kemar_TUBerlin_0deg.bin', dtype='float')
    IR0 = np.reshape(IR0, (numChans, -1))
    
    IR90 = np.fromfile('Kemar_TUBerlin_90deg.bin', dtype='float')
    IR90 = np.reshape(IR90, (numChans, -1))
    
    IR = np.append(IR0[np.newaxis,:,:], IR90[np.newaxis,:,:], axis=0)
    
    return IR

if __name__ == '__main__':
    
    ### init some constant variables ###
    
    fs = int(44100)
    tDur = 5.0
    numSampsPerChan = int(fs*tDur/2)*2
    numChans = int(2)
    numIR = int(2)
    
    ### read HRIR responses for 0 and 90deg from TU berlin Kemar ###
    
    IR = loadKemarHRIR()
    lenIR = np.size(IR,2)
    
    ### generate sine wave as input signal ###
    
    inSig = np.random.randn(numChans, numSampsPerChan)*25
    weight = np.sqrt(1.0/(np.arange(numSampsPerChan/2+1)+1))
    inSig = np.fft.irfft(np.fft.rfft(inSig, axis=1)*weight, axis=1)
        
    ### calculate outputs via standard convolution in time domain ###
    
    outSig0 = np.zeros([numChans, numSampsPerChan])
    for chanCnt in np.arange(numChans):
        outSig0[chanCnt,:] = np.convolve(inSig[chanCnt,:], IR[0,chanCnt,:])[:numSampsPerChan]
        
    outSig90 = np.zeros([numChans, numSampsPerChan])
    for chanCnt in np.arange(numChans):
        outSig90[chanCnt,:] = np.convolve(inSig[chanCnt,:], IR[1,chanCnt,:])[:numSampsPerChan]
        
    ### calculate outputs via standard Overlap Add ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = numSampsPerChan/blockLen
    numBlocksTillSwitch = 16
    
    OLAinst = OLA(IR)
    outSigOLA = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLA[:,tmpIdxLo:tmpIdxHi] = OLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            OLAinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for OLA')
            
    ### calculate outputs via standard Overlap Save ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = numSampsPerChan/blockLen
    numBlocksTillSwitch = 16
    
    OLSinst = OLS(IR)
    outSigOLS = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLS[:,tmpIdxLo:tmpIdxHi] = OLSinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            OLSinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for OLS')
    
    ### calculate outputs via standard weighted overlap add RH ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = numSampsPerChan/blockLen
    numBlocksTillSwitch = 16
    
    WOLAinst = WOLA(IR, False)
    outSigWOLA_RH = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks-1):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigWOLA_RH[:,tmpIdxLo:tmpIdxHi] = WOLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            WOLAinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for WOLA RH')
    sampsDelay = int(blockLen/4)
    outSigWOLA_RH = np.append(outSigWOLA_RH[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### calculate outputs via standard weighted overlap add SHSH ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = numSampsPerChan/blockLen
    numBlocksTillSwitch = 16
    
    WOLAinst = WOLA(IR, True)
    outSigWOLA_SHSH = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks-1):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigWOLA_SHSH[:,tmpIdxLo:tmpIdxHi] = WOLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            WOLAinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for WOLA SHSH')
    sampsDelay = int(blockLen/4)
    outSigWOLA_SHSH = np.append(outSigWOLA_SHSH[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### calculate output via TVOLAP procedure to compare it ###
    
    blockLen = 256
    numBlocks = numSampsPerChan/blockLen
    numBlocksTillSwitch = 128
    
    TVOLAPinst = TVOLAP(IR, blockLen)
    outSigTVOLAP = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigTVOLAP[:,tmpIdxLo:tmpIdxHi] = TVOLAPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            TVOLAPinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for TVOLAP')
    sampsDelay = blockLen
    outSigTVOLAP = np.append(outSigTVOLAP[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### weight and add outputs via standard convolution for comparison ###
    
    numSampsTillSwitch = blockLen*numBlocksTillSwitch
    
    fadeInPos = np.arange(numSampsTillSwitch, numSampsPerChan, numSampsTillSwitch*2)-184
    fadeOutDist = numSampsTillSwitch
    
    hannFade = np.cos(np.linspace(0, np.pi/2, blockLen+1))[:-1]**2
    hannFade = np.repeat(hannFade[np.newaxis,:], numChans, axis=0)
    
    weight0 = np.ones([numChans, np.size(outSig0,1)])
    for switchCnt in (np.arange(np.size(fadeInPos))):
        fadeInBeg = fadeInPos[switchCnt]
        fadeInEnd = fadeInBeg+blockLen
        fadeOutBeg = fadeInBeg+fadeOutDist
        fadeOutEnd = fadeOutBeg+blockLen
        weight0[:,fadeInBeg:fadeInEnd] = hannFade
        weight0[:,fadeInEnd:fadeOutBeg] = 0.0
        weight0[:, fadeOutBeg:fadeOutEnd]  = hannFade[:,::-1]
        
    weight90 = -weight0+1.0
    
    outSigConv = outSig0*weight0+outSig90*weight90
    
    ### plot results ###
    
    plotSampWidth = 1024
    widthHalf = int(plotSampWidth/2)
    xAx = np.arange(plotSampWidth).astype(float)/fs*1000
    xLo = -1.1
    xHi = -xLo
    plt.figure(figsize=(4.1, 4))
    plt.subplot(2,1,1)
    plt.plot(xAx, outSigConv.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi])
    plt.subplot(2,1,2)
    tmp = outSigTVOLAP.T[32768-widthHalf:32768+widthHalf,:]-outSigConv.T[32768-widthHalf:32768+widthHalf,:]
    plt.plot(xAx, np.sqrt(tmp**2))
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([0, 1.05]);
    plt.xlabel('Time t [ms]', fontsize=12, fontname='cmr10')
    plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
    plt.tight_layout()
    plt.ylabel(' '*37 + 'Amplitude A [NV]', fontsize=12, fontname='cmr10')
    
#--------------------Licence ---------------------------------------------
# Copyright (c) 2012-2018 Hagen Jaeger                           
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