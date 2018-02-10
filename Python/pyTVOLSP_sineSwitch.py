# -*- coding: utf-8 -*-
#
# This script evaluates TVOLSP against other known Real time procedures
# regarding their switching behaviour when convolving TU Berlin KEMAR HRTFs
# with a sinusodial signal.
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 22.01.2017 (HJ)
# Ver. 0.2 debugged and tested manually 23.01.2017 (HJ)
#--------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from pyOLA import OLA
from pyOLS import OLS
from pyWOLA import WOLA
from pyTVOLSP import TVOLSP
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
    
    f0 = 750
    inSig = np.sin(2*np.pi*f0*np.linspace(0, tDur, numSampsPerChan))
    inSig = np.repeat(inSig[np.newaxis,:], numChans, axis=0)
        
    ### calculate outputs via standard Overlap Add ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = int(numSampsPerChan/blockLen-1)
    numBlocksTillSwitch = 16
    
    OLAinst = OLA(IR)
    outSigOLA = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = int(blockCnt*blockLen)
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLA[:,tmpIdxLo:tmpIdxHi] = OLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            OLAinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for OLA')
            
    ### calculate outputs via standard Overlap Save ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = int(numSampsPerChan/blockLen-1)
    numBlocksTillSwitch = 16
    
    OLSinst = OLS(IR)
    outSigOLS = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = int(blockCnt*blockLen)
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLS[:,tmpIdxLo:tmpIdxHi] = OLSinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            OLSinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for OLS')
    
    ### calculate outputs via standard weighted overlap add RH ###
    
    blockLen = 2**(int(np.log2(lenIR-1))+1)
    numBlocks = int(numSampsPerChan/blockLen-1)
    numBlocksTillSwitch = 16
    
    WOLAinst = WOLA(IR, False)
    outSigWOLA_RH = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks-1):
        tmpIdxLo = int(blockCnt*blockLen)
        tmpIdxHi = tmpIdxLo+blockLen
        outSigWOLA_RH[:,tmpIdxLo:tmpIdxHi] = WOLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            WOLAinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for WOLA RH')
    sampsDelay = int(blockLen/4)
    outSigWOLA_RH = np.append(outSigWOLA_RH[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### calculate output via TVOLSP procedure to compare it ###
    
    blockLen = 256
    numBlocks = int(numSampsPerChan/blockLen-1)
    numBlocksTillSwitch = 128
    
    TVOLSPinst = TVOLSP(IR, blockLen)
    outSigTVOLSP = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = int(blockCnt*blockLen)
        tmpIdxHi = tmpIdxLo+blockLen
        outSigTVOLSP[:,tmpIdxLo:tmpIdxHi] = TVOLSPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            TVOLSPinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for TVOLSP')
    sampsDelay = int(blockLen/2)
    outSigTVOLSP = np.append(outSigTVOLSP[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### plot results ###
    
    plotSampWidth = 600
    widthHalf = int(plotSampWidth/2)
    xAx = np.arange(plotSampWidth).astype(float)/fs*1000
    xLo = -0.9
    xHi = -xLo
    plt.figure(figsize=(4.1, 7))
    plt.subplot(4,1,1)
    plt.plot(xAx, outSigOLA.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi])
    plt.subplot(4,1,2)
    plt.plot(xAx, outSigOLS.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi])
    plt.subplot(4,1,3)
    plt.plot(xAx, outSigWOLA_RH.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi])
    plt.subplot(4,1,4)
    plt.plot(xAx, outSigTVOLSP.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi])
    plt.xlabel('Time t [ms]', fontsize=12, fontname='cmr10')
    plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
    plt.tight_layout()
    plt.ylabel(' '*95 + 'Amplitude A [NV]', fontsize=12, fontname='cmr10')
    
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