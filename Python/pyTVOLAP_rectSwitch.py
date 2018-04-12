# -*- coding: utf-8 -*-
#
# This script evaluates TVOLAP against other known Real time procedures
# regarding their switching behaviour when convolving delta impulses
# (switching polarity) with a DC signal (output is a square wave)
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 22.01.2017 (HJ)
# Ver. 0.2 debugged and tested manually 24.01.2017 (HJ)
#--------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from pyOLA import OLA
from pyWOLA import WOLA
from pyTVOLAP import TVOLAP
import time

if __name__ == '__main__':
    
    ### init some constant variables ###
    
    fs = int(44100)
    tDur = 5.0
    numSampsPerChan = int(fs*tDur/2)*2
    numChans = int(2)
    numIR = int(2)
    
    ### generate ones as input signal and two polarity switching impulses ###
    
    lenIR = 2048
    inSig = np.ones([numChans, numSampsPerChan])
    IR = np.zeros([2, numChans, lenIR])
    IR[::2,:,0] = 1.0
    IR[1::2,:,0] = -1.0
        
    ### calculate outputs via standard Overlap Add ###
    
    blockLen = int(2**(int(np.log2(lenIR-1))+1))
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
    
    ### calculate outputs via standard weighted overlap add RH ###
    
    blockLen = int(2**(int(np.log2(lenIR-1))+1))
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
     
    ### calculate output via TVOLAP procedure to compare it ###
    
    blockLen = 256
    numBlocks = int(numSampsPerChan/blockLen-1)
    numBlocksTillSwitch = 128
    
    TVOLAPinst = TVOLAP(IR, blockLen)
    outSigTVOLAP = np.copy(inSig)
    actIR = 0
    start_time = time.time()
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = int(blockCnt*blockLen)
        tmpIdxHi = tmpIdxLo+blockLen
        outSigTVOLAP[:,tmpIdxLo:tmpIdxHi] = TVOLAPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%numBlocksTillSwitch == numBlocksTillSwitch-1:
            actIR = (actIR+1)%numIR
            TVOLAPinst.setImpResp(actIR)
    print(str(time.time() - start_time)+ ' seconds for TVOLAP')
    sampsDelay = int(blockLen/2)
    outSigTVOLAP = np.append(outSigTVOLAP[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
    ### plot results ###
    
    plotSampWidth = 1024
    widthHalf = int(plotSampWidth/2)
    xAx = np.arange(plotSampWidth).astype(float)/fs*1000
    xTicks = np.arange(0, plotSampWidth*1000/fs, plotSampWidth*1000/fs/4.001)
    xLo = -1.1
    xHi = -xLo
    plt.figure(figsize=(4.1, 2.5))
    plt.plot(xAx, outSigOLA.T[32768-widthHalf:32768+widthHalf,:])
    plt.plot(xAx, outSigWOLA_RH.T[32768-widthHalf:32768+widthHalf,:])
    plt.plot(xAx, outSigTVOLAP.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.xlim([xAx[0], xAx[-1]]); plt.ylim([xLo, xHi]);
    plt.xlabel('Time t [ms]', fontsize=12, fontname='cmr10')
    plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
    plt.xticks(xTicks, [r'$0$', r'$\frac{L}{4 \cdot fs}$', r'$\frac{L}{2 \cdot fs}$', 
                        r'$\frac{3 \cdot L}{4 \cdot fs}$', r'$L$'])
    plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1)
    
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