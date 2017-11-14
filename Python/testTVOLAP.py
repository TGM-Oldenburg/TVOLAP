# -*- coding: utf-8 -*-
#
# This script evaluates an partitioned overlap add method.
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 13.11.2017 (HJ)
# Ver. 1.0 debugged and tested manually 13.11.2017 (HJ)
#--------------------------------------------------------------------------

from pyTVOLAP import TVOLAP
import numpy as np
import soundfile as sf
import matplotlib.pyplot as plt

fs = 48000
timeDurAudio = 5
timeDurIR = 0.2
numChans = 4
numIR = 20
blockLen = 512

impResp = np.zeros([numChans,int(np.round(fs*timeDurIR))])
deltaImp = np.append(np.array([1]), np.zeros(int(np.round(fs*timeDurIR*2)/2)-1))
deltaImp = np.fft.rfft(deltaImp)

freqCuts = np.round(np.logspace(np.log10(0.01),np.log10(0.99),numChans+1)*np.size(deltaImp,0)).astype(int)
for cnt in np.arange(numChans):
    deltaImpCut = np.copy(deltaImp)
    deltaImpCut[:freqCuts[cnt]] = 0.0
    deltaImpCut[freqCuts[cnt+1]:] = 0.0
    impResp[cnt,:] = np.fft.irfft(deltaImpCut)

tmpSize = int(np.size(impResp,1)/2)    
impResp = np.concatenate((impResp[:,tmpSize:], impResp[:,:tmpSize]), axis=1)
impResp = np.repeat([impResp,], numIR, axis=0)
tmpIdx = np.arange(numChans)
for cnt in np.arange(numIR):
    impResp[cnt, :, :] = impResp[cnt, tmpIdx, :]
    tmpIdx = np.append(tmpIdx[1:], tmpIdx[0])

inSig = np.random.randn(numChans, int(np.round(fs*timeDurAudio/blockLen)*blockLen))*0.1
specLen = int(np.size(inSig,1)/2+1)
weight = np.sqrt(1.0/(np.arange(specLen)+1)*specLen).astype(np.complex)
inSig = np.fft.irfft(np.fft.rfft(inSig, axis=1)*weight, axis=1)
    
TVOLAPinst = TVOLAP(impResp, blockLen)
outSig = np.copy(inSig)
actIR = 0
for blockCnt in np.arange(int(np.round(fs*timeDurAudio)/blockLen)):
    tmpIdxLo = blockCnt*blockLen
    tmpIdxHi = tmpIdxLo+blockLen
    outSig[:,tmpIdxLo:tmpIdxHi] = TVOLAPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
    if blockCnt%50 == 49:
        actIR = (actIR+1)%numIR
        TVOLAPinst.setImpResp(actIR)

sf.write('tstOut.wav', outSig.T, fs, 'PCM_16')
    
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