# -*- coding: utf-8 -*-
#
# This script evaluates an partitioned overlap add method.
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 13.11.2017 (HJ)
# Ver. 1.0 debugged and tested manually 13.11.2017 (HJ)
#--------------------------------------------------------------------------

from pyTVOLAP import TVOLAP
import matplotlib.pyplot as plt
import numpy as np

fs = 48000
timeDurAudio = 2
freqAudio = 2
timeDurIR = 0.1
numChans = 4

blockLen = 1024
processLen = 2*blockLen
nfft = 2*processLen

impResp = np.zeros([numChans,int(np.round(fs*timeDurIR))])

deltaImp = np.concatenate((np.array([1]), np.zeros(int(np.round(fs*timeDurIR*2)/2)-1)))

deltaImp = np.fft.rfft(deltaImp)

freqCuts = np.linspace(0.01,0.99,numChans+1)
numFreqs = np.size(deltaImp,0)
for cnt in np.arange(numChans):
    deltaImpCut = np.copy(deltaImp)
    deltaImpCut[:int(numFreqs*freqCuts[cnt])] = 0.0
    deltaImpCut[int(numFreqs*freqCuts[cnt+1]):] = 0.0
    impResp[cnt,:] = np.fft.irfft(deltaImpCut)

tmpSize = int(np.size(impResp,1)/2)    
impResp = np.concatenate((impResp[:,tmpSize:], impResp[:,:tmpSize]),axis=1)

xAx = np.linspace(0, timeDurAudio, int(np.round(fs*timeDurAudio)))
inSig = np.repeat([np.sign(np.sin(2*np.pi*freqAudio*xAx)),],numChans, axis=0)

impResp = impResp[np.newaxis,:,:]

TVOLAPinst = TVOLAP(impResp, blockLen)

numBlocks = int(np.round(fs*timeDurAudio)/blockLen)

outSig = np.copy(inSig)

for blockCnt in np.arange(numBlocks):
    tmpIdxLo = blockCnt*blockLen
    tmpIdxHi = tmpIdxLo+blockLen
    outSig[:,tmpIdxLo:tmpIdxHi] = TVOLAPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
    
plt.plot(inSig[0,:])
plt.plot(outSig[0,:])

    
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