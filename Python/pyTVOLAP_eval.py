# -*- coding: utf-8 -*-
#
# This script evaluates TVOLAP against other known Real time procedures
# Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
# Version History:
# Ver. 0.1 initial create (empty) 22.01.2017 (HJ)
# Ver. 0.2 debugged and tested manually 23.01.2017 (HJ)
#--------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
#import soundfile as sf
from pyOLA import OLA
from pyOLS import OLS
#from pyWOLA import WOLA
from pyTVOLAP import TVOLAP

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
    blockLen = 256
    numBlocksTillSwitch = 128
    numBlocks = int(np.round(numSampsPerChan)/blockLen)
    
    ### read HRIR responses for 0 and 90deg from TU berlin Kemar ###
    
    IR = loadKemarHRIR()
    
    ### generate two sines as input signal ###
    
    f0 = 300
    f1 = 600
    inSig = np.sin(2*np.pi*f0*np.linspace(0, tDur, numSampsPerChan))
    inSig = np.repeat(inSig[np.newaxis,:], numChans, axis=0)
        
    ### calculate outputs via standard convolution in time domain ###
    
    lenIR = np.size(IR,2);
    outSig0 = np.zeros([numChans, numSampsPerChan])
    for chanCnt in np.arange(numChans):
        outSig0[chanCnt,:] = np.convolve(inSig[chanCnt,:], IR[0,chanCnt,:])[:numSampsPerChan]
        
    outSig90 = np.zeros([numChans, numSampsPerChan])
    for chanCnt in np.arange(numChans):
        outSig90[chanCnt,:] = np.convolve(inSig[chanCnt,:], IR[1,chanCnt,:])[:numSampsPerChan]
        
    ### calculate outputs via standard Overlap Add ###
    
    OLAinst = OLA(IR, blockLen)
    outSigOLA = np.copy(inSig)
    actIR = 0
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLA[:,tmpIdxLo:tmpIdxHi] = OLAinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%128 == 127:
            actIR = (actIR+1)%numIR
            OLAinst.setImpResp(actIR)
            
    ### calculate outputs via standard Overlap Add ###
    
    OLSinst = OLS(IR, blockLen)
    outSigOLS = np.copy(inSig)
    actIR = 0
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigOLS[:,tmpIdxLo:tmpIdxHi] = OLSinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%128 == 127:
            actIR = (actIR+1)%numIR
            OLSinst.setImpResp(actIR)
    
    ### calculate output via TVOLAP procedure to compare it ###
    
    TVOLAPinst = TVOLAP(IR, blockLen)
    outSigTVOLAP = np.copy(inSig)
    actIR = 0
    for blockCnt in np.arange(numBlocks):
        tmpIdxLo = blockCnt*blockLen
        tmpIdxHi = tmpIdxLo+blockLen
        outSigTVOLAP[:,tmpIdxLo:tmpIdxHi] = TVOLAPinst.process(inSig[:,tmpIdxLo:tmpIdxHi])
        if blockCnt%128 == 127:
            actIR = (actIR+1)%numIR
            TVOLAPinst.setImpResp(actIR)
    
    ### weight and add outputs via standard convolution for comparison ###
    
    numSampsTillSwitch = blockLen*numBlocksTillSwitch
    
    fadeInPos = np.arange(numSampsTillSwitch, numSampsPerChan, numSampsTillSwitch*2)
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
    
    plt.plot(outSigOLS.T)
    
    ### save all results as WAV files for easy playback ###
    
    #sf.write('outSig90_hardSwitch.wav', (outSig0*weight90+outSig90*weight270).T, fs)
    #sf.write('outSig90_TimeCrossfade.wav', (outSig0*weight90+outSig90*weight270).T, fs)