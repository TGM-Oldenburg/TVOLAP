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
    
    ### generate two sines as input signal ###
    
    f0 = 750
    inSig = np.sin(2*np.pi*f0*np.linspace(0, tDur, numSampsPerChan))
    inSig = np.repeat(inSig[np.newaxis,:], numChans, axis=0)
    
    ### evaluation change to ones vs. switching delta impulses ###
    
    #inSig.fill(1.0)
    #IR.fill(0)
    #IR[::2,:,0] = 1.0
    #IR[1::2,:,0] = -1.0
        
    ### calculate outputs via standard convolution in time domain ###
    
    lenIR = np.size(IR,2);
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
    sampsDelay = blockLen/2
    outSigTVOLAP = np.append(outSigTVOLAP[:,sampsDelay:], np.zeros([numChans, sampsDelay]), axis=1)
    
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
    
    plotSampWidth = 600
    widthHalf = int(plotSampWidth/2)
    xAx = np.arange(plotSampWidth).astype(float)/fs*1000
    xLo = -0.9
    xHi = -xLo
    plt.figure(figsize=(4.1, 7))
    plt.subplot(4,1,1)
    plt.plot(xAx, outSigOLA.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.ylim([xLo, xHi])
    plt.subplot(4,1,2)
    plt.plot(xAx, outSigOLS.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.ylim([xLo, xHi])
    plt.subplot(4,1,3)
    plt.plot(xAx, outSigWOLA_RH.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.ylim([xLo, xHi])
    plt.subplot(4,1,4)
    plt.plot(xAx, outSigTVOLAP.T[32768-widthHalf:32768+widthHalf,:])
    plt.grid(True); plt.ylim([xLo, xHi])
    plt.xlabel('Time t [ms]', fontsize=12, fontname='cmr10')
    plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
    plt.tight_layout()
    plt.ylabel(' '*110 + 'Amplitude A [NV]', fontsize=12, fontname='cmr10')
    
#==============================================================================
#     xAx = np.arange(500).astype(float)/fs
#     plt.figure(figsize=(4.1, 3))
#     plt.plot(xAx, outSigOLA.T[32500:33000,:])
#     plt.title('Overlap add procedure switching artifacts', fontsize=12, fontname='cmr10')
#     plt.xlabel('Time t [s]', fontsize=12, fontname='cmr10')
#     plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
#     plt.grid(True)
#     plt.tight_layout()
#     
#     plt.figure(figsize=(4.1, 3))
#     plt.plot(xAx, outSigOLS.T[32500:33000,:])
#     plt.title('Overlap save procedure switching artifacts', fontsize=12, fontname='cmr10')
#     plt.xlabel('Time t [s]', fontsize=12, fontname='cmr10')
#     plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
#     plt.grid(True)
#     plt.tight_layout()
#     
#     plt.figure(figsize=(4.1, 3))
#     plt.plot(xAx, outSigTVOLAP.T[32650:33150,:])
#     plt.title('TVOLAP procedure switching artifacts', fontsize=12, fontname='cmr10')
#     plt.xlabel('Time t [s]', fontsize=12, fontname='cmr10')
#     plt.ylabel('Amplitude A [NV]', fontsize=12, fontname='cmr10')
#     plt.grid(True)
#     plt.tight_layout()
#==============================================================================
    
    ### save all results as WAV files for easy playback ###
    
    #sf.write('outSig90_hardSwitch.wav', (outSig0*weight90+outSig90*weight270).T, fs)
    #sf.write('outSig90_TimeCrossfade.wav', (outSig0*weight90+outSig90*weight270).T, fs)