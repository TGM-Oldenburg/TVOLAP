# -*- coding: utf-8 -*-
#
# Real time transfer function processing with weighted overlap add method
# implementation of the partitioned convolution in frequency domain 
# (to prevent from perceptive noticeable audio artifacts).                                 
#                                                                             
# Author: (c) Hagen Jaeger                      April 2016 - January 2018  
# MIT Release: Nov. 2017, License see end of file                              
#--------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

class WOLA():
    
    def __init__(self, impulseResponse, synthWin=False):
        
        if not isinstance(impulseResponse, np.ndarray):
            impulseResponse = np.array(impulseResponse)
            
        if impulseResponse.ndim != 3:
            raise ValueError('ImpulseResponse response is not a 3dim matrix.')
            
        self.numIR = np.size(impulseResponse, 0)
        self.numChans = np.size(impulseResponse, 1)
        self.blockLen = np.size(impulseResponse, 2)
        self.hopSize = int(self.blockLen/2)
        
        if self.numChans > self.blockLen:
            raise ValueError('channels > samples per IR. Bad formatted matrix.')

        self.nfft = int(2*self.blockLen)
        self.IR_ID = 0

        self.convMem = np.zeros([2,self.numChans,self.blockLen])
        self.convMemCnt = 0
        
        self.inDataProcess = np.zeros([self.numChans, self.blockLen])
        self.outMem = np.zeros([self.numChans, self.hopSize])
        
        self.freqIR = np.fft.rfft(impulseResponse, n=self.nfft, axis=2)
        
        self.win = 0.5+0.5*np.cos(np.linspace(-np.pi,np.pi,self.blockLen+1))[:-1]
        self.synthWin = synthWin
        if synthWin:
            self.win = np.sqrt(self.win)
    
    def process(self, data):
        self.inDataProcess[:,:self.hopSize] = self.inDataProcess[:,self.hopSize:]
        self.inDataProcess[:,self.hopSize:] = data[:,:self.hopSize]
        freqIn = np.fft.rfft(self.inDataProcess*self.win, n=self.nfft, axis=1)
        tmpOut = np.fft.irfft(freqIn*self.freqIR[self.IR_ID,:,:], n=self.nfft, axis=1)
        outSig = (tmpOut[:,:self.blockLen]+self.convMem[0,:,:])
        self.convMem[0,:,:] = tmpOut[:,self.blockLen:]
        if self.synthWin:
            outSig *= self.win
        
        self.inDataProcess[:,:self.hopSize] = self.inDataProcess[:,self.hopSize:]
        self.inDataProcess[:,self.hopSize:] = data[:,self.hopSize:]
        freqIn = np.fft.rfft(self.inDataProcess*self.win, n=self.nfft, axis=1)
        tmpOut = np.fft.irfft(freqIn*self.freqIR[self.IR_ID,:,:], n=self.nfft, axis=1)
        outSigTmp = (tmpOut[:,:self.blockLen]+self.convMem[1,:,:])
        self.convMem[1,:,:] = tmpOut[:,self.blockLen:]
        if self.synthWin:
            outSigTmp *= self.win
        
        outSig[:,:self.hopSize] += self.outMem
        outSig[:,self.hopSize:] += outSigTmp[:,:self.hopSize]
        self.outMem = outSigTmp[:,self.hopSize:]

        return outSig
        
        
    def setImpResp(self, IR_ID):
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