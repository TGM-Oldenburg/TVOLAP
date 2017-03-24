% This script generates a .wav file which consists of 24 channels and owns
% interleaved noise bursts per channel over time (first 0.1s noise burst on
% channel 1, second 0.1s noise burst on channel 2 a.s.o.
% Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
% Version History:
% Ver. 0.1 initial create (empty) 08.11.2016 HJ
% Ver. 1.0 debugged and tested 08.11.2016 HJ

clear; close all; clc

set(0,'defaulttextinterpreter','latex');

iFs = 44100;
fTimeDurAudio = 15;
fTimeDurBurst = 0.5;

iNumChans = 24;
iNumSampsPerBurst = round(iFs*fTimeDurBurst);
iNumSampsPerChan = round(iFs*fTimeDurAudio);

mAudioData = zeros(iNumSampsPerChan, iNumChans);

vWin = tukeywin(iNumSampsPerBurst,0.5);

iChanCnt = 1;
iSampCnt = 1;
while (iSampCnt <= iNumSampsPerChan-iNumSampsPerBurst)
   
    vTmp = 0.05*randn(iNumSampsPerBurst,1);
    
    mAudioData(iSampCnt:iSampCnt+iNumSampsPerBurst-1,iChanCnt) = vTmp.*vWin;
    
    iChanCnt = iChanCnt+1;
    iSampCnt = iSampCnt+iNumSampsPerBurst;
    if iChanCnt > iNumChans, iChanCnt = 1; end
end

figure('units', 'normalized', 'position', [0.05, 0.05, 0.5, 0.85])
for iChanCnt = 1:iNumChans
    subplot(iNumChans, 1, iChanCnt)
    plot(mAudioData(:,iChanCnt));
    axis([0, iNumSampsPerChan, -1, 1])
end

audiowrite('InterleavedNoiseBurst24Ch.wav', mAudioData, iFs);

%--------------------Licence ---------------------------------------------
% Copyright (c) <2011-2017> Hagen Jaeger
% Institute for Hearing Technology and Audiology
% Jade University of Applied Sciences
% All rights reserved.

% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in 
%    the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its 
%    contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
% COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
% USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.