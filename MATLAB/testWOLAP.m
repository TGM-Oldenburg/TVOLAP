% This script evaluates an partitioned overlap add method.
% Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
% Version History:
% Ver. 0.1 initial create (empty) 08.11.2016 HJ
% Ver. 1.0 debugged and tested 08.11.2016 HJ

clear; close all; clc

set(0,'defaulttextinterpreter','latex');

iFs = 48000;
fTimeDurAudio = 1.2;
fTimeDurIR = 0.4;

iBlockLen = 1024;
iProcessLen = 2*iBlockLen;
iNfft = 2*iProcessLen;

iNumSampsIR = round(fTimeDurIR*iFs/2)*2;
fTau = fTimeDurIR/10;
fAlpha = 1-exp(-1/(fTau*iFs));
iFreq = 5;
vIR = sin(iFreq*linspace(-pi, pi, iNumSampsIR*0.5))';
vIR = vIR.*exp(-fAlpha*(1:iNumSampsIR*0.5)');
vIR = [zeros(iNumSampsIR/4,1); vIR];
vIR = vIR/sum(abs(vIR));

iNumParts = ceil(iNumSampsIR/iProcessLen);
iNumFreqMems = 2*iNumParts;
iNumSampsIR = iNumParts*iProcessLen;
vIR(end:iNumSampsIR) = 0;

iNumSampsAudio = round(fTimeDurAudio*iFs);
iFreq = 3;
vIn = square(2*pi*iFreq*((1:iNumSampsAudio)'/iFs));
vOut = zeros(size(vIn));

iNumBlocks = floor(iNumSampsAudio/iBlockLen)-1;

vInProcess = zeros(iProcessLen,1);
vConvMem = zeros(iProcessLen,2);
vOutMem = zeros(iBlockLen,1);

vWin = 0.5+0.5*cos(linspace(-pi,pi,iProcessLen+1))';
vWin = vWin(1:end-1);

mFreqIR = zeros(iProcessLen+1,iNumParts);
for iPartCnt=1:iNumParts
    iTmpIdxLo = (iPartCnt-1)*iProcessLen+1;
    iTmpIdxHi = iTmpIdxLo+iProcessLen-1;
    vTmp = fft(vIR(iTmpIdxLo:iTmpIdxHi), iNfft);
    mFreqIR(:,iPartCnt) = vTmp(1:iProcessLen+1);
end

mFreq = zeros(iProcessLen+1, iNumFreqMems);

tic;

iFreqSaveCnt = 1;
iModCnt = 1;
for iBlockCnt=1:iNumBlocks
    iTmpIdxLo = (iBlockCnt-1)*iBlockLen+1;
    iTmpIdxHi = iTmpIdxLo+iBlockLen-1;
    
    vInProcess(1:iBlockLen) = vInProcess(iBlockLen+1:end);
    vInProcess(iBlockLen+1:end) = vIn(iTmpIdxLo:iTmpIdxHi);
    
    vTmp = fft(vInProcess.*vWin, iNfft);
    mFreq(:,iFreqSaveCnt) = vTmp(1:iProcessLen+1);
    
    vFreqSum = zeros(iProcessLen+1,1);
    
    iFreqReadCnt = iFreqSaveCnt;
    for iPartCnt = 1:iNumParts
        
        vFreqSum = vFreqSum + mFreq(:,iFreqReadCnt).*mFreqIR(:,iPartCnt);
        
        iFreqReadCnt = iFreqReadCnt-2;
        if iFreqReadCnt<1, iFreqReadCnt = iFreqReadCnt+iNumFreqMems; end
    end
        
    vTmp = ifft([vFreqSum ; flipud(conj(vFreqSum(2:end-1)))], iNfft, 'symmetric');
    vOutProcess = (vTmp(1:iProcessLen)+vConvMem(:,iModCnt));
    vConvMem(:,iModCnt) = vTmp(iProcessLen+1:end);
    
    vOut(iTmpIdxLo:iTmpIdxHi) = vOutProcess(1:iBlockLen)+vOutMem;
    vOutMem = vOutProcess(iBlockLen+1:end);
    
    iModCnt = iModCnt+1;
    iFreqSaveCnt = iFreqSaveCnt+1;
    if iModCnt>2, iModCnt = 1; end
    if iFreqSaveCnt>iNumFreqMems, iFreqSaveCnt=1; end
end

disp(['WOLAP needed ' , num2str(toc*1000) , 'ms to process.']);

tic;

vOutComp = filter(vIR, 1, vIn);
vOutComp = vOutComp(1:iNumSampsAudio);

disp(['filter(b,1,x) needed ' , num2str(toc*1000) , 'ms to process.']);

vDims = [2, 2, 8, 13.8];
hFig = figure('units', 'centimeters', 'position', vDims);

hAx1 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 11.1, 6.4, 2]);
plot(hAx1, 1/iFs:1/iFs:iNumSampsAudio/iFs, vIn, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, fTimeDurAudio-0.1, min(vIn)*1.2, max(vIn)*1.2]);
title('Input signal x(t)');

hAx2 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 7.8, 6.4, 2]);
plot(hAx2, 1/iFs:1/iFs:iNumSampsIR/iFs, vIR, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, fTimeDurIR, min(vIR)*1.2, max(vIR)*1.2]);
title('Impulse response h(t)');

hAx3 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 4.5, 6.4, 2]);
plot(hAx3, 1/iFs:1/iFs:iNumSampsAudio/iFs, vOutComp, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, fTimeDurAudio-0.1, min(vOutComp)*1.2, max(vOutComp)*1.2]);
title('Output signal y(t) via direct conv.'); ylabel('\hspace{3cm} Amplitude A [NV]');

hAx4 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 1.2, 6.4, 2]);
plot(hAx4, 1/iFs:1/iFs:iNumSampsAudio/iFs, vOut, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, fTimeDurAudio-0.1, min(vOut)*1.2, max(vOut)*1.2]);
title('Input signal x(t)');
title('Output signal y(t) via conv. engine');
xlabel('Time t [s]');

set(hFig, 'PaperPosition', [0, 0, vDims(3:4)]); 
set(hFig, 'PaperSize', vDims(3:4))
print(hFig,'-dpdf','-r600','testWOLAP.pdf');

%--------------------Licence ---------------------------------------------
% Copyright (c) <2011-2016> Hagen Jaeger
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