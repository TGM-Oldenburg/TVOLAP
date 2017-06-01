% This script evaluates an partitioned overlap add method.
% Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
% Version History:
% Ver. 0.1 initial create (empty) 08.11.2016 (HJ)
% Ver. 1.0 debugged and tested 08.11.2016 (HJ)
% Ver. 2.0 Encapsulated the processing routine in WOLAP class (HJ)

clear; close all; clc

set(0,'defaulttextinterpreter','latex');

fs = 48000;
timeDurAudio = 2;
timeDurIR = 0.2;

blockLen = 1024;
processLen = 2*blockLen;
nfft = 2*processLen;

numSampsIR = round(timeDurIR*fs/2)*2;
tau = timeDurIR/10;
alpha = 1-exp(-1/(tau*fs));
fzero = 1;
stepFzero = 1;
numChans = 2;
numIR = 10;
numImpulses = 10;
switchFreq = 10;
plotChan = 2;
plotIR = 1;

impResp = zeros(numSampsIR, numChans, numIR);
for cntIR=1:numIR
    fzero = fzero + stepFzero;
    for chanCnt=1:numChans
        impRespTmp = zeros(numSampsIR,1);
        impRespTmp(numSampsIR*0.5+1:end) = ...
            (-sin(fzero*linspace(0, 2*pi, numSampsIR*0.5))').*exp(-alpha*(1:numSampsIR*0.5)');
        impRespTmp = impRespTmp/sum(abs(impRespTmp));
        impResp(:,chanCnt, cntIR) = impRespTmp;
    end
end

numSampsAudio = round(timeDurAudio*fs);
inData = zeros(numSampsAudio, numChans);
inData(floor(linspace(10, numSampsAudio-10, numImpulses)/blockLen)*blockLen+1,:) = 1;
outData = zeros(size(inData));

numBlocks = floor(numSampsAudio/blockLen)-1;

wolapInst = WOLAP(impResp, blockLen);
IR_ID = 1;

tic;
for blockCnt=1:numBlocks
    tmpIdxLo = (blockCnt-1)*blockLen+1;
    tmpIdxHi = tmpIdxLo+blockLen-1;
    outData(tmpIdxLo:tmpIdxHi,:) = wolapInst.process(inData(tmpIdxLo:tmpIdxHi,:));
    if mod(blockCnt,switchFreq) == 0
        IR_ID = mod(IR_ID,numIR)+1;
        wolapInst.setIR(IR_ID);
    end
end
disp(['WOLAP needed ' , num2str(toc*1000) , 'ms to process.']);

outDataComp = inData;
tic;
for chanCnt=1:numChans
    outDataComp(:, chanCnt) = filter(impResp(:,chanCnt,1), 1, inData(:,chanCnt));
end
outDataComp = outDataComp(1:numSampsAudio,:);
disp(['filter(b,1,x) needed ' , num2str(toc*1000) , 'ms to process.']);

dims = [2, 2, 20, 13.8];
hFig = figure('units', 'centimeters', 'position', dims);

inData = inData(:,plotChan);
outData = outData(:,plotChan);
impResp = impResp(:,plotChan, plotIR);
outDataComp = outDataComp(:,plotChan);

ax1 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 11.1, 18.4, 2]);
plot(ax1, 1/fs:1/fs:numSampsAudio/fs, inData, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, timeDurAudio-0.1, min(inData)*1.2, max(inData)*1.2]);
title('Input signal x(t)');

ax2 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 7.8, 18.4, 2]);
plot(ax2, 1/fs:1/fs:numSampsIR/fs, impResp, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, timeDurIR, min(impResp)*1.2, max(impResp)*1.2]);
title('Impulse response h(t)');

ax3 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 4.5, 18.4, 2]);
plot(ax3, 1/fs:1/fs:numSampsAudio/fs, outDataComp, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, timeDurAudio-0.1, min(outDataComp)*1.2, max(outDataComp)*1.2]);
title('Output signal y(t) via direct conv.'); ylabel('\hspace{3cm} Amplitude A [NV]');

ax4 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 1.2, 18.4, 2]);
plot(ax4, 1/fs:1/fs:numSampsAudio/fs, outData, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, timeDurAudio-0.1, min(outData)*1.2, max(outData)*1.2]);
title('Input signal x(t)');
title('Output signal y(t) via conv. engine');
xlabel('Time t [s]');

set(hFig, 'PaperPosition', [0, 0, dims(3:4)]); 
set(hFig, 'PaperSize', dims(3:4))
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