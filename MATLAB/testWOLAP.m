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

disp(min(outData)*1.2)
disp(max(outData)*1.2)

ax4 = axes('parent', hFig, 'units', 'centimeters', 'position', [1.3, 1.2, 18.4, 2]);
plot(ax4, 1/fs:1/fs:numSampsAudio/fs, outData, 'linewidth', 2, 'color', [0.1, 0.1, 0.1]);
grid on; axis([0, timeDurAudio-0.1, min(outData)*1.2, max(outData)*1.2]);
title('Input signal x(t)');
title('Output signal y(t) via conv. engine');
xlabel('Time t [s]');

%set(hFig, 'PaperPosition', [0, 0, dims(3:4)]); 
%set(hFig, 'PaperSize', dims(3:4))
%print(hFig,'-dpdf','-r600','testWOLAP.pdf');

%--------------------Licence ---------------------------------------------
% Copyright (c) 2012-2017 Hagen Jaeger                           
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
% DEALINGS IN THE SOFTWARE.  