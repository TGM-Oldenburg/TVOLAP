% This script evaluates an partitioned overlap add method.
% Author: Hagen Jaeger (c) TGM @ Jade Hochschule applied licence see EOF
% Version History:
% Ver. 0.1 initial create (empty) 08.11.2016 (HJ)
% Ver. 1.0 debugged and tested 08.11.2016 (HJ)
% Ver. 2.0 Encapsulated the processing routine in TVOLAP class (HJ)

clear; close all; clc

set(0,'defaulttextinterpreter','latex');

fs = 48000;
timeDurAudio = 5;
timeDurIR = 0.2;
numChans = 4;
numChansWav = 2;
numIR = 20;
blockLen = 512;

impResp = zeros(round(fs*timeDurIR), numChans);
deltaImp = [1; zeros(round(fs*timeDurIR*2)/2-1,1)];
deltaImp = fft(deltaImp);
deltaImp = deltaImp(1:end/2+1);

freqCuts = round(logspace(log10(0.01),log10(0.99),numChans+1)*size(deltaImp,1));
weight = sqrt(freqCuts(2)/freqCuts(1)
for cnt = 1:numChans
    deltaImpCut = numChans*deltaImp./weight^(cnt-1);
    deltaImpCut(1:freqCuts(cnt)) = 0.0;
    deltaImpCut(freqCuts(cnt+1):end) = 0.0;
    impResp(:,cnt) = real(ifft([deltaImpCut; conj(deltaImpCut(end-1:-1:2))]));
end

tmpSize = size(impResp,1)/2;
impResp = [impResp(tmpSize+1:end,:); impResp(1:tmpSize,:)];
impResp = repmat(impResp, 1, 1, numIR);
tmpIdx = 1:numChans;
for cnt = 1:numIR;
    impResp(:, :, cnt) = impResp(:, tmpIdx, cnt);
    tmpIdx = [tmpIdx(2:end), tmpIdx(1)];
end

inSig = randn(round(fs*timeDurAudio/blockLen)*blockLen, numChans)*0.1;
specLen = size(inSig,1)/2+1;
weight = repmat(sqrt(1./(1:specLen)*specLen)',1,numChans);
inSig = fft(inSig);
inSig = inSig(1:end/2+1,:).*weight;
inSig = real(ifft([inSig; conj(inSig(end-1:-1:2,:))]));

outSig = inSig;
TVOLAPInst = TVOLAP(impResp, blockLen);
IR_ID = 1;

tic;
for blockCnt = 1:floor(fs*timeDurAudio/blockLen)
    tmpIdxLo = (blockCnt-1)*blockLen+1;
    tmpIdxHi = tmpIdxLo+blockLen-1;
    outSig(tmpIdxLo:tmpIdxHi,:) = TVOLAPInst.process(inSig(tmpIdxLo:tmpIdxHi,:));
    if mod(blockCnt,50) == 49
        IR_ID = mod(IR_ID,numIR)+1;
        TVOLAPInst.setIR(IR_ID);
    end
end
disp(['TVOLAP needed ' , num2str(toc*1000) , 'ms to process.']);

outSigComp = inSig;
tic;
for chanCnt=1:numChans
    outSigComp(:, chanCnt) = filter(impResp(:,chanCnt,1), 1, inSig(:,chanCnt));
end
disp(['filter(b,1,x) needed ' , num2str(toc*1000) , 'ms to process.']);

if numChans>1
    audiowrite('outTst.wav', outSig(:,1:numChansWav), fs)
else
    audiowrite('outTst.wav', outSig, fs)
end

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