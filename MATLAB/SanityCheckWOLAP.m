% This script describes how to calculate the partitioned weighted overlap 
% add (WOLAP) method and was written for testing purposes (better / easier 
% visualization tools and debugging as in C++ code).
% Author: H. Jaeger (c) IHA @ Jade Hochschule applied licence see EOF 
% Version History:
% Ver. 0.01 initial create (empty) 17.01.2017 HJ
% Ver. 1.0 seems to work  21.01.2017 HJ

clear; close all; clc;
addpath('Ext')

%% feel free to change
sAudioName = 'HotelCaliforniaSnippet.wav'; % type in the audio file name
iAudioDeviceNr = 0; %0 or [] means standard output device
fChangeTime = 0.05; %change impulse response every fChangeTime seconds
iProcessLen = 256; %make four parts for WOLAP

%% preallocate some useful variables
iOverlapFact = 2; %50% overlap, define as 2^x (we recommend no change)
iBlockLen = iProcessLen/iOverlapFact; %audio block length
iNfft = 2*iProcessLen; %fft length of aucio (normally 2* blocklength for OLA)

%% initialize real time (blockoriented) sound in- / output (wav and mp3)
audioIn = blockRead(sAudioName, pwd, iBlockLen); %wav/mp3 read
audioOut = rtAudioIO([], iAudioDeviceNr, iBlockLen, audioIn.iFs, audioIn.iNchans);

%% generate syntethic time varying impulse responses (IIR->FIR peak filters)
fGain = -30;
fQualityFact = 20.0;
fLoFreq = 100;
fHiFreq = 20000;
iNumIR = 80;
iImpRespLen = 8192;
vCenterFreq = logspace(log10(fLoFreq), log10(fHiFreq), iNumIR+2);
vCenterFreq = [vCenterFreq, fliplr(vCenterFreq), vCenterFreq, fliplr(vCenterFreq)];
iNumIR = 2*iNumIR;
mIR = zeros(iImpRespLen, audioIn.iNchans, iNumIR);
mIR(1000,:,:) = 1.0;
iNumNotches = 9;
for kk=1:iNumIR
    for ll=0:iNumNotches
        fOm=2*pi*vCenterFreq(kk+round(iNumIR/iNumNotches*ll))/audioIn.iFs;
        fAlpha=sin(fOm)*0.5/fQualityFact;
        fTmp1 = 10^(fGain*0.025);
        fTmp2 = 1+fAlpha/fTmp1;
        vB(1) = (1+fAlpha*fTmp1)/fTmp2;
        vB(2) = (-2*cos(fOm))/fTmp2;
        vB(3) = (1-fAlpha*fTmp1)/(1+fAlpha/fTmp1);
        vA(1) = 1.0;
        vA(2) = vB(2);
        vA(3) = (1-fAlpha/fTmp1)/(1+fAlpha/fTmp1);

        mIR(:,:,kk) = filter(vB, vA, mIR(:,:,kk));
    end
end

iNumParts = iImpRespLen/iProcessLen; %calculate number of frequency parts
iNumFreqMems = 2*iNumParts; %number of frequency vector states

%% calculate a periodic sqrt hann window with right number of channels:
vTimeAx = repmat(linspace(0,2*pi,iProcessLen+1)',1,audioIn.iNchans);
vWin = 0.5-0.5*cos(vTimeAx(1:end-1,:)); %sqrt hann and periodic window

%% some preallocation stuff
iChangeRate = floor(fChangeTime*audioIn.iFs/iBlockLen); %change rate of IR
iBlockCnt = 1;
iPosCnt = 1;
mInProcess = zeros(iProcessLen, audioIn.iNchans);
mOutProcess = zeros(iProcessLen, audioIn.iNchans);
mOutMem = zeros(iBlockLen, audioIn.iNchans);
mOut = zeros(iBlockLen, audioIn.iNchans);
mConvMem = zeros(iProcessLen, audioIn.iNchans, 2);
mInFreq = zeros(iProcessLen+1, audioIn.iNchans, iNumFreqMems);
mInFreqSum = zeros(iProcessLen+1, audioIn.iNchans);

%% divide the impulse response into partitions and make the fft of all
mImpRespParts = zeros(iProcessLen+1, audioIn.iNchans, iNumParts, iNumIR);
for iIRCnt = 1:iNumIR
    for iPartCnt = 1:iNumParts
        iTmpIdxLoTime = (iPartCnt-1)*iProcessLen+1;
        iTmpIdxHiTime = iTmpIdxLoTime+iProcessLen-1;
        iTmpIdxLoFreq = (iPartCnt-1)*iNfft+1;
        iTmpIdxHiFreq = iTmpIdxLoFreq+iNfft-1;
        mTmp = fft(mIR(iTmpIdxLoTime:iTmpIdxHiTime,:,iIRCnt), iNfft, 1);
        mImpRespParts(:, :, iPartCnt, iIRCnt) = mTmp(1:iProcessLen+1,:);
    end
end

%% some plotting stuff
hf = figure('units', 'normalized', 'position', [0.05, 0.05, 0.9, 0.7], ...
    'CloseRequestFcn', 'set(hf,''currentchar'',''b''); delete(hf); bClosed = true;');
set(hf,'currentchar','a'); %set initial char (needed for the loop break)
bClosed = false;
hAxTime = axes('parent', hf, 'units', 'normalized', 'position', ...
    [0.05 0.575 0.9 0.4], 'nextPlot', 'add');
hAxFreq = axes('parent', hf, 'units', 'normalized', 'position', ...
    [0.05 0.1 0.9 0.4], 'nextPlot', 'add');
mActImpRespPlot = fft(mIR(:,1,1));
mActImpRespPlot = 20*log10(abs(mActImpRespPlot(1:iImpRespLen*0.5+1)));
hpTime = plot(hAxTime, mIR(:,1,1));
set(hAxTime, 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLim', ...
    [1, iImpRespLen], 'YLim', [-0.1, 0.1], 'YGrid', 'on');
hpFreq = plot(hAxFreq, mActImpRespPlot);
set(hAxFreq, 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLim', ...
    [5, iImpRespLen*0.5], 'YLim', [2*fGain-2, 3], 'YGrid', 'on', 'XScale', 'log');

mActImpResp = mImpRespParts(:,:,:,1);
iFreqSaveCnt = 1;
iModCnt = 1;
bLoop = true;
%% run an infinite loop (break it on any keyboard event)
while bLoop
    %copy old audio from the back of the process block to the front(50% OL)
    mInProcess(1:iBlockLen,:) = mInProcess(iBlockLen+1:end,:);
    
    %read a new audio block, write it at the end of the process block
    mInProcess(iBlockLen+1:end,:) = audioIn.readOneBlock();
    
    %make the fft of the analyze windowed (sqrt hann) process block
    mTmp = fft(mInProcess.*vWin, iNfft, 1);
    mInFreq(:,:,iFreqSaveCnt) = mTmp(1:iProcessLen+1,:);
    
    %delete sum memory, run the loop over all partitions (complex mult.)
    mInFreqSum = zeros(iProcessLen+1, audioIn.iNchans);
    
    iFreqReadCnt = iFreqSaveCnt;
    for iPartCnt = 1:iNumParts
        mInFreqSum = mInFreqSum + mInFreq(:,:,iFreqSaveCnt).*mActImpResp(:,:,iPartCnt);
        
        iFreqReadCnt = iFreqReadCnt-iOverlapFact;
        if iFreqReadCnt<1, iFreqReadCnt = iFreqReadCnt+iNumFreqMems; end
    end
    
    %calculate the inverse fft after all parts were multiplied and added
    mTmp = ifft([mInFreqSum ; flipud(conj(mInFreqSum(2:end-1,:)))], iNfft, 1, 'symmetric');
    
    %apply the overlap add (new spectrum + old remainder spectrum)
    mOutProcess = mTmp(1:iProcessLen,:)+mConvMem(:,:,iModCnt);
    mConvMem(:,:,iModCnt) = mTmp(iProcessLen+1:end,:); %save the "new" remainder spectrum
    
    %overlap old and new audio block, sqrtHann*sqrtHann = hann window
    %and with 50% overlap, this gives perfect resynthesis.
    mOut = mOutProcess(1:iBlockLen,:)+mOutMem;
    %save the remaining half processing block to the overlap memory
    mOutMem = mOutProcess(iBlockLen+1:end,:);
    %playback the calculated audio block
    audioOut.writeOneBlock(mOut);
    
    %recalculate counters
    iModCnt = iModCnt+1;
    iFreqSaveCnt = iFreqSaveCnt+1;
    if iModCnt>2, iModCnt = 1; end
    if iFreqSaveCnt>iNumFreqMems, iFreqSaveCnt=1; end
    
    %only plotting stuff from here to the end
    drawnow();
    
    iBlockCnt = iBlockCnt+1;
    if iBlockCnt>iChangeRate
        iPosCnt = iPosCnt+1;
        if iPosCnt>iNumIR, iPosCnt = 1; end
        mActImpResp = mImpRespParts(:,:,:,iPosCnt);
        mActImpRespPlot = fft(mIR(:,1,iPosCnt));
        mActImpRespPlot = 20*log10(abs(mActImpRespPlot(1:iImpRespLen*0.5+1)));
        set(hpTime, 'yData', mIR(:,1,iPosCnt));
        set(hpFreq, 'yData', mActImpRespPlot);
        iBlockCnt = 1; 
    end
    
    if ~bClosed
        figure(hf);
        if ~strcmp(get(hf, 'currentchar'), 'a'), close(hf); end
    else
        bLoop = false;
    end
end

%--------------------Licence ---------------------------------------------
% Copyright (c) <2011-2017> Hagen Jaeger
% Institute for Hearing Technology and Audiology
% Jade University of Applied Sciences 
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject
% to the following conditions:
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.