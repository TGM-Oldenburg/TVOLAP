classdef WOLAP < handle
%--------------------------------------------------------------------------
% Real time head related transfer function processing with special 
% time-variant implementation of the partitioned convolution in frequency 
% domain (to prevent from perceptive noticeable audio artifacts).                                 
%                                                                              
% Author: (c) Hagen Jaeger                      April 2016 - March 2017   
% LGPL Release: May 2017, License see end of file                              
%--------------------------------------------------------------------------
    properties(GetAccess=public, SetAccess=protected)
        blockLen;
        processLen;
        nfft;
        numParts;
        numFreqMems;
        numChans;
        numIR;
        freqSaveCnt;
        modCnt;
        win;
        freqIR;
        freqMem;
        inData;
        convMem;
        outDataMem;
        freqSumReset;
        inDataProcess;
        IR_ID;
    end
    
    methods
        function obj = WOLAP(impResp, blockLen)
            if ndims(impResp) ~= 3
                error(['impulse response is in wrong format. Has to be ' ...
                    '<length of one IR> x <channels> x <number of IRs>.'])
            end
            
            if blockLen ~= 2^round(log2(blockLen))
                error('block length has to be element of 2^n.');
            end
            
            obj.blockLen = blockLen;
            obj.processLen = obj.blockLen*2;
            obj.nfft = 2*obj.processLen;
            obj.numIR = size(impResp, 3);
            obj.numChans = size(impResp,2);
            
            obj.numParts = ceil(size(impResp,1)/obj.processLen);
            obj.numFreqMems = 2*obj.numParts;
            numSampsIR = obj.numParts*obj.processLen;
            impResp(end:numSampsIR,:,:) = 0;

            obj.convMem = zeros(obj.processLen,obj.numChans,2);
            obj.outDataMem = zeros(obj.blockLen,obj.numChans);

            obj.win = 0.5+0.5*cos(linspace(-pi,pi,obj.processLen+1))';
            obj.win = repmat(obj.win(1:end-1),1,obj.numChans);
            
            obj.freqIR = zeros(obj.processLen+1, obj.numChans,obj.numParts, obj.numIR);
            for cntIR=1:obj.numIR
                for partCnt=1:obj.numParts
                    for chanCnt=1:obj.numChans
                        tmpIdxLo = (partCnt-1)*obj.processLen+1;
                        tmpIdxHi = tmpIdxLo+obj.processLen-1;
                        tmpVec = fft(impResp(tmpIdxLo:tmpIdxHi, chanCnt, cntIR), obj.nfft);
                        obj.freqIR(:,chanCnt, partCnt, cntIR) = tmpVec(1:obj.processLen+1);
                    end
                end
            end

            obj.freqMem = zeros(obj.processLen+1, obj.numChans, obj.numFreqMems);
            
            obj.freqSumReset = zeros(obj.processLen+1, obj.numChans);
            
            obj.freqSaveCnt = 1;
            obj.modCnt = 1;
            obj.outDataMem = zeros(obj.blockLen, obj.numChans);
            obj.inDataProcess = zeros(obj.processLen, obj.numChans);
            obj.IR_ID = 1;
        end
        
        function data = process(obj, data)
            obj.inDataProcess(1:obj.blockLen,:) = obj.inDataProcess(obj.blockLen+1:end,:);
            obj.inDataProcess(obj.blockLen+1:end,:) = data;
            tmpVec = fft(obj.inDataProcess.*obj.win, obj.nfft);
            obj.freqMem(:,:,obj.freqSaveCnt) = tmpVec(1:obj.processLen+1,:);

            freqSum = obj.freqSumReset;

            freqReadCnt = obj.freqSaveCnt;
            for partCnt = 1:obj.numParts
                freqSum = freqSum + obj.freqMem(:,:,freqReadCnt).*obj.freqIR(:,:,partCnt,obj.IR_ID);

                freqReadCnt = freqReadCnt-2;
                if freqReadCnt<1, freqReadCnt = freqReadCnt+obj.numFreqMems; end
            end

            tmpVec = real(ifft([freqSum ; flipud(conj(freqSum(2:end-1,:)))], obj.nfft));
            outDataProcess = tmpVec(1:obj.processLen,:)+obj.convMem(:,:,obj.modCnt);
            obj.convMem(:,:,obj.modCnt) = tmpVec(obj.processLen+1:end,:);

            data = outDataProcess(1:obj.blockLen,:)+obj.outDataMem;
            obj.outDataMem = outDataProcess(obj.blockLen+1:end,:);

            obj.modCnt = obj.modCnt+1;
            obj.freqSaveCnt = obj.freqSaveCnt+1;
            if obj.modCnt>2, obj.modCnt = 1; end
            if obj.freqSaveCnt>obj.numFreqMems, obj.freqSaveCnt=1; end 
        end
        
        function setIR(obj, IR_ID)
            if (IR_ID >= 1) && (IR_ID <= obj.numIR)
                IR_ID = round(IR_ID);
                obj.IR_ID = IR_ID;
            else
               error('requested impulse response ID is out of range.');
            end
        end
    end
end

%------------------------------License-------------------------------------
% Copyright (c) 2012-2017 Hagen Jaeger, Uwe Simmer                            
%                                                                             
% This program is free software: you can redistribute it and/or modify        
% it under the terms of the GNU Lesser General Public License as published 
% by  the Free Software Foundation, either version 3 of the License, or           
% (at your option) any later version.                                         
%                                                                             
% This program is distributed in the hope that it will be useful,             
% but WITHOUT ANY WARRANTY; without even the implied warranty of              
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               
% GNU Lesser General Public License for more details.                         
%                                                                             
% You should have received a copy of the GNU Lesser General Public License    
% along with this program. If not, see <http://www.gnu.org/licenses/>.        
%--------------------------------------------------------------------------