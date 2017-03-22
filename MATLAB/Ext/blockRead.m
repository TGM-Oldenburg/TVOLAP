classdef blockRead < handle
    %Audio read class
    
    properties (GetAccess = 'public', SetAccess = 'private')
        sInFile;
        sInPath;
        iInFid;
        bIsInitialized = false;
        iFs;
        iBlockLen;
        iNchans;
        iDataLen;
        iBps;
        iFrameSize;
        iBitsPerSample;
        iBytesPerBlock;
        iReadLen;
        iHeaderLen;
        sReadFmt;
        fNormFact;
        fNormFactInv;
        iAudioBytes;
        iByteCnt;
    end
    
    methods
        function obj = blockRead(sInFileInt, sInPathInt, iBlockLen)
            
            if nargin < 1
                error('wavRW:Constructor:InputErr', ...
                    'At least you have to enter a audio filename (current folder).');
            end
            if nargin > 4 
                error('wavRead:Constructor:InputErr', 'Too many input arguments'); 
            end
            
            if ~strcmp(sInFileInt(end-3:end), '.wav') && ~strcmp(sInFileInt(end-3:end), '.mp3')
                sInFileInt = [sInFileInt '.wav']; 
            end
            
            if nargin < 2, sInPathInt = pwd; end
            
            if ~exist(fullfile(sInPathInt,sInFileInt), 'file');
                error('AudioIO:FileNotFound', ...
                    'Input file %s does not exist. Wrong Filename, -path?', ...
                    fullfile(sInPathInt,sInFileInt));
            end
            
            obj.sInFile = sInFileInt;
            obj.sInPath = sInPathInt;
            
            if strcmpi(sInFileInt(end-3:end),'.mp3')
                if ~exist(fullfile(obj.sInPath,[obj.sInFile(1:end-4) '.wav']), 'file')
                    obj.convMp3ToWav(obj.sInFile, obj.sInPath);
                else
                    obj.sInFile = [obj.sInFile(1:end-4) '.wav'];
                end
            end
           
            if nargin < 3, iBlockLen = 1024; disp('BlockLen set to 1024'); end
            
            obj.iBlockLen = iBlockLen;
                
            obj.iInFid = fopen(fullfile(obj.sInPath,obj.sInFile), 'r');
            wavReadHeader(obj);
            
            obj.iFrameSize = floor((obj.iBitsPerSample+7)/8)*obj.iNchans;
            obj.iBps = obj.iFs*obj.iFrameSize;
            obj.sReadFmt = ['bit' num2str(obj.iBitsPerSample)];
            obj.iReadLen = obj.iBlockLen*obj.iNchans;
            obj.iHeaderLen = 44;

            obj.fNormFact = 2^(-obj.iBitsPerSample+1);
            obj.fNormFactInv = 2^(obj.iBitsPerSample-1);
            obj.iByteCnt = 0;
            obj.iBytesPerBlock = obj.iBlockLen*obj.iFrameSize;
                
            obj.bIsInitialized = true;
        end
        
        function delete(obj)
            if obj.bIsInitialized
                fclose(obj.iInFid);
            end
        end
        
        function wavReadHeader(obj)
            if ~strcmp('RIFF', fread(obj.iInFid, 4, '*char')')
                error('wavReadHeader:RIFF', ...
                    'RIFF header not found. Wrong or corrupted file.');
            end
            
            obj.iDataLen = fread(obj.iInFid, 1, 'uint32');
                        
            if ~strcmp('WAVE', fread(obj.iInFid, 4, '*char')')
                error('wavReadHeader:WAVE', ...
                    'WAVE header not found. Wrong or corrupted file.');
            end
            
            tmp = '____';
            fileIdx = 12;
            
            while (~strcmp(' tmf', tmp) && fileIdx < obj.iDataLen)
                fileIdx = fileIdx+1;
                tmp(2:4) = tmp(1:3);
                tmp(1) = fread(obj.iInFid, 1, '*char');
            end
            
            if fileIdx >= obj.iDataLen
                error('wavReadHeader:fmt', ...
                    'fmt header not found. Wrong or corrupted file.')
            end
            
            tmp = fread(obj.iInFid, 1, 'uint32');
            
            if tmp < 16
                error('wavReadHeader:fmtLen', ...
                    'Format section too short. Unsupported or corrupted file.');
            elseif tmp > 16
                error('wavReadHeader:fmtLen', ...
                    'Format section too long. Unsupported WAVE file format.');
            end
            
            if fread(obj.iInFid, 1, 'int16') ~= 1
                error('wavReadHeader:PCM', ...
                    'This is no linear PCM .wav file.');
            end
            
            obj.iNchans = fread(obj.iInFid, 1, 'uint16');
            
            obj.iFs = fread(obj.iInFid, 1, 'uint32');
            
            obj.iBps = fread(obj.iInFid, 1, 'uint32');
            
            obj.iFrameSize = fread(obj.iInFid, 1, 'uint16');
            
            obj.iBitsPerSample = fread(obj.iInFid, 1, 'uint16');
            
            obj.sReadFmt = ['bit' num2str(obj.iBitsPerSample)];
            
            tmp = '____';
            
            while (~strcmp('atad', tmp) && fileIdx < obj.iDataLen)
                fileIdx = fileIdx+1;
                tmp(2:4) = tmp(1:3);
                tmp(1) = fread(obj.iInFid, 1, '*char');
            end
            
            if fileIdx >= obj.iDataLen
                error('wavReadHeader:data', ...
                    'data header not found. Unsupported or corrupted file.')
            end
            
            obj.iAudioBytes = fread(obj.iInFid, 1, 'uint32');
            
            obj.fNormFact = 2^(-obj.iBitsPerSample+1);
            
            obj.iReadLen = obj.iNchans*obj.iBlockLen;
            
            obj.iHeaderLen = ftell(obj.iInFid);
        end
        
        function mIn = readOneBlock(obj)
            mIn = reshape(fread(obj.iInFid, obj.iReadLen, obj.sReadFmt)', ...
                obj.iNchans, [])'*obj.fNormFact;
            obj.iByteCnt = obj.iByteCnt + obj.iBytesPerBlock;
            
            if obj.iByteCnt >= obj.iAudioBytes-obj.iBytesPerBlock
                obj.iByteCnt = obj.iBytesPerBlock;
                fseek(obj.iInFid, 44+obj.iBytesPerBlock, -1);
            end
        end
        
        function iFs = getSamplingRate(obj)
            iFs = obj.iFs;
        end
        
        function iNchans = getNrOfChans(obj)
            iNchans = obj.iNchans;
        end
        
        function setPositionNorm(obj, fNewPosition)
            if fNewPosition>1
                fNewPosition = 1;
                warning('wavRead:setPositionInPercent', 'invalid position input');
            elseif fNewPosition<0
                fNewPosition = 0;
                warning('wavRead:setPositionInPercent', 'invalid position input');
            end
            
            obj.iByteCnt = obj.iHeaderLen+obj.iAudioBytes*fNewPosition;
            obj.iByteCnt = round(obj.iByteCnt/obj.iFrameSize)*obj.iFrameSize;
            if obj.iByteCnt >= (obj.iAudioBytes-obj.iBytesPerBlock)
                obj.iByteCnt = obj.iAudioBytes-obj.iBytesPerBlock;
            end
            fseek(obj.iInFid, obj.iByteCnt, -1);
        end
        
        function fPosition = getPositionNorm(obj)
            fPosition = (obj.iByteCnt-obj.iHeaderLen)/obj.iAudioBytes;
        end
        
        function convMp3ToWav(obj, inFileName, inFilePath, outFileName, outFilePath)
        %convert a mp3 file
            if nargin < 1
                error('Mp3ToWav:InvalidInput', ...
                    'At least, you have to enter the name of the mp3 file (current dir)');
            end
            if length(inFileName) < 5, inFileName = [inFileName '.mp3']; end
            if inFileName(end-3) ~= '.', inFileName = [inFileName '.mp3']; end
            if ~strcmpi(inFileName(end-3:end),'.mp3')
                error('Mp3ToWav:InvalidFile', ...
                'the input file is not a mp3-file. Wrong file ending.');
            end
            if nargin < 3, inFilePath = pwd; end
            if nargin < 4, outFileName = [inFileName(1:end-4) '.wav']; end
            if nargin < 5, outFilePath = pwd; end
            if ~exist(fullfile(inFilePath, inFileName), 'file')
                error('Mp3ToWav:FileNotFound', ...
                    ['The inputted file' char(10) fullfile(inFilePath, inFileName)...
                    char(10) 'does not exist.'])
            end
            
            if ispc
                cmd = ['"' fullfile('Ext','mpg123.exe') '" -q -w "' ...
                    fullfile(outFilePath,outFileName) '" "'...
                    fullfile(inFilePath,inFileName) '"'];
            elseif ismac
                cmd = ['"' fullfile('Ext','mpg123.maci64') '" -q -w "' ...
                fullfile(outFilePath,outFileName) '" "'...
                fullfile(inFilePath,inFileName) '"'];
            else
                cmd = ['"' fullfile('Ext','mpg123.gnlx86') '" -q -w "' ...
                fullfile(outFilePath,outFileName) '" "'...
                fullfile(inFilePath,inFileName) '"'];
            end
            system(cmd);
            
           obj.sInFile = outFileName;
           obj.sInPath = outFilePath;
           
        end
    end   
end

