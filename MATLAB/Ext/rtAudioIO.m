classdef rtAudioIO < handle
    %Audio In-Output class
    
    properties (GetAccess = 'public', SetAccess = 'private')
        sInDevName;
        iInDevNr;
        sOutDevName;
        iOutDevNr;
        iFs;
        iBlockLen;
        iNchans;
        mSoundInfo;
        mSoundInputsIdx;
        mSoundOutputsIdx;
    end
    
    methods (Access = public)
        function obj = rtAudioIO(iInDevNr, iOutDevNr, iBlockLen, iFs, iNchans)
            
            tmp = msound('verbose', 0);

            if tmp.IsOpenForReading || tmp.IsOpenForWriting
                msound('close');
            end

            if nargin < 1, iInDevNr = 0; end
            if nargin < 2, iOutDevNr = 0; end
            if nargin < 3, iBlockLen = 1024; disp('BlockLen set to 1024'); end
            if nargin < 4, iFs = 44100; disp('Sampling frequency set to 44100Hz'); end
            if nargin < 5, iNchans = 1; disp('Nr of chans set to 1'); end
            
            obj.mSoundInfo = msound('deviceInfo');
            
            obj.sInDevName = 'N/A';
            obj.sOutDevName = 'N/A';
            
            if sum([obj.mSoundInfo.inputs]) < 1
                iInDevNr = [];
            end
            if sum([obj.mSoundInfo.outputs]) < 1
                iOutDevNr = [];
            end
            
            if isempty(iInDevNr) && isempty(iOutDevNr)
                error('rtAudioIO:Constructor:NoDevices', ...
                    'No audio devices avilable. Init aborted.');
            end
            
            obj.mSoundInputsIdx = find([obj.mSoundInfo.inputs] > 0);
            obj.mSoundOutputsIdx = find([obj.mSoundInfo.outputs] > 0);
            
            if iInDevNr ~= 0
                obj.iInDevNr = obj.mSoundInputsIdx(iInDevNr);
            else
                obj.iInDevNr = iInDevNr;
            end
            
            if iOutDevNr ~= 0
                obj.iOutDevNr = obj.mSoundOutputsIdx(iOutDevNr);
            else
                obj.iOutDevNr = iOutDevNr;
            end
            
            obj.iFs = iFs;
            obj.iBlockLen = iBlockLen;
            obj.iNchans = iNchans;

            obj.deviceChanger();
        end
        
        function delete(obj)
            obj.iFs = 0;
            msound('close');
        end
        
        function sInDevName = changeInDevice(obj, iInDevNr)
            if iInDevNr < 0
                warning(['Input device number cannot be smaller than zero. ' ...
                    'Changed to zero (default device).'])
                iInDevNr = 0;
            end
            
            if iInDevNr > length(obj.mSoundInfo)
                warning(['Invalid Input device number (too big). ' ...
                    'Changed to zero (default device).'])
                iInDevNr = 0;
            end
            
            if iInDevNr ~= 0
                obj.iOutDevNr = obj.mSoundOutputsIdx(iInDevNr);
            end
            
            obj.deviceChanger();
            
            sInDevName = obj.sInDevName;
        end
        
        function sOutDevName = changeOutDevice(obj, iOutDevNr)
            if iOutDevNr < 0
                warning(['Input device number cannot be smaller than zero. ' ...
                    'Changed to zero (default device).'])
                iOutDevNr = 0;
            end
            
            if iOutDevNr > length(obj.mSoundInfo)
                warning(['Invalid Input device number (too big). ' ...
                    'Changed to zero (default device).'])
                iOutDevNr = 0;
            end
            
            if iOutDevNr ~= 0
                obj.iOutDevNr = obj.mSoundOutputsIdx(iOutDevNr);
            end
            
            obj.deviceChanger();
            
            sOutDevName = obj.sOutDevName;
        end
        
        function [cInDevNames, cOutDevNames] = getAllDeviceNames(obj)    
            if ~isempty(obj.iInDevNr)
                cInDevNames = {obj.mSoundInfo(obj.mSoundInputsIdx).name};
                cInDevNames = strcat(': ', cInDevNames);
                cInDevNames = strcat({obj.mSoundInfo(obj.mSoundInputsIdx).api}, cInDevNames);
            else
                cInDevNames = 'N/A';
            end
            
            if ~isempty(obj.iOutDevNr)
                cOutDevNames = {obj.mSoundInfo(obj.mSoundOutputsIdx).name};
                cOutDevNames = strcat(': ', cOutDevNames);
                cOutDevNames = strcat({obj.mSoundInfo(obj.mSoundOutputsIdx).api}, cOutDevNames);
            else
                cOutDevNames = 'N/A';
            end
        end
        
    end
    
    methods (Access = private)

        function deviceChanger(obj)
            if ~isempty(obj.iInDevNr)
                if obj.iInDevNr > 0, obj.sInDevName = obj.mSoundInfo(obj.iInDevNr).name; 
                else obj.sInDevName = 'Default'; end
            end
            
            if ~isempty(obj.iOutDevNr)
                if obj.iOutDevNr > 0, obj.sOutDevName = obj.mSoundInfo(obj.iOutDevNr).name; 
                else obj.sOutDevName = 'Default'; end
            end
            
            if ~isempty(obj.iInDevNr) && ~isempty(obj.iOutDevNr)
                msound('openRW', [obj.iInDevNr obj.iOutDevNr], obj.iFs, ...
                    obj.iBlockLen, obj.iNchans);
            elseif ~isempty(obj.iInDevNr)
                msound('openRead', obj.iInDevNr, obj.iFs, obj.iBlockLen, obj.iNchans);
            elseif ~isempty(obj.iOutDevNr)
                msound('openWrite', obj.iOutDevNr, obj.iFs, obj.iBlockLen, obj.iNchans);
            end
        end
    end
    
    methods(Static = true)
        function mIn = readOneBlock() 
            mIn = msound('getSamples'); 
        end
        
        function writeOneBlock(mOut) 
            msound('putSamples', mOut); 
        end
    end
end