% =============================================================================
function [kSpace, varargout] = fastReadVD13D(tFilename,varargin)
% -----------------------------------------------------------------------------
% fastReadVD13D  Read a VD13D measurement data file with the property that all 
%                  lines for a given measurement have the same number of 
%                  samples.  This differs slightly from the VB reader 
%                  where only one measurement was stored per file.  VD can
%                  store multiple measurements.
%
% kSpace = fastReadVD13D() 
%    prompts user with gui and returns kSpace cell array, each cell storing
%               the kSpace from a particular measurement.
%
% [kSpace kOther] = fastReadVD13D(filename, ReadOther) 
%    reads filename and returns kSpace as well as kOther if it exists where
%       Other can be PhaseCor, RTFeedback, NoiseAdj.  These also are now
%       cell arrays, one element per measurement
%
% kSpace = fastReadVD13D(filename, option1, option2, option3, ....)
%    reads filename and returns kspace with options applied where
%    options may be:
%       ResetFFTScale: Turns off application of FFTScaling factors (default off)
%       ReadOneChannel, ReadOneChannelIndex - turn on and read one zero-based
%          indexed coil (default: OFF)
%       FillToFullFourier: Output array will have full FTlen for proper 
%          voxel sizes (default: ON)
%       ShiftDCToMatrixCenter : default OFF
% 
% To Invoke the ReadOther and option# options, define the spefic variable and 
%    assign the value 1 (ex. FillToFullFourier=1) before adding to function call.
% 
%>> FillToFullFourier = 1;
%>> kSpaceFull = fastReadVD13D(filename, FillToFullFourier); 
% -----------------------------------------------------------------------------

   % MDH constants
   ulMDH_ACQEND = 2^uint32(0);
   ulMDH_RTFEEDBACK = 2^uint32(1);
   ulMDH_SYNCDATA = 2^uint32(5);
   ulMDH_RAWDATACORRECTION = 2^uint32(10);
   ulMDH_PHASCOR= 2^uint32(21);
   ulMDH_NOISEADJSCAN = 2^uint32(25);
   ulMDH_REFLECT = 2^uint32(24);

   % Some default sizes
   ullMDHScanSize = uint64(192);
   ullMDHChanSize = uint64(32);
   ullMDHPMUSize = uint64(60);
   ullMaxBytesPerRead = 2^uint64(27);

   %% Definition of OFF and ON constants
   OFF = 0;
   ON = 1;

   %% Silent mode - no text output to stdout or gui
   sCon.sUser.bSilent = OFF;

   %% Permit running in non-gui environment
   sCon.sUser.bNoGUI = ON; 

   %% Read PMU data
   sCon.sUser.bReadPMU = OFF;

   %% bCon_ApplyRawDataCorrection should be equal to 1 to force the correction
   %%    of lines with low gain in TSE like sequences
   sCon.sUser.bApplyRawDataCorrection = OFF; 

   %% bCon_ResetFFTScale should be equal to 1 to reset FFTscale and DataCorrection
   %% for each coil to 1
   sCon.sUser.bResetFFTScale = OFF; 

   %% sCon.sUser.bFillToFullFourier should be equal to 1 to fill k-space 
   %% to complete size in the cases of partial Fourier acquisitions
   %% k-space is filled by zeros
   sCon.sUser.bFillToFullFourier = OFF;

   %% Shift k-space DC to center of output matrix
   %%    Shifts the output k-space such that the DC occurs at n/2+1 for all dimensions
   sCon.sUser.bShiftDCToMatrixCenter = OFF;

   %% Additional Flags
   sCon.sUser.bReadPhaseCor = OFF;
   sCon.sUser.bReadRTFeedback = OFF;
   sCon.sUser.bReadNoiseAdj = OFF;
   sCon.sUser.bReadRelSliceNumber = OFF;
   sCon.sUser.bFoundPhaseCor = 0;
   sCon.sUser.bFoundRTFeedback = 0;
   sCon.sUser.bFoundNoiseAdj = 0;
   sCon.sUser.bFoundRelSliceNumber = 0;

   % --------------------------------------------------------------------------
   % Handle inputs
   % --------------------------------------------------------------------------
   sCon.sUser.dIArgOut = 1;

   if nargin > 1

      for dIArg=1:nargin-1

         if strcmp(inputname(dIArg+1),'ApplyRawDataCorrection')
            sCon.sUser.bApplyRawDataCorrection = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'NoGUI')
            sCon.sUser.bNoGUI = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'Silent')
            sCon.sUser.bSilent = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'ResetFFTScale')
            sCon.sUser.bResetFFTScale = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'FillToFullFourier')
            sCon.sUser.bFillToFullFourier = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'ShiftDCToMatrixCenter')
            sCon.sUser.bShiftDCToMatrixCenter = (varargin{dIArg} ~= 0);
         end
         if strcmp(inputname(dIArg+1),'ReadPMU')
            sCon.sUser.bReadPMU = (varargin{dIArg} ~= 0);
            sCon.sUser.iPMUvarargout = sCon.sUser.dIArgOut;
            sCon.sUser.dIArgOut = sCon.sUser.dIArgOut + 1;
         end
         if strcmp(inputname(dIArg+1),'ReadPhaseCor')
            sCon.sUser.bReadPhaseCor = (varargin{dIArg} ~= 0);
            sCon.sUser.iPhaseCorvarargout = sCon.sUser.dIArgOut;
            sCon.sUser.dIArgOut = sCon.sUser.dIArgOut + 1;
         end
         if strcmp(inputname(dIArg+1),'ReadRTFeedback')
            sCon.sUser.bReadRTFeedback = (varargin{dIArg} ~= 0);
            sCon.sUser.iRTFeedbackvarargout = sCon.sUser.dIArgOut;
            sCon.sUser.dIArgOut = sCon.sUser.dIArgOut + 1;
         end
         if strcmp(inputname(dIArg+1),'ReadNoiseAdj')
            sCon.sUser.bReadNoiseAdj = (varargin{dIArg} ~= 0);
            sCon.sUser.iNoiseAdjvarargout = sCon.sUser.dIArgOut;
            sCon.sUser.dIArgOut = sCon.sUser.dIArgOut + 1;
         end
         if strcmp(inputname(dIArg+1),'ReadRelSliceNumber')
            sCon.sUser.bReadRelSliceNumber = (varargin{dIArg} ~= 0);
            sCon.sUser.iRelSliceNumbervarargout = sCon.sUser.dIArgOut;
            sCon.sUser.dIArgOut = sCon.sUser.dIArgOut + 1;
         end

      end

   end

   % GUI to define tFilename for measurement file
   if nargin == 0
      if ~sCon.sUser.bNoGUI
         tFilename = uigetfile('*.dat','Select File to Read');
      else
         if ~sCon.sUser.bSilent
            'Must provide file name to read'
         end
         return
      end
   end
   if exist(tFilename, 'file') ~= 2
      if ~sCon.sUser.bSilent
         if ~sCon.sUser.bNoGUI
            hError = errordlg(['File does not exist: ' tFilename], 'Error');
         else
            ['File does not exist: ', tFilename]
         end
      end
      kSpace = 0;
      for dIArg=1:nargout-1
         varargout{dIArg} = 0;
      end
      return
   end

   if nargin > 0 && nargout < sCon.sUser.dIArgOut
      if ~sCon.sUser.bSilent
         if sCon.sUser.bNoGUI
            'You must have an equal number of kSpace output variables for all the kSpace you expect to read.'
         else
            hError = errordlg(['You must have an equal number of kSpace output variables for all the kSpace you expect to read.'], 'Insufficient Outputs');
         end
      end
      kSpace = 0;
      for dIArg=1:nargout-1
         varargout{dIArg} = 0;
      end
      return;
   end

   text = ['Filename=%s\n', ...
           'ApplyRawDataCorrection=%d\n', ...
           'ResetFFTScale=%d\n', ...
           'FillToFullFourier=%d\n', ...
           'ShiftDCToMatrixCenter=%d\n', ...
           'ReadPMU=%d\n', ...
           'ReadRelSliceNumber=%d\n', ...
           'ReadPhaseCor=%d\n', ...
           'ReadRTFeedback=%d\n', ...
           'ReadNoiseAdj=%d\n'];
   tMsg = sprintf(text, ...
                  tFilename, ...
                  sCon.sUser.bApplyRawDataCorrection, ...
                  sCon.sUser.bResetFFTScale, ...
                  sCon.sUser.bFillToFullFourier, ...
                  sCon.sUser.bShiftDCToMatrixCenter, ...
                  sCon.sUser.bReadPMU, ...
                  sCon.sUser.bReadRelSliceNumber, ...
                  sCon.sUser.bReadPhaseCor, ...
                  sCon.sUser.bReadRTFeedback, ...
                  sCon.sUser.bReadNoiseAdj);
   if ~sCon.sUser.bSilent
      if sCon.sUser.bNoGUI
         tMsg
      else
         uiwait(msgbox(tMsg, 'File Will Be Read With These Parameters','modal'));
      end
   end

   sFileStat = dir(tFilename);
   ullFileSize = uint64(sFileStat.bytes);

   if ~(sCon.sUser.bNoGUI || sCon.sUser.bSilent)
      hProg = waitbar(0, 'Extracting Protocol', 'Name', 'Progress', 'Resize', 'on');
   end

   fid = fopen(tFilename,'r','ieee-le');

   % --------------------------------------------------------------------------
   % First check if VB or VD
   % --------------------------------------------------------------------------
   sMrParcRaidFileHeader.ulSize  = fread(fid,1,'uint32=>uint32');
   sMrParcRaidFileHeader.ulNMeas = fread(fid,1,'uint32=>uint32');
   if not(and(sMrParcRaidFileHeader.ulSize < 10000, sMrParcRaidFileHeader.ulNMeas <= 64))
      if ~sCon.sUser.bSilent
         if ~sCon.sUser.bNoGUI
            hError = errordlg(['Appears to be VB data: ' tFilename], 'Error');
         else
            ['Appears to be VB data: ', tFilename]
         end
      end
      kSpace = 0;
      for dIArg=1:nargout-1
         varargout{dIArg} = 0;
      end
      return
   end

   % --------------------------------------------------------------------------
   % Extract measurement parameters
   % --------------------------------------------------------------------------
   for dIMeas = 1:sMrParcRaidFileHeader.ulNMeas

      sMrParcRaidFileEntry{dIMeas}.ulMeasID = fread(fid,1,'uint32=>uint32');
      sMrParcRaidFileEntry{dIMeas}.ulFileID = fread(fid,1,'uint32=>uint32');
      sMrParcRaidFileEntry{dIMeas}.ullOffset = fread(fid,1,'uint64=>uint64');
      sMrParcRaidFileEntry{dIMeas}.ullSize = fread(fid,1,'uint64=>uint64');
      sMrParcRaidFileEntry{dIMeas}.cPatName = fread(fid,64,'uchar=>char');
      sMrParcRaidFileEntry{dIMeas}.cProtName = fread(fid,64,'uchar=>char');
      sMrParcRaidFileEntry{dIMeas}.sCon = sCon;

   end

   % --------------------------------------------------------------------------
   % Loop through measurements
   % --------------------------------------------------------------------------

   for dIMeas = sMrParcRaidFileHeader.ulNMeas:sMrParcRaidFileHeader.ulNMeas

      % Skip to start of this measurement
      fseek(fid,sMrParcRaidFileEntry{dIMeas}.ullOffset,'bof');

      % Protocol extraction
      sDim.dMeasHeaderSize = fread(fid,1,'int32');
      sDim.ullMeasHeaderSize = uint64(sDim.dMeasHeaderSize);
      cInfo = fread(fid,sDim.dMeasHeaderSize-4,'uchar=>char');
      cInfo = cInfo';
      sProtocol.cLong = cInfo;

      textStart = 'MeasYaps';
      textEnd = 'Phoenix';
      indexStart = strfind(cInfo,textStart) + length(textStart) + 5;
      indexEnd = strfind(cInfo,textEnd) - 3;
      sProtocol.cShort = cInfo(indexStart:indexEnd);
      clear cInfo;
      cInfo = sProtocol.cShort;

      %% Some Acquisition Parameters ASCCONV section
      cSearch = 'sKSpace.ucDimension[^=]+=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      dScanDimension = str2num(t{1}{1});
      if dScanDimension == 4
         sMrParcRaidFileEntry{dIMeas}.sCon.sFile.bFlag3D = ON;
      else
         sMrParcRaidFileEntry{dIMeas}.sCon.sFile.bFlag3D = OFF;
      end

      % PeFTLen is more appropriate here
      cSearch = 'iPEFTLen[^{}]+{\s*([0-9]+)\s*}';
      t = regexp(sProtocol.cLong,cSearch,'tokens');
      sDim.sRecon.dLIN = str2double(t{1}{1});

      cSearch = 'MeasUID[^{}]+{\s*([0-9]+)\s*}';
      t = regexp(sProtocol.cLong,cSearch,'tokens');
      sMrParcRaidFileEntry{dIMeas}.sCon.sFile.tMeasUID = t{1}{1};

      cSearch = 'sKSpace.lPartitions[^=]+=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      sDim.sRecon.dPAR = str2num(t{1}{1});

      cSearch = 'sSliceArray.lSize[^=]+=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      sDim.sRecon.dSLC = str2num(t{1}{1});

      cSearch = '.lRxChannelConnected\s*=\s*';
      t = regexp(cInfo,cSearch,'tokens');
      sDim.sRaw.dCHA = length(t);
      cSearch = 'AdjustSeq%/AdjCoilSensSeq';
      t = strfind(cInfo,cSearch) + length(cSearch);
      if isempty(t) == 0
         sDim.sRaw.dCHA = sDim.sRaw.dCHA-1;
      end

      cSearch = 'lContrasts\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
         sDim.sRaw.dECO = 1;
      else
         sDim.sRaw.dECO = str2num(t{1}{1});
      end

      cSearch = 'lSets\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
         sDim.sRaw.dSET = 1;
      else
         sDim.sRaw.dSET = str2num(t{1}{1});
      end

      cSearch = 'lAverages\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
         sDim.sRaw.dACQ = 1;
      else
         sDim.sRaw.dACQ = str2num(t{1}{1});
      end

      cSearch = 'lRepetitions\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
          sDim.sRaw.dREP = 1;
      else
          sDim.sRaw.dREP = str2num(t{1}{1});
      end

      cSearch = 'sPhysioImaging.lPhases\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
          sDim.sRaw.dPHS = 1;
      else
          sDim.sRaw.dPHS = str2num(t{1}{1});
      end

      cSearch = 'sFastImaging.lEPIFactor\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      sMrParcRaidFileEntry{dIMeas}.sCon.sFile.dEPIFactor = str2num(t{1}{1});

      cSearch = 'sFastImaging.lTurboFactor\s*=\s*([0-9]+)\s*';
      t = regexp(cInfo,cSearch,'tokens');
      if isempty(t) == 1
          sMrParcRaidFileEntry{dIMeas}.sCon.sFile.dTurboFactor = 1;
      else
          sMrParcRaidFileEntry{dIMeas}.sCon.sFile.dTurboFactor = str2num(t{1}{1});
      end

      clear cInfo;
      cInfo = sProtocol.cLong;
      cSearch = 'ParamLong."iMaxNoOfRxChannels';
      t = strfind(cInfo,cSearch) + length(cSearch);
      if isempty(t) == 0
         ss = sscanf(cInfo(t:t+35),'%s');
         sDim.sRaw.dMaxNoOfRxChannels = str2num(ss(isstrprop(ss,'digit')));
         if sDim.sRaw.dMaxNoOfRxChannels < sDim.sRaw.dCHA
            sDim.sRaw.dCHA = sDim.sRaw.dMaxNoOfRxChannels;
         end
      end

      cSearch = 'ParamLong."lNoOfPhaseCorrScans';
      t = strfind(cInfo,cSearch) + length(cSearch);
      if isempty(t) == 0
         ss = sscanf(cInfo(t:t+35),'%s');
         if isempty(str2num(ss(isstrprop(ss,'digit')))) == 0
            sDim.sRaw.dNPhCorScan = str2num(ss(isstrprop(ss,'digit')));
         end
      end
      if sMrParcRaidFileEntry{dIMeas}.sCon.sFile.dTurboFactor > 1
         sDim.sRaw.dNPhCorScan = 1;
      end
      if sMrParcRaidFileEntry{dIMeas}.sCon.sFile.dEPIFactor > 1
         sDim.sRaw.dNPhCorScan = 1;
      end

      %% YAPS Parameters from ASCCONV
      cSearch = 'iNoOfFourierColumns';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.dNoOfFourierColumns = str2num(ss(isstrprop(ss,'digit')));

      cSearch = 'flReadoutOSFactor';
      t = strfind(cInfo,cSearch) + length(cSearch);
      sDim.sRaw.dOSFactorRO = sscanf(cInfo(t:t+10),'%f');
      sDim.sRaw.dOSFactorRO = 2;
      sDim.sRecon.dNoOfNonOSColumns = round(sDim.sRecon.dNoOfFourierColumns/sDim.sRaw.dOSFactorRO);

      cSearch = 'NoOfFourierLines';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.dNoOfFourierLines = str2num(ss(isstrprop(ss,'digit')));

      cSearch = 'iNoOfFourierPartitions';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.dNoOfFourierPartitions = str2num(ss(isstrprop(ss,'digit')));

      cSearch = 'iRoFTLength';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.dRoFTLength = str2num(ss(isstrprop(ss,'digit')));

      cSearch = 'iPEFTLength';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.dPEFTLength = str2num(ss(isstrprop(ss,'digit')));

      cSearch = 'i3DFTLength';
      t = strfind(cInfo,cSearch) + length(cSearch);
      ss = sscanf(cInfo(t:t+35),'%s');
      sDim.sRecon.d3DFTLength = str2num(ss(isstrprop(ss,'digit')));
      if sDim.sRecon.dNoOfFourierPartitions == 1
         sDim.sRecon.d3DFTLength = 1;
      end

      %% Raw Data Correction Factors (NOTE: CURRENTLY OFF!!!!)
      mCoilSelect = containers.Map;
      cSearch = '{\s*{\s*{\s*"[^"]+"[^\n]+';
      coilarray = regexp(cInfo, cSearch, 'match');
      if ~isempty(coilarray)
         coilarray = coilarray{1};
         tElement = '{\s*{\s*"(?<name>[^"]+)"\s*}\s*{\s*(?<fft>[\d\.]+)\s*}\s*{\s*(?<re>[\d\.-]+)\s*}\s*{\s*(?<im>[\d\.-]+)\s*}\s*}';
         mCoilElements = regexp(coilarray, tElement,'names');
         if length(mCoilElements) == sDim.sRaw.dCHA
            for c=1:sDim.sRaw.dCHA
               tName = mCoilElements(c).name;
               sCoilSelect.dFFTScale = sscanf(mCoilElements(c).fft,'%f');
               sCoilSelect.dRawDataCorrectionFactor = complex(sscanf(mCoilElements(c).re,'%f'),sscanf(mCoilElements(c).im,'%f'));
               mCoilSelect(tName) = sCoilSelect;
            end
         end
         if (length(keys(mCoilSelect)) ~= sDim.sRaw.dCHA) && sCon.sUser.bApplyRawDataCorrection
            if ~sCon.sUser.bSilent
               if ~sCon.sUser.bNoGUI
                  hError = errordlg('Non-unique channel names in CoilSelect', 'Error');
               else
                  'Non-unique channel names in CoilSelect'
               end
            end
            kSpace = 0;
            for dIArg=1:nargout-1
               varargout{dIArg} = 0;
            end
            return
         end
      end

      if sCon.sUser.bReadRelSliceNumber

         tSearch = 'ParamLong."relSliceNumber">\s*[^n]*{';
         dIStart = regexp(cInfo,tSearch);
         if ~isempty(dIStart)
            dIStart = dIStart(1) + strfind(cInfo(dIStart(1):dIStart(1) + length(tSearch)+20),'{');;
            dIStart = dIStart(1);
            dIEnd = strfind(cInfo(dIStart(1):end),'}');
            if ~isempty(dIEnd)
               dIEnd = dIStart + dIEnd(1) - 2;
               tSearch = '([^0-9\s]*[\-0-9]+)(\s+)';
               cRelSliceNum = regexp(cInfo(dIStart:dIEnd),tSearch,'tokens');
               if ~isempty(cRelSliceNum)
                  dRelativeSliceNumber = zeros(1,length(cRelSliceNum));
                  for i=1:length(cRelSliceNum)
                     dRelativeSliceNumber(i) = str2double(cRelSliceNum{i}{1});
                  end
                  RelativeSliceNumber{dIMeas} = dRelativeSliceNumber;
                  sCon.sUser.bFoundRelSliceNumber = 1;
               end
            end
         end
      end

      clear cInfo sCoilSelect;
      cInfo = sProtocol.cShort;

      %% FFT Correction Factors
      for c=1:sDim.sRaw.dCHA
         cSearch = 'asList\[';
         cSearch = sprintf('%s%d',cSearch,c-1);
         cSearch = [cSearch,'\]\.sCoilElementID.tElement\s*=\s*"([^"]+)"\s*'];
         t = regexp(cInfo,cSearch,'tokens');
         sCoilSelectMeas.tElement = t{1}{1};
         cSearch = 'asList\[';
         cSearch = sprintf('%s%d',cSearch,c-1);
         cSearch = [cSearch,'\]\.lRxChannelConnected\s*=\s*([0-9.]+)\s*'];
         t = regexp(cInfo,cSearch,'tokens');
         sCoilSelectMeas.lRxChannel = str2num(t{1}{1});
         cSearch = 'asList\[';
         cSearch = sprintf('%s%d',cSearch,c-1);
         cSearch = [cSearch,'\]\.lADCChannelConnected\s*=\s*([0-9.]+)\s*'];
         t = regexp(cInfo,cSearch,'tokens');
         lADCChannel = str2num(t{1}{1});
         if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bResetFFTScale == 1 
            sCoilSelectMeas.flFFTCorrectionFactor = 1.0; 
         else
            cSearch = 'aFFT_SCALE\[';
            cSearch = sprintf('%s%d',cSearch,c-1);
            cSearch = [cSearch,'\]\.flFactor\s*=\s*([0-9.]+)\s*'];
            t = regexp(cInfo,cSearch,'tokens');
            sCoilSelectMeas.flFFTCorrectionFactor = str2double(t{1}{1});
         end
         keySet(c) = lADCChannel;
         valueSet(c) = sCoilSelectMeas;
      end
      mCoilSelectMeas = containers.Map(keySet,arrayfun(@(x) ({x}), valueSet));

      clear keySet valueSet cInfo;

      if sDim.sRecon.dNoOfFourierLines == sDim.sRecon.dLIN && sDim.sRecon.dNoOfFourierPartitions == sDim.sRecon.dPAR
         sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bFillToFullFourier = OFF;
      end

      if sMrParcRaidFileEntry{dIMeas}.sCon.sFile.bFlag3D == OFF
         sDim.sRecon.dPAR = sDim.sRecon.dNoOfFourierPartitions;
      end

      clear sProtocol;

      if sDim.sRecon.dNoOfFourierLines > sDim.sRecon.dLIN
          sDim.sRecon.dLIN = sDim.sRecon.dNoOfFourierLines;
      end

      if sDim.sRecon.dNoOfFourierPartitions > sDim.sRecon.dPAR
          sDim.sRecon.dPAR = sDim.sRecon.dNoOfFourierPartitions;
      end


      % --------------------------------------------------------------------------
      % First pass to chunk this measurement into blocks with the same scan lengths
      % --------------------------------------------------------------------------
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.5, hProg, 'Measuring number of scans');
      end

      % With readout and PMU being different sizes, 
      %    we now have to chunk in blocks of same sized
      %    lines.  That requires scanning all line lengths
      %    up front.
      ullSizeAllScansInBytes = sMrParcRaidFileEntry{dIMeas}.ullSize - sDim.ullMeasHeaderSize;
      st = fseek(fid,sMrParcRaidFileEntry{dIMeas}.ullOffset + sDim.ullMeasHeaderSize,'bof');
      % Read the first line to get the DMA length and line type
      cMDHScan = fread(fid,ullMDHScanSize,'uchar=>uchar');
      ullDMALength = cMDHScan(1:4);
      ullDMALength(4) = 0;
      ulEvalMask1 = typecast(cMDHScan(41:44),'uint32');
      bHasACQEND = (bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND);
      ullDMALength = uint64(typecast(ullDMALength,'uint32'));
      ullIRelativePosition = uint64(1);
      dIBlock = 1;
      sMDHBlock(dIBlock).ullOffset = sMrParcRaidFileEntry{dIMeas}.ullOffset + sDim.ullMeasHeaderSize;
      sMDHBlock(dIBlock).ullLength = ullDMALength;
      sMDHBlock(dIBlock).ullNScan = uint64(0);
      sMDHBlock(dIBlock).bIsSYNCDATA = (bitand(ulEvalMask1,ulMDH_SYNCDATA) == ulMDH_SYNCDATA);
      sMDHBlock(dIBlock).bIsACQEND = (bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND);
      % NOTE: We make a very simple assumption that the same number of
      %    samples per scan will hold for scans of the same ullDMALength
      sMDHBlock(dIBlock).ullNSamplesInScan = uint64(typecast(cMDHScan(49:50),'uint16'));
      % Similarly for the number of channels (which mean nothing in SYNCDATA scans)
      sMDHBlock(dIBlock).usNChannelUsed = typecast(cMDHScan(51:52),'uint16');
      ullNScanInMeas = uint64(0);
      while ((ullIRelativePosition < ullSizeAllScansInBytes) && ~bHasACQEND)
         if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
            cWinMsg = 'Pass 1 of 3: Parsing blocks';
            waitbar(double(ullIRelativePosition)/double(ullSizeAllScansInBytes), hProg, cWinMsg);
         end
         st = fseek(fid,sMDHBlock(dIBlock).ullOffset + sMDHBlock(dIBlock).ullNScan * ullDMALength,'bof');
         ullBytesRemaining = ullSizeAllScansInBytes-ullIRelativePosition+uint64(1);
         if ullDMALength < 2048
            ullNScanToRead=uint64(10);  % We don't usually get a lot short lines
         else
            ullNScanToRead = max(idivide(ullMaxBytesPerRead, ullDMALength,'fix'),uint64(1));
         end   
         if (ullNScanToRead*ullDMALength > ullBytesRemaining)
            ullNScanToRead = idivide(ullBytesRemaining,ullDMALength,'fix');
            if ullNScanToRead < 1
               cMDHScan = fread(fid,ullMDHScanSize,'uchar=>uchar');
               ullDMALengthNext = cMDHScan(1:4);
               ullDMALengthNext(4) = 0;
               ullDMALengthNext = uint64(typecast(ullDMALengthNext,'uint32'));
               ullNScanToRead = uint64(1);
               if ullDMALengthNext == ullDMALength
                  error_message
               else
                  ulEvalMask1 = typecast(cMDHScan(41:44),'uint32');
                  dIBlock = dIBlock + 1;
                  sMDHBlock(dIBlock).ullOffset = sMDHBlock(dIBlock-1).ullOffset + sMDHBlock(dIBlock-1).ullNScan * ullDMALength;
                  sMDHBlock(dIBlock).ullNScan = uint64(0);
                  ullDMALength = ullDMALengthNext;
                  sMDHBlock(dIBlock).ullLength = ullDMALength;
                  sMDHBlock(dIBlock).bIsSYNCDATA = (bitand(ulEvalMask1,ulMDH_SYNCDATA) == ulMDH_SYNCDATA);
                  sMDHBlock(dIBlock).bIsACQEND = (bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND);
                  sMDHBlock(dIBlock).ullNSamplesInScan = uint64(typecast(cMDHScan(49:50),'uint16'));
                  sMDHBlock(dIBlock).usNChannelUsed = typecast(cMDHScan(51:52),'uint16');
                  st = fseek(fid,sMDHBlock(dIBlock).ullOffset + sMDHBlock(dIBlock).ullNScan * ullDMALength,'bof');
               end
            end
         end
         cData = fread(fid,ullDMALength*ullNScanToRead,'uchar=>uchar');
         for iStart=1:ullDMALength:ullDMALength*ullNScanToRead
            ullDMALengthNext = cData(iStart:iStart+3);
            ullDMALengthNext(4) = 0;
            ullDMALengthNext = uint64(typecast(ullDMALengthNext,'uint32'));
            ulEvalMask1 = typecast(cData(iStart+40:iStart+43),'uint32');
            bHasACQEND = bHasACQEND || sMDHBlock(dIBlock).bIsACQEND;
            if ullDMALengthNext == ullDMALength
               sMDHBlock(dIBlock).ullNScan = sMDHBlock(dIBlock).ullNScan + 1;
               ullIRelativePosition = ullIRelativePosition + ullDMALength;
            else
               dIBlock = dIBlock + 1;
               sMDHBlock(dIBlock).ullOffset = sMDHBlock(dIBlock-1).ullOffset + sMDHBlock(dIBlock-1).ullNScan * ullDMALength;
               sMDHBlock(dIBlock).ullNScan = uint64(0);
               sMDHBlock(dIBlock).bIsSYNCDATA = (bitand(ulEvalMask1,ulMDH_SYNCDATA) == ulMDH_SYNCDATA);
               sMDHBlock(dIBlock).bIsACQEND = (bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND);
               sMDHBlock(dIBlock).ullNSamplesInScan = uint64(typecast(cData(iStart+48:iStart+49),'uint16'));
               sMDHBlock(dIBlock).usNChannelUsed = typecast(cData(iStart+50:iStart+51),'uint16');
               ullDMALength = ullDMALengthNext;
               sMDHBlock(dIBlock).ullLength = ullDMALength;
               break;
            end
            if bHasACQEND
               break;
            end
         end % for over block of scans with same dmalength
      end % while over full file
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         cWinMsg = 'Pass 1 of 3: Parsing blocks';
         waitbar(double(ullIRelativePosition)/double(ullSizeAllScansInBytes), hProg, cWinMsg);
      end
      if ((ullIRelativePosition < ullSizeAllScansInBytes) && bHasACQEND)         
         if ~sCon.sUser.bSilent
            if ~sCon.sUser.bNoGUI
               hWarn = warndlg('ACQ_END before end of file.');
            else
               'ACQ_END before end of file.'
            end
         end
      end
      if ~bHasACQEND
         if ~sCon.sUser.bSilent
            if ~sCon.sUser.bNoGUI
               hWarn = warndlg('No ACQ_END encountered.  Data possibly incomplete or corrupted.');
            else
               'No ACQ_END encountered.  Data possibly incomplete or corrupted.'
            end
         end
      end
      dNBlockSYNCDATA = 0;
      dNBlockACQEND = 0;
      dNBlockScan=0;
      for dIBlock=1:length(sMDHBlock)
         if sMDHBlock(dIBlock).bIsACQEND
            dNBlockACQEND = dNBlockACQEND + 1;
            ullBlockNScanACQEND(dNBlockACQEND) = sMDHBlock(dIBlock).ullNScan;
         elseif sMDHBlock(dIBlock).bIsSYNCDATA
            dNBlockSYNCDATA = dNBlockSYNCDATA + 1;
            ullBlockNScanSYNCDATA(dNBlockSYNCDATA) = sMDHBlock(dIBlock).ullNScan;
         else
            dNBlockScan = dNBlockScan + 1;
            ullBlockNScan(dNBlockScan) = sMDHBlock(dIBlock).ullNScan;
            ullDMALengthScan(dNBlockScan) = sMDHBlock(dIBlock).ullLength;
            ullNSamplesInScan(dNBlockScan) = sMDHBlock(dIBlock).ullNSamplesInScan;
         end
      end
      if dNBlockSYNCDATA > 0
         ullNScanSYNCDATA = sum(ullBlockNScanSYNCDATA(:));
      end
      if dNBlockACQEND > 0
         ullNScanACQEND = sum(ullBlockNScanACQEND(:));
      end
      if dNBlockScan < 1
         if ~sCon.sUser.bSilent
            if ~sCon.sUser.bNoGUI
               hError = errordlg('No readouts detected.', 'Error');
            else
               'No readouts detected.'
            end
         end
         kSpace = 0;
         for dIArg=1:nargout-1
            varargout{dIArg} = 0;
         end
         return
      end
      ullNScanInMeas = uint64(sum(ullBlockNScan(:)));
      if (length(unique(ullNSamplesInScan)) ~= 1) || (length(unique(ullDMALengthScan)) ~= 1)
         if ~sCon.sUser.bSilent
            if ~sCon.sUser.bNoGUI
               hError = errordlg('Multiple different samples per scan detected.', 'Error');
            else
               'Multiple different samples per scan detected.'
            end
         end
         kSpace = 0;
         for dIArg=1:nargout-1
            varargout{dIArg} = 0;
         end
         return
      end
      ullNSamplesInScan = unique(ullNSamplesInScan);
      ullNSamplesInScan = ullNSamplesInScan(1);
      % Read the first real scan channel header to get scan length in bytes
      for dIBlock=1:length(sMDHBlock)
         if sMDHBlock(dIBlock).bIsACQEND || sMDHBlock(dIBlock).bIsSYNCDATA
            continue
         else
            st = fseek(fid,sMDHBlock(dIBlock).ullOffset,'bof');
            cMDHScan = fread(fid,ullMDHScanSize,'uchar=>uchar');
            cMDHChan = fread(fid,ullMDHChanSize,'uchar=>uchar');
            ullSizeChanInBytes = uint64(bitshift(typecast(cMDHChan(1:4),'uint32'),-8));
            break
         end
      end

      % --------------------------------------------------------------------------
      % Second pass, Part 1 to extract PMU
      % --------------------------------------------------------------------------
      if (dNBlockSYNCDATA > 0) && sCon.sUser.bReadPMU
         sPMU.sTimestamp.dNColumn = 0;
         sPMU.sTimestamp.dRaw = reshape([],2,0);
         sPMU.sTimestamp.dData = reshape([],2,0);
         sWaveform.tName = '';
         sWaveform.dNColumn = 0;
         sWaveform.dNRow = 0;
         sWaveform.dPeriod = 0;
         sWaveform.ulData = [];
         sWaveform.ulTrigger = [];
         cNames = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
         sPMU.sWaveforms = repmat(sWaveform,1,length(cNames));
         for dIType=1:length(cNames)
            sPMU.sWaveforms(dIType).tName = cNames(dIType);
         end
         for dIBlock = 1:length(sMDHBlock)

            % Skip SYNCDATA
            if ~sMDHBlock(dIBlock).bIsSYNCDATA 
               continue
            end

            st = fseek(fid,sMDHBlock(dIBlock).ullOffset,'bof');
            ullNScanToRead = sMDHBlock(dIBlock).ullNScan;
            cData = fread(fid,[sMDHBlock(dIBlock).ullLength ullNScanToRead],'uchar=>uchar');
        
            offset = ullMDHScanSize + ullMDHPMUSize;

            % PMU Header
            ulTimestamp = cData(offset+1:offset+8,:);
            ulTimestamp = reshape(typecast(ulTimestamp(:),'uint32'),2,ullNScanToRead);
            ulCounter = cData(offset+9:offset+12,:);
            ulCounter = typecast(ulCounter(:),'uint32');
            ulDuration = cData(offset+13:offset+16,:);
            ulDuration = typecast(ulDuration(:),'uint32');

            sPMU.sTimestamp.dRaw(:,sPMU.sTimestamp.dNColumn+1:sPMU.sTimestamp.dNColumn + ullNScanToRead) = ulTimestamp;
            sPMU.sTimestamp.dData(:,sPMU.sTimestamp.dNColumn+1:sPMU.sTimestamp.dNColumn + ullNScanToRead) = double(ulTimestamp) * 2.5d-3;
            sPMU.sTimestamp.dNColumn = sPMU.sTimestamp.dNColumn + ullNScanToRead;

            ulPMUData = cData(offset+17:end,:);
            ullNPMUPoints = idivide(sMDHBlock(dIBlock).ullLength - offset - uint64(16),uint64(4),'fix');
            ulPMUData = reshape(typecast(ulPMUData(:),'uint32'),ullNPMUPoints,ullNScanToRead);

            dPMUPoint = 1;
            while dPMUPoint <= ullNPMUPoints

               ulType = bitand(1.0*ulPMUData(dPMUPoint,:),hex2dec('01FF0000'),'uint32');
               if length(unique(ulType)) ~= 1
                  if ~sCon.sUser.bSilent
                     if sCon.sUser.bNoGUI
                        'Unexpected PMU type data'
                     else
                        hError = errordlg('Unexpected PMU type data', 'PMU Problem');
                     end
                  end
                  kSpace = 0;
                  for dIArg=1:nargout-1
                     varargout{dIArg} = 0;
                  end
                  return;
               else
                  ulType = unique(ulType);
               end
               % Check for end of data magic
               if ulType == hex2dec('01FF0000')
                  break 
               end
               ulType = bitand(1.0*ulType,hex2dec('000F0000'));
               ulType = bitshift(ulType,-16);
               ulPeriod = ulPMUData(dPMUPoint+1,:);
               ulNPeriod = idivide(ulDuration(:),ulPeriod(:),'fix');
               if (length(unique(ulNPeriod)) ~= 1) || (length(unique(ulPeriod(:))) ~= 1)
                  if ~sCon.sUser.bSilent
                     if sCon.sUser.bNoGUI
                        'Unexpected PMU period or point data'
                     else
                        hError = errordlg('Unexpected PMU period or point data', 'PMU Problem');
                     end
                  end
                  kSpace = 0;
                  for dIArg=1:nargout-1
                     varargout{dIArg} = 0;
                  end
                  return;
               else
                  ulNPeriod = unique(ulNPeriod);
                  ulPeriod = unique(ulPeriod);
               end
               
               ulPMUPoint = uint32(dPMUPoint);
               ulPMUThisData = ulPMUData(ulPMUPoint+2:ulPMUPoint+2+ulNPeriod-1,:);
               % Waveforms are the lower 12 bits
               ulWaveform = bitand(1.0*ulPMUThisData,hex2dec('00000fff'),'uint32');
               % Triggers in higher bits
               ulTriggers = bitand(1.0*ulPMUThisData,hex2dec('fffff000'),'uint32');

               if sPMU.sWaveforms(ulType).dPeriod == 0
                  sPMU.sWaveforms(ulType).dPeriod = (double(ulPeriod) * 100d-6);
               else
                  if sPMU.sWaveforms(ulType).dPeriod ~= (double(ulPeriod)*100d-6)
                     if ~sCon.sUser.bSilent
                        if sCon.sUser.bNoGUI
                           'Unexpected PMU period'
                        else
                           hError = errordlg('Unexpected PMU period', 'PMU Problem');
                        end
                     end
                     kSpace = 0;
                     for dIArg=1:nargout-1
                        varargout{dIArg} = 0;
                     end
                     return;

                  end
               end
               sPMU.sWaveforms(ulType).ulData(1:ulNPeriod,sPMU.sWaveforms(ulType).dNColumn+1:sPMU.sWaveforms(ulType).dNColumn + ullNScanToRead) = ulWaveform;
               sPMU.sWaveforms(ulType).ulTrigger(1:ulNPeriod,sPMU.sWaveforms(ulType).dNColumn+1:sPMU.sWaveforms(ulType).dNColumn + ullNScanToRead) = ulTriggers;
               sPMU.sWaveforms(ulType).dNRow(sPMU.sWaveforms(ulType).dNColumn+1:sPMU.sWaveforms(ulType).dNColumn + ullNScanToRead) = ulNPeriod;
               sPMU.sWaveforms(ulType).dNColumn = sPMU.sWaveforms(ulType).dNColumn + ullNScanToRead;

               % Need to jump our pointer forward
               dPMUPoint = dPMUPoint + 2 + double(ulNPeriod);

            end

         end % for over blocks

         % Consolidate values
         sWaveformOut.tName = '';
         sWaveformOut.dNPoint = 0;
         sWaveformOut.dTime = [];
         sWaveformOut.ulData = [];
         sWaveformOut.ulTrigger = [];
         cNames = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
         sPMUOut.sTimestamp = sPMU.sTimestamp;
         sPMUOut.sWaveforms = repmat(sWaveformOut,1,length(cNames));
         for dIType=1:length(cNames)
            sPMUOut.sWaveforms(dIType).tName = cNames(dIType);
            dNRowTotal = sum(sPMU.sWaveforms(dIType).dNRow(:));
            sPMUOut.sWaveforms(dIType).ulData = zeros(1,dNRowTotal,'uint32');
            sPMUOut.sWaveforms(dIType).ulTrigger = zeros(1,dNRowTotal,'uint32');
            sPMUOut.sWaveforms(dIType).dTime = zeros(1,dNRowTotal);
            iStart = 1;
            sPMUOut.sWaveforms(dIType).dNPoint = dNRowTotal;
            for iCol=1:sPMU.sWaveforms(dIType).dNColumn

               nRow = sPMU.sWaveforms(dIType).dNRow(iCol);
               timeRelative = (0:(nRow-1)) * sPMU.sWaveforms(dIType).dPeriod;
               timeStart = sPMU.sTimestamp.dData(2,iCol);
               timeValues = timeStart + timeRelative;
               waveValues = sPMU.sWaveforms(dIType).ulData(1:nRow,iCol);
               triggers = sPMU.sWaveforms(dIType).ulTrigger(1:nRow,iCol);

               sPMUOut.sWaveforms(dIType).dTime(iStart:iStart+nRow-1) = timeValues(:);
               sPMUOut.sWaveforms(dIType).ulData(iStart:iStart+nRow-1) = waveValues(:);
               sPMUOut.sWaveforms(dIType).ulTrigger(iStart:iStart+nRow-1) = triggers(:);

               iStart = iStart + nRow;

            end
         end

         clear sPMU

         PMUOut{dIMeas} = sPMUOut;

      end % if we're loading PMU

      % --------------------------------------------------------------------------
      % Second pass, Part 2 to assemble dimensions and read order for readouts
      % --------------------------------------------------------------------------
      ullNScanRemainInMeas = ullNScanInMeas;
      ullIScanRead = uint64(0);
      cDimName          =        {'cha', 'lin', 'slc', 'par', 'acq', 'eco', 'phs', 'rep', 'set'};
      usDimPHASCOR       = uint16([    0,     0,     0,     0,     0,     0,     0,     0,     0]);
      dRankPHASCOR       =        [    0,     0,     0,     0,     0,     0,     0,     0,     0];
      usDimRTFEEDBACK    = uint16([    0,     0,     0,     0,     0,     0,     0,     0,     0]);
      dRankRTFEEDBACK    =        [    0,     0,     0,     0,     0,     0,     0,     0,     0];
      usDimNOISEADJSCAN  = uint16([    0,     0,     0,     0,     0,     0,     0,     0,     0]);
      dRankNOISEADJSCAN  =        [    0,     0,     0,     0,     0,     0,     0,     0,     0];
      usDimNormal        = uint16([    0,     0,     0,     0,     0,     0,     0,     0,     0]);
      dRankNormal        =        [    0,     0,     0,     0,     0,     0,     0,     0,     0];
      dNScanPHASCOR = 0;
      dNScanNOISEADJSCAN = 0;
      dNScanRTFEEDBACK = 0;
      dNScanNormal = 0;
      dNSegPHASCOR = 0;
      dNSegRTFEEDBACK = 0;
      dNSegNOISEADJSCAN = 0;
      dNSegNormal = 0;
      time0 = tic;

      for dIBlock = 1:length(sMDHBlock)

         % Skip SYNCDATA
         if sMDHBlock(dIBlock).bIsSYNCDATA || sMDHBlock(dIBlock).bIsACQEND
            continue
         end

         ullMaxScanPerRead = max(idivide(ullMaxBytesPerRead,sMDHBlock(dIBlock).ullLength,'fix'),uint64(1));
         ullNScanToRead = min(ullMaxScanPerRead,sMDHBlock(dIBlock).ullNScan);

         st = fseek(fid,sMDHBlock(dIBlock).ullOffset,'bof');
         ullNScanRemainInBlock = sMDHBlock(dIBlock).ullNScan;

         while ullNScanRemainInBlock > 0

            dTimeElapsed = toc(time0);
            rate = dTimeElapsed / (double(ullIScanRead) + 1);
            timeRemain = double(2 * ullNScanInMeas - ullIScanRead) * rate;
            if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
               cWinMsg = sprintf('Pass 2 of 3\nScanning dimensions from MDH\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(ullIScanRead)/double(ullNScanInMeas), hProg, cWinMsg);
            end
            cData = fread(fid,[sMDHBlock(dIBlock).ullLength ullNScanToRead],'uchar=>uchar');
            ulEvalMask1 = cData(41:44,:);
            ulEvalMask1 = typecast(ulEvalMask1(:),'uint32');
            %usCutOffDataPre = cData(81:82,:);
            %usCutOffDataPre = typecast(usCutOffDataPre(:),'uint16');
            %usCutOffDataPost = cData(83:84,:);
            %usCutOffDataPost = typecast(usCutOffDataPost(:),'uint16');
            %bIsACQEND = bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND;
            %bIsREFLECT = bitand(ulEvalMask1,ulMDH_REFLECT) == ulMDH_REFLECT;
            bIsNOISEADJSCAN = bitand(ulEvalMask1,ulMDH_NOISEADJSCAN) == ulMDH_NOISEADJSCAN;
            bIsPHASCOR = bitand(ulEvalMask1,ulMDH_PHASCOR) == ulMDH_PHASCOR;
            bIsRTFEEDBACK = bitand(ulEvalMask1,ulMDH_RTFEEDBACK) == ulMDH_RTFEEDBACK;
            %bIsRAWDATACORRECTION = bitand(ulEvalMask1,ulMDH_RAWDATACORRECTION) == ulMDH_RAWDATACORRECTION;
            usIndexLIN = cData(53:54,:);
            usIndexLIN = typecast(usIndexLIN(:),'uint16');
            usIndexACQ = cData(55:56,:);
            usIndexACQ = typecast(usIndexACQ(:),'uint16');
            usIndexSLC = cData(57:58,:);
            usIndexSLC = typecast(usIndexSLC(:),'uint16');
            usIndexPAR = cData(59:60,:);
            usIndexPAR = typecast(usIndexPAR(:),'uint16');
            usIndexECO = cData(61:62,:);
            usIndexECO = typecast(usIndexECO(:),'uint16');
            usIndexPHS = cData(63:64,:);
            usIndexPHS = typecast(usIndexPHS(:),'uint16');
            usIndexREP = cData(65:66,:);
            usIndexREP = typecast(usIndexREP(:),'uint16');
            usIndexSET = cData(67:68,:);
            usIndexSET = typecast(usIndexSET(:),'uint16');
            usIndexSEG = cData(69:70,:);
            usIndexSEG = typecast(usIndexSEG(:),'uint16');
            usChannelID = zeros(2,sMDHBlock(dIBlock).usNChannelUsed,ullNScanToRead,'uint8');
            for dICHA = 1:sMDHBlock(dIBlock).usNChannelUsed
               ullOffset  = uint64(ullMDHScanSize) + uint64(dICHA-1) * ullSizeChanInBytes;
               usChannelID(:,dICHA,:) = reshape(cData(ullOffset+25:ullOffset+26,:),[2,1,ullNScanToRead]);
            end
            usChannelID = reshape(typecast(usChannelID(:),'uint16'),sMDHBlock(dIBlock).usNChannelUsed,ullNScanToRead);
            usChannelIDUniq = unique(usChannelID(:));
            usChannelIDUniqNum = length(usChannelIDUniq);
            if (usChannelIDUniqNum ~= sMDHBlock(dIBlock).usNChannelUsed) || (sum(usChannelIDUniq == usChannelID(:,1)) ~= sMDHBlock(dIBlock).usNChannelUsed)
               if ~sCon.sUser.bSilent
                  if sCon.sUser.bNoGUI
                     ['Mismatch between nChannel and coils found: ' int2str(usChannelIDUniqNum) ' ' int2str(sMDHBlock(dIBlock).usNChannelUsed)]
                  else
                     hError = errordlg(['Mismatch between nChannel and coils found: ' int2str(usChannelIDUniqNum) ' ' int2str(sMDHBlock(dIBlock).usNChannelUsed)], 'Dimension Conflict');
                  end
               end
               kSpace = 0;
               for dIArg=1:nargout-1
                  varargout{dIArg} = 0;
               end
               return;
            end
            dIndexNOISEADJSCAN = find(bIsNOISEADJSCAN);
            if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadNoiseAdj == ON && numel(dIndexNOISEADJSCAN) > 0 
               dNScanNOISEADJSCAN = dNScanNOISEADJSCAN + numel(dIndexNOISEADJSCAN);
               dNSegNOISEADJSCAN = max(dNSegNOISEADJSCAN, max(usIndexSEG(dIndexNOISEADJSCAN(:))));
               usDimNOISEADJSCAN = max(usDimNOISEADJSCAN, ...
                                     [sMDHBlock(dIBlock).usNChannelUsed-1, ...
                                      max(usIndexLIN(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexSLC(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexPAR(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexACQ(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexECO(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexPHS(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexREP(dIndexNOISEADJSCAN(:))), ...
                                      max(usIndexSET(dIndexNOISEADJSCAN(:)))]);
               if ullNScanToRead > 1
                  bChange = [dIndexNOISEADJSCAN(2:end) > 0 ...
                            usIndexLIN(dIndexNOISEADJSCAN(2:end)) ~= usIndexLIN(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexSLC(dIndexNOISEADJSCAN(2:end)) ~= usIndexSLC(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexPAR(dIndexNOISEADJSCAN(2:end)) ~= usIndexPAR(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexACQ(dIndexNOISEADJSCAN(2:end)) ~= usIndexACQ(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexECO(dIndexNOISEADJSCAN(2:end)) ~= usIndexECO(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexPHS(dIndexNOISEADJSCAN(2:end)) ~= usIndexPHS(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexREP(dIndexNOISEADJSCAN(2:end)) ~= usIndexREP(dIndexNOISEADJSCAN(1:end-1)) ...
                            usIndexSET(dIndexNOISEADJSCAN(2:end)) ~= usIndexSET(dIndexNOISEADJSCAN(1:end-1))];
                  dRankNOISEADJSCAN = dRankNOISEADJSCAN + sum(bChange, 1);
               end
               if sMDHBlock(dIBlock).usNChannelUsed < 2
                  dRankNOISEADJSCAN(1)=0;
               end
            end
            dIndexPHASCOR = find(bIsPHASCOR);
            if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadPhaseCor == ON && sDim.sRaw.dNPhCorScan > 0 && numel(dIndexPHASCOR) > 0
               dNScanPHASCOR = dNScanPHASCOR + numel(dIndexPHASCOR);
               dNSegPHASCOR = max(dNSegPHASCOR, max(usIndexSEG(dIndexPHASCOR(:))));
               usDimPHASCOR = max(usDimPHASCOR, ...
                                [sMDHBlock(dIBlock).usNChannelUsed-1, ...
                                 max(usIndexLIN(dIndexPHASCOR(:))), ...
                                 max(usIndexSLC(dIndexPHASCOR(:))), ...
                                 max(usIndexPAR(dIndexPHASCOR(:))), ...
                                 max(usIndexACQ(dIndexPHASCOR(:))), ...
                                 max(usIndexECO(dIndexPHASCOR(:))), ...
                                 max(usIndexPHS(dIndexPHASCOR(:))), ...
                                 max(usIndexREP(dIndexPHASCOR(:))), ...
                                 max(usIndexSET(dIndexPHASCOR(:)))]);
               if ullNScanToRead > 1
                  bChange = [dIndexPHASCOR(2:end) > 0 ...
                            usIndexLIN(dIndexPHASCOR(2:end)) ~= usIndexLIN(dIndexPHASCOR(1:end-1)) ...
                            usIndexSLC(dIndexPHASCOR(2:end)) ~= usIndexSLC(dIndexPHASCOR(1:end-1)) ...
                            usIndexPAR(dIndexPHASCOR(2:end)) ~= usIndexPAR(dIndexPHASCOR(1:end-1)) ...
                            usIndexACQ(dIndexPHASCOR(2:end)) ~= usIndexACQ(dIndexPHASCOR(1:end-1)) ...
                            usIndexECO(dIndexPHASCOR(2:end)) ~= usIndexECO(dIndexPHASCOR(1:end-1)) ...
                            usIndexPHS(dIndexPHASCOR(2:end)) ~= usIndexPHS(dIndexPHASCOR(1:end-1)) ...
                            usIndexREP(dIndexPHASCOR(2:end)) ~= usIndexREP(dIndexPHASCOR(1:end-1)) ...
                            usIndexSET(dIndexPHASCOR(2:end)) ~= usIndexSET(dIndexPHASCOR(1:end-1))];
                  dRankPHASCOR = dRankPHASCOR + sum(bChange, 1);
               end
               if sMDHBlock(dIBlock).usNChannelUsed < 2
                  dRankPHASCOR(1)=0;
               end
            end
            dIndexRTFEEDBACK = find(bIsRTFEEDBACK);
            if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadRTFeedback == ON && numel(dIndexRTFEEDBACK) > 0
               dNScanRTFEEDBACK = dNScanRTFEEDBACK + numel(dIndexRTFEEDBACK);
               dNSegRTFEEDBACK = max(dNSegRTFEEDBACK, max(usIndexSEG(dIndexRTFEEDBACK(:))));
               usDimRTFEEDBACK = max(usDimRTFEEDBACK, ...
                                   [sMDHBlock(dIBlock).usNChannelUsed-1, ...
                                    max(usIndexLIN(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexSLC(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexPAR(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexACQ(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexECO(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexPHS(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexREP(dIndexRTFEEDBACK(:))), ...
                                    max(usIndexSET(dIndexRTFEEDBACK(:)))]);
               if ullNScanToRead > 1
                  bChange = [dIndexRTFEEDBACK(2:end) > 0 ...
                            usIndexLIN(dIndexRTFEEDBACK(2:end)) ~= usIndexLIN(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexSLC(dIndexRTFEEDBACK(2:end)) ~= usIndexSLC(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexPAR(dIndexRTFEEDBACK(2:end)) ~= usIndexPAR(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexACQ(dIndexRTFEEDBACK(2:end)) ~= usIndexACQ(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexECO(dIndexRTFEEDBACK(2:end)) ~= usIndexECO(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexPHS(dIndexRTFEEDBACK(2:end)) ~= usIndexPHS(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexREP(dIndexRTFEEDBACK(2:end)) ~= usIndexREP(dIndexRTFEEDBACK(1:end-1)) ...
                            usIndexSET(dIndexRTFEEDBACK(2:end)) ~= usIndexSET(dIndexRTFEEDBACK(1:end-1))];
                  dRankRTFEEDBACK = dRankRTFEEDBACK + sum(bChange, 1);
               end
               if sMDHBlock(dIBlock).usNChannelUsed < 2
                  dRankRTFEEDBACK(1)=0;
               end
            end
            dIndexNormal = find(bIsNOISEADJSCAN == 0 & bIsRTFEEDBACK == 0 & bIsPHASCOR == 0);
            if numel(dIndexNormal) > 0 
               if dNScanNormal == 0
                  %KSpaceCentreColumn = typecast(cData(85:86,dIndexNormal(1)),'uint16');
                  KSpaceCentreLineNo = typecast(cData(97:98,dIndexNormal(1)),'uint16');
                  KSpaceCentrePartitionNo = typecast(cData(99:100,dIndexNormal(1)),'uint16');
               end
               dNScanNormal = dNScanNormal + numel(dIndexNormal);
               dNSegNormal = max(dNSegNormal, max(usIndexSEG(dIndexNormal(:))));
               usDimNormal = max(usDimNormal, ...
                               [sMDHBlock(dIBlock).usNChannelUsed-1, ...
                                max(usIndexLIN(dIndexNormal(:))), ...
                                max(usIndexSLC(dIndexNormal(:))), ...
                                max(usIndexPAR(dIndexNormal(:))), ...
                                max(usIndexACQ(dIndexNormal(:))), ...
                                max(usIndexECO(dIndexNormal(:))), ...
                                max(usIndexPHS(dIndexNormal(:))), ...
                                max(usIndexREP(dIndexNormal(:))), ...
                                max(usIndexSET(dIndexNormal(:)))]);
               if ullNScanToRead > 1
                  bChange = [dIndexNormal(2:end) > 0 ...
                            usIndexLIN(dIndexNormal(2:end)) ~= usIndexLIN(dIndexNormal(1:end-1)) ...
                            usIndexSLC(dIndexNormal(2:end)) ~= usIndexSLC(dIndexNormal(1:end-1)) ...
                            usIndexPAR(dIndexNormal(2:end)) ~= usIndexPAR(dIndexNormal(1:end-1)) ...
                            usIndexACQ(dIndexNormal(2:end)) ~= usIndexACQ(dIndexNormal(1:end-1)) ...
                            usIndexECO(dIndexNormal(2:end)) ~= usIndexECO(dIndexNormal(1:end-1)) ...
                            usIndexPHS(dIndexNormal(2:end)) ~= usIndexPHS(dIndexNormal(1:end-1)) ...
                            usIndexREP(dIndexNormal(2:end)) ~= usIndexREP(dIndexNormal(1:end-1)) ...
                            usIndexSET(dIndexNormal(2:end)) ~= usIndexSET(dIndexNormal(1:end-1))];
                  dRankNormal = dRankNormal + sum(bChange, 1);
               end
               if sMDHBlock(dIBlock).usNChannelUsed < 2
                  dRankNormal(1)=0;
               end
            end

            ullIScanRead = ullIScanRead + ullNScanToRead;
            ullNScanRemainInBlock = sMDHBlock(dIBlock).ullNScan - ullIScanRead;
            ullNScanRemainInMeas = ullNScanInMeas - ullIScanRead;
            ullNScanToRead = min(ullNScanToRead,ullNScanRemainInBlock);

         end % while over subsets of blocked lines

      end % for over blocks of non-syncdata

      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.5, hProg, 'MDH Scan complete');
      end
      clear cData bChange;
      clear ulEvalMask1 usCutOffDataPre usCutOffDataPost ;
      clear bIsACQEND bIsREFLECT bIsNOISEADJSCAN bIsPHASCOR bIsRTFEEDBACK bIsRAWDATACORRECTON;
      clear usChannelID usIndexLIN usIndexACQ usIndexSLC usIndexPAR usIndexECO usIndexPHS usIndexREP usIndexSET usIndexSEG;
      clear dIndexNOISEADJSCAN dIndexPHASCOR dIndexRTFEEDBACK dIndexNormal;
      clear ullIScanRead ullNScanRemainInMeas ullNScanToRead;

      usDimPHASCOR = usDimPHASCOR .* uint16(dRankPHASCOR > 1) + 1;
      usDimNOISEADJSCAN = usDimNOISEADJSCAN .* uint16(dRankNOISEADJSCAN > 1) + 1;
      usDimRTFEEDBACK = usDimRTFEEDBACK .* uint16(dRankRTFEEDBACK > 1) + 1;
      usDimNormal = usDimNormal .* uint16(dRankNormal > 1) + 1;
      dNSegPHASCOR = dNSegPHASCOR + 1;
      dNSegNOISEADJSCAN = dNSegNOISEADJSCAN + 1;
      dNSegRTFEEDBACK = dNSegRTFEEDBACK + 1;
      dNSegNormal = dNSegNormal + 1;

      dNChaPHASCOR      = usDimPHASCOR(1);
      dNChaRTFEEDBACK   = usDimRTFEEDBACK(1);
      dNChaNOISEADJSCAN = usDimNOISEADJSCAN(1);
      dNChaNormal       = usDimNormal(1);
      dNLinPHASCOR      = usDimPHASCOR(2);
      dNLinRTFEEDBACK   = usDimRTFEEDBACK(2);
      dNLinNOISEADJSCAN = usDimNOISEADJSCAN(2);
      dNLinNormal       = usDimNormal(2);
      if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bFillToFullFourier == OFF
         if (dNLinPHASCOR > 1)
            dNLinPHASCOR      = max(dNLinPHASCOR,sDim.sRecon.dNoOfFourierLines);
         end
         if (dNLinRTFEEDBACK > 1)
            dNLinRTFEEDBACK   = max(dNLinRTFEEDBACK,sDim.sRecon.dNoOfFourierLines);
         end
         if (dNLinNOISEADJSCAN > 1)
            dNLinNOISEADJSCAN = max(dNLinNOISEADJSCAN,sDim.sRecon.dNoOfFourierLines);
         end
         if (dNLinNormal > 1)
            dNLinNormal       = max(dNLinNormal,sDim.sRecon.dNoOfFourierLines);
         end
      else
         if (dNLinPHASCOR > 1)
            dNLinPHASCOR      = max(dNLinPHASCOR,sDim.sRecon.dLIN);
         end
         if (dNLinRTFEEDBACK > 1)
            dNLinRTFEEDBACK   = max(dNLinRTFEEDBACK,sDim.sRecon.dLIN);
         end
         if (dNLinNOISEADJSCAN > 1)
            dNLinNOISEADJSCAN = max(dNLinNOISEADJSCAN,sDim.sRecon.dLIN);
         end
         if dNLinNormal > 1
            dNLinNormal       = max(dNLinNormal,sDim.sRecon.dLIN);
         end
      end
      dNSlcPHASCOR      = usDimPHASCOR(3);
      dNSlcRTFEEDBACK   = usDimRTFEEDBACK(3);
      dNSlcNOISEADJSCAN = usDimNOISEADJSCAN(3);
      dNSlcNormal       = usDimNormal(3);
      if dNSlcPHASCOR > 1
         dNSlcPHASCOR      = max(dNSlcPHASCOR,sDim.sRecon.dSLC);
      end
      if dNSlcRTFEEDBACK > 1
         dNSlcRTFEEDBACK   = max(dNSlcRTFEEDBACK,sDim.sRecon.dSLC);
      end
      if dNSlcNOISEADJSCAN > 1
         dNSlcNOISEADJSCAN = max(dNSlcNOISEADJSCAN,sDim.sRecon.dSLC);
      end
      if dNSlcNormal > 1
         dNSlcNormal       = max(dNSlcNormal,sDim.sRecon.dSLC);
      end
      dNParPHASCOR      = usDimPHASCOR(4);
      dNParRTFEEDBACK   = usDimRTFEEDBACK(4);
      dNParNOISEADJSCAN = usDimNOISEADJSCAN(4);
      dNParNormal       = usDimNormal(4);
      if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bFillToFullFourier == OFF
         if dNParPHASCOR > 1
            dNParPHASCOR      = max(dNParPHASCOR,sDim.sRecon.dNoOfFourierPartitions);
         end
         if dNParRTFEEDBACK > 1
            dNParRTFEEDBACK   = max(dNParRTFEEDBACK,sDim.sRecon.dNoOfFourierPartitions);
         end
         if dNParNOISEADJSCAN > 1
            dNParNOISEADJSCAN = max(dNParNOISEADJSCAN,sDim.sRecon.dNoOfFourierPartitions);
         end
         if dNParNormal > 1
            dNParNormal       = max(dNParNormal,sDim.sRecon.dNoOfFourierPartitions);
         end
      else
         if dNParPHASCOR > 1
            dNParPHASCOR      = max(dNParPHASCOR,sDim.sRecon.dPAR);
         end
         if dNParRTFEEDBACK > 1
            dNParRTFEEDBACK   = max(dNParRTFEEDBACK,sDim.sRecon.dPAR);
         end
         if dNParNOISEADJSCAN > 1
            dNParNOISEADJSCAN = max(dNParNOISEADJSCAN,sDim.sRecon.dPAR);
         end
         if dNParNormal > 1
            dNParNormal       = max(dNParNormal,sDim.sRecon.dPAR);
         end
      end
      dNAcqPHASCOR      = usDimPHASCOR(5);
      dNAcqRTFEEDBACK   = usDimRTFEEDBACK(5);
      dNAcqNOISEADJSCAN = usDimNOISEADJSCAN(5);
      dNAcqNormal       = usDimNormal(5);
      if dNAcqPHASCOR > 1
         dNAcqPHASCOR      = max(dNAcqPHASCOR,sDim.sRaw.dACQ);
      end
      if dNAcqRTFEEDBACK > 1
         dNAcqRTFEEDBACK   = max(dNAcqRTFEEDBACK,sDim.sRaw.dACQ);
      end
      if dNAcqNOISEADJSCAN > 1
         dNAcqNOISEADJSCAN = max(dNAcqNOISEADJSCAN,sDim.sRaw.dACQ);
      end
      if dNAcqNormal > 1
         dNAcqNormal       = max(dNAcqNormal,sDim.sRaw.dACQ);
      end
      dNEcoPHASCOR      = usDimPHASCOR(6);
      dNEcoRTFEEDBACK   = usDimRTFEEDBACK(6);
      dNEcoNOISEADJSCAN = usDimNOISEADJSCAN(6);
      dNEcoNormal       = usDimNormal(6);
      if dNEcoPHASCOR > 1
         dNEcoPHASCOR      = max(dNEcoPHASCOR,sDim.sRaw.dECO);
      end
      if dNEcoRTFEEDBACK > 1
         dNEcoRTFEEDBACK   = max(dNEcoRTFEEDBACK,sDim.sRaw.dECO);
      end
      if dNEcoNOISEADJSCAN > 1
         dNEcoNOISEADJSCAN = max(dNEcoNOISEADJSCAN,sDim.sRaw.dECO);
      end
      if dNEcoNormal > 1
         dNEcoNormal       = max(dNEcoNormal,sDim.sRaw.dECO);
      end
      dNPhsPHASCOR      = usDimPHASCOR(7);
      dNPhsRTFEEDBACK   = usDimRTFEEDBACK(7);
      dNPhsNOISEADJSCAN = usDimNOISEADJSCAN(7);
      dNPhsNormal       = usDimNormal(7);
      if dNPhsPHASCOR > 1
         dNPhsPHASCOR      = max(dNPhsPHASCOR,sDim.sRaw.dPHS);
      end
      if dNPhsRTFEEDBACK > 1
         dNPhsRTFEEDBACK   = max(dNPhsRTFEEDBACK,sDim.sRaw.dPHS);
      end
      if dNPhsNOISEADJSCAN > 1
         dNPhsNOISEADJSCAN = max(dNPhsNOISEADJSCAN,sDim.sRaw.dPHS);
      end
      if dNPhsNormal > 1
         dNPhsNormal       = max(dNPhsNormal,sDim.sRaw.dPHS);
      end
      dNRepPHASCOR      = usDimPHASCOR(8);
      dNRepRTFEEDBACK   = usDimRTFEEDBACK(8);
      dNRepNOISEADJSCAN = usDimNOISEADJSCAN(8);
      dNRepNormal       = usDimNormal(8);
      if dNRepPHASCOR > 1
         dNRepPHASCOR      = max(dNRepPHASCOR,sDim.sRaw.dREP);
      end
      if dNRepRTFEEDBACK > 1
         dNRepRTFEEDBACK   = max(dNRepRTFEEDBACK,sDim.sRaw.dREP);
      end
      if dNRepNOISEADJSCAN > 1
         dNRepNOISEADJSCAN = max(dNRepNOISEADJSCAN,sDim.sRaw.dREP);
      end
      if dNRepNormal > 1
         dNRepNormal       = max(dNRepNormal,sDim.sRaw.dREP);
      end
      dNSetPHASCOR      = usDimPHASCOR(9);
      dNSetRTFEEDBACK   = usDimRTFEEDBACK(9);
      dNSetNOISEADJSCAN = usDimNOISEADJSCAN(9);
      dNSetNormal       = usDimNormal(9);
      if dNSetPHASCOR > 1
         dNSetPHASCOR      = max(dNSetPHASCOR,sDim.sRaw.dSET);
      end
      if dNSetRTFEEDBACK > 1
         dNSetRTFEEDBACK   = max(dNSetRTFEEDBACK,sDim.sRaw.dSET);
      end
      if dNSetNOISEADJSCAN > 1
         dNSetNOISEADJSCAN = max(dNSetNOISEADJSCAN,sDim.sRaw.dSET);
      end
      if dNSetNormal > 1
         dNSetNormal       = max(dNSetNormal,sDim.sRaw.dSET);
      end

      dDimOutPHASCOR = uint64([dNChaPHASCOR dNLinPHASCOR dNSlcPHASCOR dNParPHASCOR dNAcqPHASCOR dNEcoPHASCOR dNPhsPHASCOR dNRepPHASCOR dNSetPHASCOR]);
      dDimOutNOISEADJSCAN = uint64([dNChaNOISEADJSCAN dNLinNOISEADJSCAN dNSlcNOISEADJSCAN dNParNOISEADJSCAN dNAcqNOISEADJSCAN dNEcoNOISEADJSCAN dNPhsNOISEADJSCAN dNRepNOISEADJSCAN dNSetNOISEADJSCAN]);
      dDimOutRTFEEDBACK = uint64([dNChaRTFEEDBACK dNLinRTFEEDBACK dNSlcRTFEEDBACK dNParRTFEEDBACK dNAcqRTFEEDBACK dNEcoRTFEEDBACK dNPhsRTFEEDBACK dNRepRTFEEDBACK dNSetRTFEEDBACK]);
      dDimOutNormal = uint64([dNChaNormal dNLinNormal dNSlcNormal dNParNormal dNAcqNormal dNEcoNormal dNPhsNormal dNRepNormal dNSetNormal]);

      [dRankPHASCOR, orderPHASCOR] = sort(dRankPHASCOR, 2, 'descend');
      dDimOutPHASCOR = dDimOutPHASCOR(orderPHASCOR);
      dCumProdPHASCOR = [uint64(1) fastReadVD13DCumprod(dDimOutPHASCOR(1:end-1))];
      [dRankNOISEADJSCAN, orderNOISEADJSCAN] = sort(dRankNOISEADJSCAN, 2, 'descend');
      dDimOutNOISEADJSCAN = dDimOutNOISEADJSCAN(orderNOISEADJSCAN);
      dCumProdNOISEADJSCAN = [uint64(1) fastReadVD13DCumprod(dDimOutNOISEADJSCAN(1:end-1))];
      [dRankRTFEEDBACK, orderRTFEEDBACK] = sort(dRankRTFEEDBACK, 2, 'descend');
      dDimOutRTFEEDBACK = dDimOutRTFEEDBACK(orderRTFEEDBACK);
      dCumProdRTFEEDBACK = [uint64(1) fastReadVD13DCumprod(dDimOutRTFEEDBACK(1:end-1))];
      [dRankNormal, orderNormal] = sort(dRankNormal, 2, 'descend');
      dDimOutNormal = dDimOutNormal(orderNormal);
      dCumProdNormal = [uint64(1) fastReadVD13DCumprod(dDimOutNormal(1:end-1))];

      if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bFillToFullFourier == OFF
         ullNCol = ullNSamplesInScan;
      else
         ullNCol = max(ullNSamplesInScan, uint64(sDim.sRecon.dRoFTLength));
      end

      clear usDimPHASCOR usDimRTFEEDBACK usDimNOISEADJSCAN usDimNormal;
      clear dRankPHASCOR dRankRTFEEDBACK dRankNOISEADJSCAN dRankNormal;
      clear dNChaPHASCOR dNLinPHASCOR dNSlcPHASCOR dNParPHASCOR dNAcqPHASCOR dNEcoPHASCOR dNPhsPHASCOR dNRepPHASCOR dNSetPHASCOR;
      clear dNChaNOISEADJSCAN dNLinNOISEADJSCAN dNSlcNOISEADJSCAN dNParNOISEADJSCAN dNAcqNOISEADJSCAN dNEcoNOISEADJSCAN dNPhsNOISEADJSCAN dNRepNOISEADJSCAN dNSetNOISEADJSCAN;
      clear dNChaRTFEEDBACK dNLinRTFEEDBACK dNSlcRTFEEDBACK dNParRTFEEDBACK dNAcqRTFEEDBACK dNEcoRTFEEDBACK dNPhsRTFEEDBACK dNRepRTFEEDBACK dNSetRTFEEDBACK;
      clear dNChaNormal dNLinNormal dNSlcNormal dNParNormal dNAcqNormal dNEcoNormal dNPhsNormal dNRepNormal dNSetNormal;

      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(1.0, hProg, 'Output dimensions computed.');
      end

      % --------------------------------------------------------------------------
      % Third pass to read the data
      % --------------------------------------------------------------------------
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.0, hProg, 'Allocating arrays');
      end

      % Allocate arrays
      if dNScanPHASCOR > 0
         sCon.sUser.bFoundPhaseCor = 1;
         kPhaCor{dIMeas}.data = complex(zeros(ullNCol, fastReadVD13DProd(dDimOutPHASCOR), dNSegPHASCOR, 'single'));
      end
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.25, hProg, 'Allocating arrays');
      end
      if dNScanNOISEADJSCAN > 0
         sCon.sUser.bFoundNoiseAdj = 1;
         kNoiseAdj{dIMeas}.data = complex(zeros(ullNCol, fastReadVD13DProd(dDimOutNOISEADJSCAN), 'single'));
         if dNSegNOISEADJSCAN > 1
            kNoiseAdj{dIMeas}.seg = zeros(1,fastReadVD13DProd(dDimOutNOISEADJSCAN),'uint16');
         end
      end
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.50, hProg, 'Allocating arrays');
      end
      if dNScanRTFEEDBACK > 0
         sCon.sUser.bFoundRTFeedback = 1;
         kRTFeedback{dIMeas}.data = complex(zeros(ullNCol, fastReadVD13DProd(dDimOutRTFEEDBACK), 'single'));
         if dNSegRTFEEDBACK > 1
            kRTFeedback{dIMeas}.seg = zeros(1,fastReadVD13DProd(dDimOutRTFEEDBACK),'uint16');
         end
      end
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.75, hProg, 'Allocating arrays');
      end
      if dNScanNormal > 0
         kSpace{dIMeas}.data = complex(zeros(ullNCol, fastReadVD13DProd(dDimOutNormal), 'single'));
         if dNSegNormal > 1
            kSpace{dIMeas}.seg = zeros(1,fastReadVD13DProd(dDimOutNormal),'uint16');
         end
      end
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(1.0, hProg, 'Allocating arrays');
      end

      ullNScanRemainInMeas = ullNScanInMeas;
      ullIScanRead = uint64(0);
      time0 = tic;

      for dIBlock=1:length(sMDHBlock)

         % Skip SYNCDATA
         if sMDHBlock(dIBlock).bIsSYNCDATA || sMDHBlock(dIBlock).bIsACQEND
            continue
         end

         ullMaxScanPerRead = max(idivide(ullMaxBytesPerRead,sMDHBlock(dIBlock).ullLength,'fix'),uint64(1));
         ullNScanToRead = min(ullMaxScanPerRead,sMDHBlock(dIBlock).ullNScan);

         st = fseek(fid,sMDHBlock(dIBlock).ullOffset,'bof');
         ullNScanRemainInBlock = sMDHBlock(dIBlock).ullNScan;

         while ullNScanRemainInBlock > 0

            dTimeElapsed = toc(time0);
            rate = dTimeElapsed / (double(ullIScanRead) + 1);
            timeRemain = double(ullNScanInMeas - ullIScanRead) * rate;
            if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
               cWinMsg = sprintf('Pass 3 of 3\nReading scans\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(ullIScanRead)/double(ullNScanInMeas), hProg, cWinMsg);
            end
            cData = fread(fid,[sMDHBlock(dIBlock).ullLength ullNScanToRead],'uchar=>uchar');
            ulEvalMask1 = cData(41:44,:);
            ulEvalMask1 = typecast(ulEvalMask1(:),'uint32');
            usCutOffDataPre = cData(81:82,:);
            usCutOffDataPre = typecast(usCutOffDataPre(:),'uint16');
            usCutOffDataPost = cData(83:84,:);
            usCutOffDataPost = typecast(usCutOffDataPost(:),'uint16');
            bIsACQEND = bitand(ulEvalMask1,ulMDH_ACQEND) == ulMDH_ACQEND;
            bIsNOISEADJSCAN = bitand(ulEvalMask1,ulMDH_NOISEADJSCAN) == ulMDH_NOISEADJSCAN;
            bIsPHASCOR = bitand(ulEvalMask1,ulMDH_PHASCOR) == ulMDH_PHASCOR;
            bIsRTFEEDBACK = bitand(ulEvalMask1,ulMDH_RTFEEDBACK) == ulMDH_RTFEEDBACK;
            bIsRAWDATACORRECTION = bitand(ulEvalMask1,ulMDH_RAWDATACORRECTION) == ulMDH_RAWDATACORRECTION;
            usIndexLIN = cData(53:54,:);
            usIndexLIN = uint64(typecast(usIndexLIN(:),'uint16'));
            usIndexACQ = cData(55:56,:);
            usIndexACQ = uint64(typecast(usIndexACQ(:),'uint16'));
            usIndexSLC = cData(57:58,:);
            usIndexSLC = uint64(typecast(usIndexSLC(:),'uint16'));
            usIndexPAR = cData(59:60,:);
            usIndexPAR = uint64(typecast(usIndexPAR(:),'uint16'));
            usIndexECO = cData(61:62,:);
            usIndexECO = uint64(typecast(usIndexECO(:),'uint16'));
            usIndexPHS = cData(63:64,:);
            usIndexPHS = uint64(typecast(usIndexPHS(:),'uint16'));
            usIndexREP = cData(65:66,:);
            usIndexREP = uint64(typecast(usIndexREP(:),'uint16'));
            usIndexSET = cData(67:68,:);
            usIndexSET = uint64(typecast(usIndexSET(:),'uint16'));
            usIndexSEG = cData(69:70,:);
            usIndexSEG = uint64(typecast(usIndexSEG(:),'uint16'));
            usChannelID = zeros(2,sMDHBlock(dIBlock).usNChannelUsed,ullNScanToRead,'uint8');
            for dICHA = 1:sMDHBlock(dIBlock).usNChannelUsed
               ullOffset  = uint64(ullMDHScanSize) + uint64(dICHA-1) * ullSizeChanInBytes;
               usChannelID(:,dICHA,:) = reshape(cData(ullOffset+25:ullOffset+26,:),[2,1,ullNScanToRead]);
            end
            usChannelID = reshape(typecast(usChannelID(:),'uint16'),sMDHBlock(dIBlock).usNChannelUsed,ullNScanToRead);
            ullSizeSamplesInBytes = sMDHBlock(dIBlock).ullNSamplesInScan * uint64(8);
            for dICHA = 1:sMDHBlock(dIBlock).usNChannelUsed

               usChannelIDOut = squeeze(usChannelID(dICHA,:));
               ullOffset  = uint64(ullMDHScanSize) + uint64(dICHA-1) * ullSizeChanInBytes + ullMDHChanSize + uint64(1);

               tempArr = cData(ullOffset:ullOffset+ullSizeSamplesInBytes-1,:);
               tempArr = reshape(typecast(tempArr(:),'single'), [2*sMDHBlock(dIBlock).ullNSamplesInScan ullNScanToRead]);
               tempArr = complex(tempArr(1:2:end,:),tempArr(2:2:end,:));

               % Handle all the reflections at once
               bIsREFLECT = bitand(ulEvalMask1,ulMDH_REFLECT) == ulMDH_REFLECT;
               indexReflect = find(bIsREFLECT);
               if ~isempty(indexReflect)
                  tempArr(:,indexReflect) = flipud(tempArr(:,indexReflect));
               end

               % Handle Pre-cutoff
               indexPre = find(usCutOffDataPre > 0 & ~bIsNOISEADJSCAN);
               if ~isempty(indexPre)
                  usCutOffDataPreUniq = unique(usCutOffDataPre(indexPre));
                  for iUniq = 1:numel(usCutOffDataPreUniq)
                     indexCut = find(usCutOffDataPre(indexPre) == usCutOffDataPreUniq(iUniq));
                     tempArr(1:usCutOffDataPreUniq(iUniq),indexPre(indexCut)) = ...
                        complex(zeros(usCutOffDataPreUniq(iUniq),numel(indexCut),'single'));
                  end % loop over uniq cuts
               end

               % Handle Post-cutoff
               indexPost = find(usCutOffDataPost > 0 & ~bIsNOISEADJSCAN);
               if ~isempty(indexPost)
                  usCutOffDataPostUniq = unique(usCutOffDataPost(indexPost));
                  for iUniq = 1:numel(usCutOffDataPostUniq)
                     indexCut = find(usCutOffDataPost(indexPost) == usCutOffDataPostUniq(iUniq));
                     tempArr(1:usCutOffDataPostUniq(iUniq),indexPost(indexCut)) = ...
                        complex(zeros(usCutOffDataPostUniq(iUniq),numel(indexCut),'single'));
                  end % loop over uniq cuts
               end

               % Phase cor scans
               dIndexPHASCOR = find(bIsPHASCOR & ~bIsACQEND);
               if ~isempty(dIndexPHASCOR)
                  if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadPhaseCor == ON && sDim.sRaw.dNPhCorScan > 0
                     for kScan = 1:numel(dIndexPHASCOR)
                        jScan = dIndexPHASCOR(kScan);
                        usIndexCHA = find(usChannelIDUniq == usChannelIDOut(jScan)) - 1;
                        usIndexOut = [usIndexCHA usIndexLIN(jScan) usIndexSLC(jScan) usIndexPAR(jScan) usIndexACQ(jScan) usIndexECO(jScan) usIndexPHS(jScan) usIndexREP(jScan) usIndexSET(jScan)];
                        temp = tempArr(:,jScan);
                        usIndexOutPhaCor = mod(usIndexOut(orderPHASCOR),dDimOutPHASCOR);
                        iScanPhaCor = sum(usIndexOutPhaCor(:) .* dCumProdPHASCOR(:)) + 1;
                        kPhaCor{dIMeas}.data(1:sMDHBlock(dIBlock).ullNSamplesInScan,iScanPhaCor,usIndexSEG(jScan)+1) = temp;
                     end
                  end
               end

               % RT Feedback
               dIndexRTFEEDBACK = find(bIsRTFEEDBACK & ~bIsACQEND);
               if ~isempty(dIndexRTFEEDBACK)
                  if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadRTFeedback == ON 
                     for kScan = 1:numel(dIndexRTFEEDBACK)
                        jScan = dIndexRTFEEDBACK(kScan);
                        usIndexCHA = find(usChannelIDUniq == usChannelIDOut(jScan)) - 1;
                        usIndexOut = [usIndexCHA usIndexLIN(jScan) usIndexSLC(jScan) usIndexPAR(jScan) usIndexACQ(jScan) usIndexECO(jScan) usIndexPHS(jScan) usIndexREP(jScan) usIndexSET(jScan)];
                        temp = tempArr(:,jScan);
                        usIndexOutRTFeedback = mod(usIndexOut(orderRTFEEDBACK),dDimOutRTFEEDBACK);
                        iScanRTFeedback = sum(usIndexOutRTFeedback(:) .* dCumProdRTFEEDBACK(:)) + 1;
                        kRTFeedback{dIMeas}.data(1:sMDHBlock(dIBlock).ullNSamplesInScan,iScanRTFeedback) = temp;
                        if dNSegRTFEEDBACK > 1
                           kRTFeedback{dIMeas}.seg(1,iScanRTFeedback) = usIndexSEG(jScan) + 1;
                        end
                     end
                  end
               end

               % Noise adj scans
               dIndexNOISEADJSCAN = find(bIsNOISEADJSCAN & ~bIsACQEND);
               if ~isempty(dIndexNOISEADJSCAN)
                  if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bReadNoiseAdj == ON 
                     for kScan = 1:numel(dIndexNOISEADJSCAN)
                        jScan = dIndexNOISEADJSCAN(kScan);
                        %usIndexCHA = find(usChannelIDUniq == usChannelIDOut(jScan)) - 1;
                        usIndexOutNoiseAdj = mod(usIndexOut(orderNOISEADJSCAN),dDimOutNOISEADJSCAN);
                        iScanNoiseAdj = sum(usIndexOutNoiseAdj(:) .* dCumProdNOISEADJSCAN(:)) + 1;
                        kNoiseAdj{dIMeas}.data(1:sMDHBlock(dIBlock).ullNSamplesInScan,iScanNoiseAdj) = temp;
                        if dNSegNOISEADJSCAN > 1
                           kNoiseAdj{dIMeas}.seg(1,iScanNoiseAdj) = usIndexSEG(jScan) + 1;
                        end
                     end
                  end
               end

               % Normal lines
               dIndexNormal = find(~(bIsPHASCOR | bIsRTFEEDBACK | bIsNOISEADJSCAN | bIsACQEND ));
               for kScan = 1:numel(dIndexNormal)
                  jScan = dIndexNormal(kScan);
                  keyCHACoilSelectMeas = usChannelIDOut(jScan) + 1;
                  keyCHACoilSelect = mCoilSelectMeas(keyCHACoilSelectMeas).tElement;
                  usIndexCHA = find(usChannelIDUniq == (keyCHACoilSelectMeas-1)) - 1;
                  usIndexOut = [usIndexCHA usIndexLIN(jScan) usIndexSLC(jScan) usIndexPAR(jScan) usIndexACQ(jScan) usIndexECO(jScan) usIndexPHS(jScan) usIndexREP(jScan) usIndexSET(jScan)];
                  temp = tempArr(:,jScan);
                  usIndexOutNormal = mod(usIndexOut(orderNormal),dDimOutNormal);
                  iScanNormal = sum(usIndexOutNormal(:) .* dCumProdNormal(:)) + 1;
                  if bIsRAWDATACORRECTION(jScan) && sCon.sUser.bApplyRawDataCorrection
                     temp = mCoilSelect(keyCHACoilSelect).dRawDataCorrectionFactor*temp;
                  end
                  kSpace{dIMeas}.data(1:sMDHBlock(dIBlock).ullNSamplesInScan,iScanNormal) = mCoilSelectMeas(keyCHACoilSelectMeas).flFFTCorrectionFactor*temp;
                  if dNSegNormal > 1
                     kSpace{dIMeas}.seg(1,iScanNormal) = usIndexSEG(jScan) + 1;
                  end
               end %endfor over this block of lines

            end % endfor over channels

            ullIScanRead = ullIScanRead + ullNScanToRead;
            ullNScanRemainInBlock = sMDHBlock(dIBlock).ullNScan - ullIScanRead;
            ullNScanRemainInMeas = ullNScanInMeas - ullIScanRead;
            ullNScanToRead = min(ullNScanToRead,ullNScanRemainInBlock);

         end % while over sub-blocks of lines
      end % for over sMDHBlocks
      fclose(fid);
      clear cData;
      clear dIndexNOISEADJSCAN dIndexPHASCOR dIndexRTFEEDBACK dIndexNormal;
      clear ulEvalMask1 usCutOffDataPre usCutOffDataPost tempArr;
      clear bIsACQEND bIsREFLECT bIsNOISEADJSCAN bIsPHASCOR bIsRTFEEDBACK bIsRAWDATACORRECTON;
      clear usChannelID usChannelIDOut usIndexLIN usIndexACQ usIndexSLC usIndexPAR usIndexECO usIndexPHS usIndexREP usIndexSET usIndexSEG;
      clear ullIScanRead ullNScanRemainInMeas ullNScanToRead;
      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(0.5, hProg, 'All lines read from file.  Reshaping.');
      end

      % Reshape the results
      if dNScanPHASCOR > 0
         indexDimNonZeroPHASCOR = find(dDimOutPHASCOR > 1);
         kPhaCor{dIMeas}.data = reshape(kPhaCor{dIMeas}.data, [ullNCol dDimOutPHASCOR(indexDimNonZeroPHASCOR) dNSegPHASCOR]);
      end
      if dNScanNOISEADJSCAN > 0
         indexDimNonZeroNOISEADJSCAN = find(dDimOutNOISEADJSCAN > 1);
         kNoiseAdj{dIMeas}.data = reshape(kNoiseAdj{dIMeas}.data, [ullNCol dDimOutNOISEADJSCAN(indexDimNonZeroNOISEADJSCAN)]);
         if dNSegNOISEADJSCAN > 1
            kNoiseAdj{dIMeas}.seg = squeeze(reshape(kNoiseAdj{dIMeas}.seg, [1 dDimOutNOISEADJSCAN(indexDimNonZeroNOISEADJSCAN)]));
         end
      end
      if dNScanRTFEEDBACK > 0
         indexDimNonZeroRTFEEDBACK = find(dDimOutRTFEEDBACK > 1);
         kRTFeedback{dIMeas}.data = reshape(kRTFeedback{dIMeas}.data, [ullNCol dDimOutRTFEEDBACK(indexDimNonZeroRTFEEDBACK)]);
         if dNSegRTFEEDBACK > 1
            kRTFeedback{dIMeas}.seg = squeeze(reshape(kRTFeedback{dIMeas}.seg, [1 dDimOutRTFEEDBACK(indexDimNonZeroRTFEEDBACK)]));
         end
      end
      if dNScanNormal > 0
         indexDimNonZeroNormal = find(dDimOutNormal > 1);
         kSpace{dIMeas}.data = reshape(kSpace{dIMeas}.data, [ullNCol dDimOutNormal(indexDimNonZeroNormal)]);
         if dNSegNormal > 1
            kSpace{dIMeas}.seg = squeeze(reshape(kSpace{dIMeas}.seg, [1 dDimOutNormal(indexDimNonZeroNormal)]));
         end
      end

      if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bShiftDCToMatrixCenter == ON
         if dNScanPHASCOR > 0
            allShift = dDimOutPHASCOR * 0;
            colShift = idivide(ullNCol - ullNSamplesInScan, uint64(2), 'fix');
            if dDimOutPHASCOR(2) > 1
               linShift = uint64(dDimOutPHASCOR(2)/2) - uint64(KSpaceCentreLineNo);
               allShift(2) = linShift;
            end
            if dDimOutPHASCOR(4) > 1
               parShift = uint64(dDimOutPHASCOR(4)/2) - uint64(KSpaceCentrePartitionNo);
               allShift(4) = parShift;
            end
            allShift = allShift(orderPHASCOR);
            allShift = [colShift allShift(indexDimNonZeroPHASCOR) 0];
            kPhaCor{dIMeas}.data = circshift(kPhaCor{dIMeas}.data,allShift);
         end
         if dNScanNOISEADJSCAN > 0
            allShift = dDimOutNOISEADJSCAN * 0;
            colShift = idivide(ullNCol - ullNSamplesInScan, uint64(2), 'fix');
            if dDimOutNOISEADJSCAN(2) > 1
               linShift = uint64(dDimOutNOISEADJSCAN(2)/2) - uint64(KSpaceCentreLineNo);
               allShift(2) = linShift;
            end
            if dDimOutNOISEADJSCAN(4) > 1
               parShift = uint64(dDimOutNOISEADJSCAN(4)/2) - uint64(KSpaceCentrePartitionNo);
               allShift(4) = parShift;
            end
            allShift = allShift(orderNOISEADJSCAN);
            allShift = [colShift allShift(indexDimNonZeroNOISEADJSCAN)];
            kNoiseAdj{dIMeas}.data = circshift(kNoiseAdj{dIMeas}.data,allShift);
         end
         if dNScanRTFEEDBACK > 0
            allShift = dDimOutRTFEEDBACK * 0;
            colShift = idivide(ullNCol - ullNSamplesInScan, uint64(2), 'fix');
            if dDimOutRTFEEDBACK(2) > 1
               linShift = uint64(dDimOutRTFEEDBACK(2)/2) - uint64(KSpaceCentreLineNo);
               allShift(2) = linShift;
            end
            if dDimOutRTFEEDBACK(4) > 1
               parShift = uint64(dDimOutRTFEEDBACK(4)/2) - uint64(KSpaceCentrePartitionNo);
               allShift(4) = parShift;
            end
            allShift = allShift(orderRTFEEDBACK);
            allShift = [colShift allShift(indexDimNonZeroRTFEEDBACK)];
            kRTFeedback{dIMeas}.data = circshift(kRTFeedback{dIMeas}.data,allShift);
         end
         if dNScanNormal > 0
            allShift = dDimOutNormal * 0;
            colShift = idivide(ullNCol - ullNSamplesInScan, uint64(2), 'fix');
            if dDimOutNormal(2) > 1
               linShift = uint64(dDimOutNormal(2)/2) - uint64(KSpaceCentreLineNo);
               allShift(2) = linShift;
            end
            if dDimOutNormal(4) > 1
               parShift = uint64(dDimOutNormal(4)/2) - uint64(KSpaceCentrePartitionNo);
               allShift(4) = parShift;
            end
            allShift = allShift(orderNormal);
            allShift = [colShift allShift(indexDimNonZeroNormal)];
            kSpace{dIMeas}.data = circshift(kSpace{dIMeas}.data,allShift);
         end
      end

      if ~(sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI || sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent)
         waitbar(1.0, hProg, 'kSpace complete.');
         close(hProg);
      end

      if dNScanNormal > 0
         tMsg = sprintf('Normal kSpace dimensions:\n   col=%d', ullNCol);
         cDimNameNormal = cDimName(orderNormal);
         for iDim=1:numel(indexDimNonZeroNormal)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        cDimNameNormal{indexDimNonZeroNormal(iDim)}, ...
                                        dDimOutNormal(indexDimNonZeroNormal(iDim))));
         end
         if (dNSegNormal > 1) 
            tMsg = strcat(tMsg, sprintf('\n   seg=%d', dNSegNormal));
         end
         kSpace{dIMeas}.dim = ['col', cDimNameNormal(indexDimNonZeroNormal)];
      end

      if dNScanPHASCOR > 0
         tMsg = strcat(tMsg, sprintf('\nPhase Cor dimensions:\n   col=%d', ullNCol));
         cDimNamePHASCOR = cDimName(orderPHASCOR);
         for iDim=1:numel(indexDimNonZeroPHASCOR)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        cDimNamePHASCOR{indexDimNonZeroPHASCOR(iDim)}, ...
                                        dDimOutPHASCOR(indexDimNonZeroPHASCOR(iDim))));
         end
         tMsg = strcat(tMsg, sprintf('\n   seg=%d', dNSegPHASCOR));
         kPhaCor{dIMeas}.dim = ['col', cDimNamePHASCOR(indexDimNonZeroPHASCOR), 'seg'];
      end

      if ~sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bSilent
         if sMrParcRaidFileEntry{dIMeas}.sCon.sUser.bNoGUI
            tMsg
         else
            uiwait(msgbox(tMsg, 'Output kSpace Dimensions','modal'));
         end
      end

      clear dim sProtocol mCoilSelect mCoilSelectMeas

   end % endfor over all measurements in file
   for dIArg=1:nargout
      varargout{dIArg} = 0;
   end
   if sCon.sUser.bReadPMU
      varargout{sCon.sUser.iPMUvarargout} = PMUOut;
   end
   if and(sCon.sUser.bFoundPhaseCor,sCon.sUser.bReadPhaseCor)
      varargout{sCon.sUser.iPhaseCorvarargout} = kPhaCor;
   end
   if and(sCon.sUser.bFoundNoiseAdj,sCon.sUser.bReadNoiseAdj)
      varargout{sCon.sUser.iNoiseAdjvarargout} = kNoiseADj;
   end
   if and(sCon.sUser.bFoundRTFeedback,sCon.sUser.bReadRTFeedback)
      varargout{sCon.sUser.iRTFeedbackvarargout} = kRTFeedback;
   end
   if and(sCon.sUser.bFoundRelSliceNumber,sCon.sUser.bReadRelSliceNumber)
      varargout{sCon.sUser.iRelSliceNumbervarargout} = RelativeSliceNumber;
   end

end

% =============================================================================
function output = fastReadVD13DCumprod(input)
% -----------------------------------------------------------------------------
% fastReadVD13DCumprod  Compute the cumulative product for an input array.
%                          Apparently before matlab 2013, cumprod only supports
%                          floating cumprod and I need one for integers!
%
% output = fastReadVD13DCUMPROD(input) returns array of same type and structure
%                                        as input with cumulative product 
%                                        computed from first to last element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = ones(size(input),class(input(1)));
   output(1) = input(1);
   for p=2:numel(input)
      output(p) = output(p-1) * input(p);
   end

end

% =============================================================================
function output = fastReadVD13DProd(input)
% -----------------------------------------------------------------------------
% fastReadVD13DCumprod  Compute the product for an input array.
%                          Apparently before matlab 2013, prod only supports
%                          floating prod and I need one for integers!
%
% output = fastReadVD13DPROD(input) returns array of same type and structure
%                                     as input with product 
%                                     computed from first to last element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = input(1);
   for p=2:numel(input)
      output = output * input(p);
   end

end



