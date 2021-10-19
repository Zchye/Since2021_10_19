classdef hNRSchedulingLogger < handle
%hNRSchedulingLogger Scheduler logging and visualization
%   The class implements logging mechanism as well as visualization of the
%   logs. The following two types of visualizations are shown:
%    (i) Display of CQI values for UEs over the bandwidth
%   (ii) Display of resource grid assignment to UEs. This 2D time-frequency
%        grid shows the RB allocation to the UEs in the previous slot for
%        symbol based scheduling and previous frame for slot based
%        scheduling. HARQ process for the assignments is also shown
%        alongside the UE's RNTI

%   Copyright 2019-2021 The MathWorks, Inc.

    properties
        %NCellID Cell id to which the logging and visualization object belongs
        NCellID (1, 1) {mustBeInteger, mustBeInRange(NCellID, 0, 1007)} = 1;

        % NumUEs Count of UEs
        NumUEs

        % NumHARQ Number of HARQ processes
        % The default value 16 HARQ processes
        NumHARQ (1, 1) {mustBeInteger, mustBeInRange(NumHARQ, 1, 16)} = 16;

        % NumFrames Number of frames in simulation
        NumFrames

        % SchedulingType Type of scheduling (slot based or symbol based)
        % Value 0 means slot based and value 1 means symbol based. The
        % default value is 0
        SchedulingType (1, 1) {mustBeInteger, mustBeInRange(SchedulingType, 0, 1)} = 0;

        % DuplexMode Duplexing mode
        % Frequency division duplexing (FDD) or time division duplexing (TDD)
        % Value 0 means FDD and 1 means TDD. The default value is 0
        DuplexMode (1, 1) {mustBeInteger, mustBeInRange(DuplexMode, 0, 1)} = 0;

        % ColumnIndexMap Mapping the column names of logs to respective column indices
        % It is a map object
        ColumnIndexMap

        % GrantColumnIndexMap Mapping the column names of scheduling logs to respective column indices
        % It is a map object
        GrantLogsColumnIndexMap

        % NumRBs Number of resource blocks
        % A vector of two elements. First element represents number of
        % PDSCHRBs and second element represents PUSCHRBS
        NumRBs = zeros(2, 1);
        
        % Bandwidth Carrier bandwidth
        % A vector of two elements. First element represents downlink and
        % second element represents uplink bandwidth respectively
        Bandwidth

        % RBGSizeConfig Type of RBG table to use
        % Flag used in determining the RBGsize. Value 1 represents
        % (configuration-1 RBG table) or 2 represents (configuration-2 RBG
        % table) as defined in 3GPP TS 38.214 Section 5.1.2.2.1. The
        % default value is 1
        RBGSizeConfig = 1;

        % SchedulingLog Symbol-by-symbol log of the simulation
        % In FDD mode first element contains downlink scheduling
        % information and second element contains uplink scheduling
        % information. In TDD mode first element contains scheduling
        % information of both downlink and uplink
        SchedulingLog = cell(2, 1);

        % GrantLog Log of the scheduling grants 
        % It also contains the parameters for scheduling decisions
        GrantLog

        % IsLogReplay Flag to decide the type of post-simulation visualization
        % whether to show plain replay of the resource assignment during
        % simulation or of the selected slot (or frame). During the
        % post-simulation visualization, setting the value to 1 just
        % replays the resource assignment of the simulation frame-by-frame
        % (or slot-by-slot). Setting value to 0 gives the option to select
        % a particular frame (or slot) to see the way resources are
        % assigned in the chosen frame (or slot)
        IsLogReplay
        
        % Resource grid information related properties
        % EnableResourceGridVisualization Switch to turn on/off the resource grid visualization (resource-grid occupancy)
        EnableResourceGridVisualization = true;
        % RGMaxRBsToDisplay Max number of RBs displayed in resource grid visualization
        RGMaxRBsToDisplay = 20
        % RGMaxSlotsToDisplay Max number of slots displayed in resource grid visualization
        RGMaxSlotsToDisplay = 10

        % CQI information related properties
        % EnableCQIGridVisualization Switch to turn on/off the CQI grid visualization
        EnableCQIGridVisualization = true;
        % CVMaxRBsToDisplay Max number of RBs to be displayed in CQI visualization
        CVMaxRBsToDisplay = 20
        % CVMaxUEsToDisplay Max number of UEs to be displayed in CQI visualization
        CVMaxUEsToDisplay = 10

        % PeakDataRateDL Theoretical peak data rate in the downlink
        % direction
        PeakDataRateDL

        % PeakDataRateUL Theoretical peak data rate in the uplink direction
        PeakDataRateUL
    end
    
    properties (GetAccess = public, SetAccess = private)
        % UEIdList RNTIs of UEs in a cell as row vector
        UEIdList
    end

    properties (Constant)
        %NumSym Number of symbols in a slot
        NumSym = 14;

        % NominalRBGSizePerBW Nominal RBG size table
        % It is for the specified bandwidth in accordance with 
        % 3GPP TS 38.214, Section 5.1.2.2.1
        NominalRBGSizePerBW = [
            36   2   4
            72   4   8
            144  8   16
            275  16  16
            ];

        % Duplexing mode related constants
        % FDDDuplexMode Frequency division duplexing mode
        FDDDuplexMode = 0;
        % TDDDuplexMode Time division duplexing mode
        TDDDuplexMode = 1;

        % Constants related to scheduling type
        % SymbolBased Symbol based scheduling
        SymbolBased = 1;
        % SlotBased Slot based scheduling
        SlotBased = 0;

        % Constants related to downlink and uplink information. These
        % constants are used for indexing logs and identifying plots
        % DownlinkIdx Index for all downlink information
        DownlinkIdx = 1;
        % UplinkIdx Index for all downlink information
        UplinkIdx = 2;

        % ColorCoding Mapping of a range of CQI values to particular color
        ColorCoding = {[0.85 0.32 0.09], [0.85 0.32 0.09], [0.88 0.50 0.09], [0.88 0.50 0.09], ...
            [0.93 0.69 0.13], [0.93 0.69 0.13], [0.98 0.75 0.26], [0.98 0.75 0.26], ...
            [0.98 0.82 0.14], [0.98 0.82 0.14], [0.8 0.81 0.16], [0.8 0.81 0.16], ...
            [0.68 0.71 0.18], [0.68 0.71 0.18], [0.46 0.67 0.18], [0.46 0.67 0.18]}
    end

    properties (Access = private)
        % NumSlotsFrame Number of slots in 10ms time frame
        NumSlotsFrame

        % CurrSlot Current slot in the frame
        CurrSlot

        % CurrFrame Current frame
        CurrFrame

        % CurrSymbol Current symbol in the slot
        CurrSymbol

        % NumLogs Number of logs to be created based on number of links
        NumLogs

        % SymbolInfo Information about how each symbol (UL/DL/Guard) is allocated
        SymbolInfo

        % SlotInfo Information about how each slot (UL/DL/Guard) is allocated
        SlotInfo

        % PlotIds
        PlotIds

        % GrantCount Keeps track of count of grants sent
        GrantCount = 0
        
        % RBItemsList Items List for RBs drop down for DL and UL
        RBItemsList = cell(2, 1);

        % Resource grid information related properties
        % ResourceGrid In FDD mode first element contains downlink resource
        % grid allocation status and second element contains uplink
        % resource grid allocation status. In TDD mode first element
        % contains resource grid allocation status for downlink and uplink.
        % Each element is a 2D resource grid of N-by-P matrix where 'N' is
        % the number of slot or symbols and 'P' is the number of RBs in the
        % bandwidth to store how UEs are assigned different time-frequency
        % resources.
        ResourceGrid = cell(2, 1);
        % RGVisualizationFigHandle Handle of the resource grid visualization
        RGVisualizationFigHandle
        % ResourceGridReTxInfo First element contains transmission status
        % in downlink and second element contains transmission status in
        % uplink for FDD mode. In TDD mode first element contains
        % transmission status for both downlink and uplink. Each element is
        % a 2D resource grid of N-by-P matrix where 'N' is the number of
        % slot or symbols and 'P' is the number of RBs in the bandwidth to
        % store type:new-transmission or retransmission.
        ResourceGridReTxInfo = cell(2, 1);
        % ResourceGridHarqInfo In FDD mode first element contains downlink
        % HARQ information and second element contains uplink HARQ
        % information. In TDD mode first element contains HARQ information
        % for downlink and uplink. Each element is a 2D resource grid of
        % N-by-P matrix where 'N' is the number of slot or symbols and 'P'
        % is the number of RBs in the bandwidth to store the HARQ process
        ResourceGridHarqInfo
        % ResourceGridTextHandles Text handles to display of the RNTI for the RBs
        ResourceGridTextHandles
        % ResourceGridInfo Text information for ResourceGridTextHandles.
        % First element contains text related to downlink and second
        % element contains text related to uplink for FDD mode. In TDD mode
        % first element contains text related to both downlink and uplink.
        ResourceGridInfo = cell(2, 1)
        % RBGSize Number of RBs in an RBG. First element represents RBG
        % size for PDSCHRBs and second element represents RBG size for
        % PUSCHRBS
        RBGSize = zeros(2, 1);
        % RVCurrView Type of scheduler scheduling information displayed in
        % CQI Visualization. Value 1 represents downlink and value 2
        % represents uplink
        RVCurrView = 1
        % RGTxtHandle UI control handle to display the frame number in resource grid visualization
        RGTxtHandle
        % RGSlotTxtHandle UI control handle to display the slot number in resource grid visualization
        RGSlotTxtHandle
        % RGLowerRBIndex Index of the first RB displayed in resource grid visualization
        RGLowerRBIndex = 0
        % RGUpperRBIndex Index of the last RB displayed in resource grid visualization
        RGUpperRBIndex
        % RGLowerSlotIndex Index of the first slot displayed in resource grid visualization
        RGLowerSlotIndex = 0
        % RGUpperSlotIndex Index of the last slot displayed in resource grid visualization
        RGUpperSlotIndex

        % CQI information related properties
        % CQIInfo First element contains downlink CQI information and
        % second element contains uplink CQI information. Each element is
        % a N-by-P matrix where 'N' is the number of UEs and 'P' is the
        % number of RBs in the bandwidth. A matrix element at position (i,
        % j) corresponds to CQI value for UE with RNTI 'i' at RB 'j'
        CQIInfo = cell(2, 1);
        % CQIVisualizationGridHandles Handles to display UE CQIs on the RBs of the bandwidth
        CQIVisualizationGridHandles
        % CQIVisualizationFigHandle Handle of the CQI visualization
        CQIVisualizationFigHandle
        % CVCurrView Type of channel quality displayed in CQI
        % Visualization. Value 1 represents downlink and value 2 represents
        % uplink
        CVCurrView = 1
        % CVLowerUEIndex Index of the first UE to be displayed in CQI visualization
        CVLowerUEIndex = 0
        % CVUpperUEIndex Index of the last UE to be displayed in CQI visualization
        CVUpperUEIndex
        % CVLowerRBIndex Index of the first RB to be displayed in CQI visualization
        CVLowerRBIndex = 0
        % CVUpperRBIndex Index of the last RB to be displayed in CQI visualization
        CVUpperRBIndex
        % CVTxtHandle UI control handle to display the frame number in CQI visualization
        CVTxtHandle
        
        % LogInterval
        LogInterval

        % StepSize
        StepSize
        
        %UEMetricsUL UE metrics for each slot in the UL direction
        % It is an array of size N-by-3 where N is the number of UEs in
        % each cell. Each column of the array contains the following
        % metrics: throughput bytes transmitted, goodput bytes transmitted,
        % and pending buffer amount bytes.
        UEMetricsUL
        
        %UEMetricsDL UE metrics for each slot in the DL direction
        % It is an array of size N-by-3 where N is the number of UEs in
        % each cell. Each column of the array contains the following
        % metrics: throughput bytes transmitted, goodput bytes transmitted,
        % and pending buffer amount bytes.
        UEMetricsDL
        
        %UplinkChannelQuality Current channel quality for the UEs in uplink
        % It is an array of size M-by-N where M and N represents the number
        % of UEs in each cell and the number of RBs respectively.
        UplinkChannelQuality
        
        %DownlinkChannelQuality Current channel quality for the UEs in downlink
        % It is an array of size M-by-N where M and N represents the number
        % of UEs in each cell and the number of RBs respectively.
        DownlinkChannelQuality
        
        %HARQProcessStatusUL HARQ process status for each UE in UL
        % It is an array of size M-by-N where M and N represents the number
        % of UEs and number of HARQ processes for each UE respectively. Each
        % element stores the last received new data indicator (NDI) values 
        % in the uplink
        HARQProcessStatusUL
        
        %HARQProcessStatusDL HARQ process status for each UE in DL
        % It is an array of size M-by-N where M and N represents the number
        % of UEs and number of HARQ processes for each UE respectively. Each
        % element stores the last received new data indicator (NDI) values 
        % in the downlink
        HARQProcessStatusDL
        
        % PeakDLSpectralEfficiency Theoretical peak spectral efficiency in
        % the downlink direction
        PeakDLSpectralEfficiency
        
        % PeakULSpectralEfficiency Theoretical peak spectral efficiency in
        % the uplink direction
        PeakULSpectralEfficiency
    end

    methods
        function obj = hNRSchedulingLogger(simParameters,  varargin)
            %hNRSchedulingLogger Construct scheduling log and visualization object
            %
            % OBJ = hNRSchedulingLogger(SIMPARAMETERS) Create scheduling
            % information logging and visualization object. It creates
            % figures for visualizing both downlink and uplink information.
            %
            % OBJ = hNRSchedulingLogger(SIMPARAMETERS, FLAG) Create scheduling
            % information logging and visualization object.
            %
            % OBJ = hNRSchedulingLogger(SIMPARAMETERS, FLAG, ISLOGREPLAY)
            % Create scheduling information logging and visualization
            % object.
            %
            % SIMPARAMETERS - It is a structure and contain simulation
            % configuration information.
            %
            % NumFramesSim      - Simulation time in terms of number of 10 ms frames
            % NumUEs            - Number of UEs
            % NCellID           - Cell identifier
            % DuplexMode        - Duplexing mode (FDD or TDD)
            % SchedulingType    - Slot-based or symbol-based scheduling
            % NumHARQ           - Number of HARQ processes
            % NumRBs            - Number of resource blocks in PUSCH and PDSCH bandwidth
            % SCS               - Subcarrier spacing
            % DLBandwidth       - Downlink bandwidth
            % ULBandwidth       - Uplink bandwidth
            % DLULPeriodicity   - Duration of the DL-UL pattern in ms (for 
            %                     TDD mode)
            % NumDLSlots        - Number of full DL slots at the start of 
            %                     DL-UL pattern (for TDD mode)
            % NumDLSyms         - Number of DL symbols after full DL slots 
            %                     in the DL-UL pattern (for TDD mode)
            % NumULSlots        - Number of full UL slots at the end of 
            %                     DL-UL pattern (for TDD mode)
            % NumULSyms         - Number of UL symbols before full UL slots
            %                     in the DL-UL pattern (for TDD mode)
            % RBGSizeConfig     - Configuration table to use for
            %                     determining the RBG size (value 1
            %                     represents table-1 and value 2 represent
            %                     table-2)
            % TTIGranularity    - Minimum TTI granularity in terms of
            %                     number of symbols (for symbol-based scheduling)
            % CQIVisualization  - Flag to enable or disable CQI
            %                     visualization
            % RBVisualization   - Flag to enable or disable resource block
            %                     grid visualization
            %
            % If FLAG = 0, Visualize downlink information.
            % If FLAG = 1, Visualize uplink information.
            % If FLAG = 2, Visualize  downlink and uplink information.
            %
            % ISLOGREPLAY = true Replays the resource assignment of the simulation
            % frame-by-frame (or slot-by-slot).
            % ISLOGREPLAY = false Gives the option to select a particular frame
            % (or slot) to see the way resources are assigned in the chosen
            % frame (or slot).

            % Validate number of frames in simulation
            obj.NumFrames = simParameters.NumFramesSim;

            if isfield(simParameters , 'NCellID')
                obj.NCellID = simParameters.NCellID;
            end

            obj.NumUEs = simParameters.NumUEs;
            obj.UEIdList = 1:obj.NumUEs;

            if isfield(simParameters, 'NumHARQ')
                obj.NumHARQ = simParameters.NumHARQ;
            end
            if isfield(simParameters, 'SchedulingType')
                obj.SchedulingType = simParameters.SchedulingType;
            end
            obj.ColumnIndexMap = containers.Map('KeyType','char','ValueType','double');
            obj.GrantLogsColumnIndexMap = containers.Map('KeyType','char','ValueType','double');
            obj.NumSlotsFrame = (10 * simParameters.SCS) / 15; % Number of slots in a 10 ms frame
            % Validate the RB visualization flag
            if isfield(simParameters, 'RBVisualization')
                if islogical(simParameters.RBVisualization)
                    % To support true/false
                    validateattributes(simParameters.RBVisualization, {'logical'}, {'nonempty', 'scalar'}, 'simParameters.RBVisualization', 'RBVisualization');
                else
                    % To support 0/1
                    validateattributes(simParameters.RBVisualization, {'numeric'}, {'nonempty', 'integer', 'scalar', '>=', 0, '<=', 1}, 'simParameters.RBVisualization', 'RBVisualization');
                end
                obj.EnableResourceGridVisualization = simParameters.RBVisualization;
            end
            
            % Validate the CQI visualization
            if isfield(simParameters, 'CQIVisualization')
                if islogical(simParameters.CQIVisualization)
                    % To support true/false
                    validateattributes(simParameters.CQIVisualization, {'logical'}, {'nonempty', 'scalar'}, 'simParameters.CQIVisualization', 'CQIVisualization');
                else
                    % To support 0/1
                    validateattributes(simParameters.CQIVisualization, {'numeric'}, {'nonempty', 'integer', 'scalar', '>=', 0, '<=', 1}, 'simParameters.CQIVisualization', 'CQIVisualization');
                end
                obj.EnableCQIGridVisualization = simParameters.CQIVisualization;
            end

            % Symbol duration for the given numerology
            symbolDuration = 1e-3/(14*(simParameters.SCS/15)); % Assuming normal cyclic prefix
            
            % Validate the Tx antenna and Rx antenna count at gNB and UE
            % respectively
            if ~isfield(simParameters, 'GNBTxAnts')
                simParameters.GNBTxAnts = 1;
            elseif ~ismember(simParameters.GNBTxAnts, [1,2,4,8,16,32,64,128,256,512,1024])
                error('nr5g:hNRSchedulingLogger:InvalidAntennaSize',...
                    'Number of gNB Tx antennas (%d) must be a member of [1,2,4,8,16,32,64,128,256,512,1024].', simParameters.GNBTxAnts);                                
            end
            if ~isfield(simParameters, 'UERxAnts')
                simParameters.UERxAnts = 1;
            % Validate the number of receiver antennas on UEs
            elseif ~ismember(simParameters.UERxAnts, [1,2,4,8,16])              
                error('nr5g:hNRUEPhy:InvalidAntennaSize',...
                    'Number of UE Rx antennas (%d) must be a member of [1,2,4,8,16].', simParameters.UERxAnts(rnti));
            end

            % Maximum number of transmission layers for each UE in DL
            numLayersDL = min(simParameters.GNBTxAnts*ones(simParameters.NumUEs, 1), simParameters.UERxAnts);
            % Maximum number of transmission layers for each UE in UL (assuming SISO)
            numLayersUL = 1;
            % Verify Duplex mode and update the properties
            if isfield(simParameters, 'DuplexMode')
                obj.DuplexMode = simParameters.DuplexMode;
            end
            
            if isfield(simParameters, 'NumDLSlots')
                numDLSlots = simParameters.NumDLSlots;
            else
                numDLSlots = 2;
            end
            if isfield(simParameters, 'NumDLSyms')
                numDLSyms = simParameters.NumDLSyms;
            else
                numDLSyms = 8;
            end
            if isfield(simParameters, 'NumULSlots')
                numULSlots = simParameters.NumULSlots;
            else
                numULSlots = 2;
            end
            if isfield(simParameters, 'NumULSyms')
                numULSyms = simParameters.NumULSyms;
            else
                numULSyms = 4;
            end
            if isfield(simParameters, 'DLULPeriodicity')
                dlulPeriodicity = simParameters.DLULPeriodicity;
            else
                dlulPeriodicity = 5;
            end
            
            if obj.DuplexMode == obj.TDDDuplexMode % TDD
                obj.NumLogs = 1;
                obj.RVCurrView = 1; % Only one view for resource grid
                % Number of DL symbols in one DL-UL pattern
                numDLSymbols = numDLSlots*14 + numDLSyms;
                % Number of UL symbols in one DL-UL pattern
                numULSymbols = numULSlots*14 + numULSyms;
                % Number of symbols in one DL-UL pattern
                numSymbols = dlulPeriodicity*(simParameters.SCS/15)*14;
                % Normalized scalar considering the downlink symbol
                % allocation in the frame structure
                scaleFactorDL = numDLSymbols/numSymbols;
                % Normalized scalar considering the uplink symbol allocation
                % in the frame structure
                scaleFactorUL = numULSymbols/numSymbols;
            else % FDD
                obj.NumLogs = 2;
                % Normalized scalars in the DL and UL directions are 1 for
                % FDD mode
                scaleFactorDL = 1;
                scaleFactorUL = 1;
            end
            
            obj.UEMetricsUL = zeros(simParameters.NumUEs, 3);
            obj.UEMetricsDL = zeros(simParameters.NumUEs, 3);
            
            % Store current UL and DL CQI values on the RBs for the UEs.
            obj.UplinkChannelQuality = zeros(simParameters.NumUEs, simParameters.NumRBs);
            obj.DownlinkChannelQuality = zeros(simParameters.NumUEs, simParameters.NumRBs);
            
            % Store the last received new data indicator (NDI) values for UL and DL HARQ
            % processes.
            obj.HARQProcessStatusUL = zeros(simParameters.NumUEs, obj.NumHARQ);
            obj.HARQProcessStatusDL = zeros(simParameters.NumUEs, obj.NumHARQ);

            if isfield(simParameters, 'DLBandwidth')
                obj.Bandwidth(obj.DownlinkIdx) = simParameters.DLBandwidth;
            end
            if isfield(simParameters, 'ULBandwidth')
                obj.Bandwidth(obj.UplinkIdx) = simParameters.ULBandwidth;
            end
            % Calculate uplink and downlink peak data rates as per 3GPP TS
            % 37.910. The number of layers used for the peak DL data rate
            % calculation is taken as the average of maximum layers
            % possible for each UE. The maximum layers possible for each UE
            % is min(gNBTxAnts, ueRxAnts)
            % Determine the plots
            if isempty(varargin) || (nargin >= 2  && varargin{1} == 2)
                % Downlink & Uplink
                obj.PlotIds = [obj.DownlinkIdx obj.UplinkIdx];
                % Average of the peak DL throughput values for each UE
                obj.PeakDataRateDL = 1e-6*(sum(numLayersDL)/simParameters.NumUEs)*scaleFactorDL*8*(948/1024)*(simParameters.NumRBs*12)/symbolDuration;
                obj.PeakDataRateUL = 1e-6*scaleFactorUL*numLayersUL*8*(948/1024)*(simParameters.NumRBs*12)/symbolDuration;
                % Calculate uplink and downlink peak spectral efficiency
                obj.PeakDLSpectralEfficiency = 1e6*obj.PeakDataRateDL/obj.Bandwidth(obj.DownlinkIdx);
                obj.PeakULSpectralEfficiency = 1e6*obj.PeakDataRateUL/obj.Bandwidth(obj.UplinkIdx);
            elseif varargin{1} == 0 % Downlink
                obj.PlotIds = obj.DownlinkIdx;
                obj.PeakDataRateDL = 1e-6*(sum(numLayersDL)/simParameters.NumUEs)*scaleFactorDL*8*(948/1024)*(simParameters.NumRBs*12)/symbolDuration;
                % Calculate downlink peak spectral efficiency
                obj.PeakDLSpectralEfficiency = 1e6*obj.PeakDataRateDL/obj.Bandwidth(obj.DownlinkIdx);
            else % Uplink
                obj.PlotIds = obj.UplinkIdx;
                obj.PeakDataRateUL = 1e-6*scaleFactorUL*numLayersUL*8*(948/1024)*(simParameters.NumRBs*12)/symbolDuration;
                % Calculate uplink peak spectral efficiency
                obj.PeakULSpectralEfficiency = 1e6*obj.PeakDataRateUL/obj.Bandwidth(obj.UplinkIdx);
                obj.CVCurrView = obj.UplinkIdx;
                obj.RVCurrView = obj.UplinkIdx;
            end

            % Initialize number of RBs, RBG size, CQI and metrics properties
            if isfield(simParameters, 'RBGSizeConfig')
                obj.RBGSizeConfig = simParameters.RBGSizeConfig;
            end
            for idx = 1: numel(obj.PlotIds)
                logIdx = obj.PlotIds(idx);
                obj.NumRBs(logIdx) = simParameters.NumRBs; % Number of RBs in DL/UL
                % Calculate the RBGSize
                rbgSizeIndex = min(find(obj.NumRBs(logIdx) <= obj.NominalRBGSizePerBW(:, 1), 1));
                if obj.RBGSizeConfig == 1
                    obj.RBGSize(logIdx) = obj.NominalRBGSizePerBW(rbgSizeIndex, 2);
                else
                    obj.RBGSize(logIdx) = obj.NominalRBGSizePerBW(rbgSizeIndex, 3);
                end
                obj.CQIInfo{logIdx} = zeros(obj.NumUEs, obj.NumRBs(logIdx)); % DL/UL channel quality
            end

            if obj.SchedulingType % Symbol based scheduling
                gridLength = obj.NumSym;
            else % Slot based scheduling
                gridLength = obj.NumSlotsFrame;
            end

            % Initialize the scheduling logs and resources grid related
            % properties
            for idx=1:min(obj.NumLogs,numel(obj.PlotIds))
                plotId = obj.PlotIds(idx);
                if obj.DuplexMode == obj.FDDDuplexMode
                    logIdx = plotId; % FDD
                else
                    logIdx = idx; % TDD
                end
                % Construct the log format
                obj.SchedulingLog{logIdx} = constructLogFormat(obj, logIdx, simParameters);
                obj.ResourceGrid{logIdx} = zeros(gridLength, obj.NumRBs(plotId));
                obj.ResourceGridReTxInfo{logIdx} = zeros(gridLength, obj.NumRBs(plotId));
                obj.ResourceGridHarqInfo{logIdx} = zeros(gridLength, obj.NumRBs(plotId));
                obj.ResourceGridInfo{logIdx} = strings(gridLength, obj.NumRBs(plotId));
            end

            % Check if it is post simulation analysis
            if nargin == 3
                obj.IsLogReplay = varargin{2};
            end

            % Construct the grant log format
            obj.GrantLog = constructGrantLogFormat(obj, simParameters);

            % Create the visualization for cell of interest
            if ~isfield(simParameters, 'CellOfInterest') || obj.NCellID == simParameters.CellOfInterest

                % Using the screen width and height, calculate the figure width
                % and height
                resolution = get(0, 'ScreenSize');
                screenWidth = resolution(3);
                screenHeight = resolution(4);
                figureWidth = screenWidth * 0.90;
                figureHeight = screenHeight * 0.85;

                if(obj.EnableCQIGridVisualization) % Create CQI visualization
                    constructCQIGridVisualization(obj, screenWidth, screenHeight, figureWidth, figureHeight);
                end

                if(obj.EnableResourceGridVisualization) % Create resource grid visualization
                    constructResourceGridVisualization(obj, screenWidth, screenHeight, figureWidth, figureHeight);
                end
            end

            if ~isempty(obj.IsLogReplay) && obj.SchedulingType == obj.SlotBased
                % Post simulation log visualization and slot based scheduling
                obj.StepSize = 1;
                obj.LogInterval = 1;
            else
                % Live visualization
                obj.LogInterval = obj.NumSym;
                if obj.SchedulingType % Symbol based scheduling
                    obj.StepSize = 1;
                else % Slot based scheduling
                    obj.StepSize = obj.NumSym;
                end
            end
        end
        
        function [dlMetrics, ulMetrics, cellMetrics] = getMACMetrics(obj, firstSlot, lastSlot, rntiList)
            %getMACMetrics Returns the MAC metrics
            %
            % [DLMETRICS, ULMETRICS] = getMACMetrics(OBJ, FIRSTSLOT,
            % LASTSLOT, RNTILIST) Returns the MAC metrics of the UE with
            % specified RNTI within the cell for both uplink and downlink direction
            %
            % FIRSTSLOT - Represents the starting slot number for
            % querying the metrics
            %
            % LASTSLOT -  Represents the ending slot for querying the metrics
            % 
            % RNTILIST - Radio network temporary identifiers of the UEs
            %
            % ULMETRICS and DLMETRICS are array of structures with following properties
            %
            %   RNTI - Radio network temporary identifier of the UE
            %
            %   TxBytes - Total number of bytes transmitted (newTx and reTx combined)
            %    
            %   NewTxBytes - Number of bytes transmitted (only newTx)
            %
            %   BufferStatus - Current buffer status of the UE
            %
            %   AssignedRBCount - Number of resource blocks assigned to the UE
            %
            %   RBsScheduled - Total number resource blocks scheduled
            %
            % CELLMETRICS is an array structure with following properties and
            % contains cell wide metrics in downlink and uplink
            % respectively
            %
            %   DLTxBytes - Total number of bytes transmitted (newTx and
            %   reTx combined) in downlink
            %    
            %   DLNewTxBytes - Number of bytes transmitted (only newTx) in
            %   downlink
            %
            %   DLRBsScheduled - Total number resource blocks scheduled in
            %   downlink
            %
            %   ULTxBytes - Total number of bytes transmitted (newTx and
            %   reTx combined) in uplink
            %    
            %   ULNewTxBytes - Number of bytes transmitted (only newTx) in uplink
            %
            %   ULRBsScheduled - Total number resource blocks scheduled in uplink

            % Calculate the actual log start and end index
            stepLogStartIdx = (firstSlot-1) * obj.LogInterval + 1;
            stepLogEndIdx = lastSlot*obj.LogInterval;
            
            % Create structure for both DL and UL
            outStruct = struct('RNTI', 0, 'TxBytes', 0, ...
                'NewTxBytes', 0, 'BufferStatus', 0, ...
                'AssignedRBCount', 0, 'RBsScheduled', 0);
            outputStruct = repmat(outStruct, [numel(rntiList) 2]);
            assignedRBsStep = zeros(obj.NumUEs, 2);
            macTxStep = zeros(obj.NumUEs, 2);
            macNewTxStep = zeros(obj.NumUEs, 2);
            bufferStatus = zeros(obj.NumUEs, 2);
            
            % Update the DL and UL metrics properties 
            for idx = 1:min(obj.NumLogs, numel(obj.PlotIds))
                plotId = obj.PlotIds(idx);
                % Determine scheduling log index
                if obj.DuplexMode == obj.FDDDuplexMode
                    schedLogIdx = plotId;
                else
                    schedLogIdx = 1;
                end
                
                numULSyms = 0;
                numDLSyms = 0;
                
                % Read the information of each slot and update the metrics
                % properties
                for i = stepLogStartIdx:obj.StepSize:stepLogEndIdx
                    slotLog = obj.SchedulingLog{schedLogIdx}(i, :);
                    rbgAssignment = slotLog{obj.ColumnIndexMap('RBG Allocation Map')};
                    throughputBytes =  slotLog{obj.ColumnIndexMap('Throughput Bytes')};
                    goodputBytes =  slotLog{obj.ColumnIndexMap('Goodput Bytes')};
                    ueBufferStatus = slotLog{obj.ColumnIndexMap('Buffer status of UEs (In bytes)')};
                    if(obj.DuplexMode == obj.TDDDuplexMode)
                        switch (slotLog{obj.ColumnIndexMap('Type')})
                            case 'UL'
                                linkIdx = 2; % Uplink information index
                                numULSyms = numULSyms + 1;
                            case 'DL'
                                linkIdx = 1; % Downlink information index
                                numDLSyms = numDLSyms + 1;
                            otherwise
                                continue;
                        end
                    else
                        linkIdx = plotId;
                    end
                    
                    % Calculate the RBs allocated to an UE
                    for ueIdx = 1 : obj.NumUEs
                        numRBGs = sum(rbgAssignment(ueIdx, :));
                        if rbgAssignment(ueIdx, end) % If RBG is allocated
                            % If the last RBG of BWP is assigned, then it might not
                            % have same number of RBs as other RBG.
                            if(mod(obj.NumRBs(plotId), obj.RBGSize(plotId)) == 0)
                                numRBs = numRBGs * obj.RBGSize(plotId);
                            else
                                lastRBGSize = mod(obj.NumRBs(plotId), obj.RBGSize(plotId));
                                numRBs = (numRBGs - 1) * obj.RBGSize(plotId) + lastRBGSize;
                            end
                        else
                            numRBs = numRBGs * obj.RBGSize(plotId);
                        end
                        
                        assignedRBsStep(ueIdx, linkIdx) = assignedRBsStep(ueIdx, linkIdx) + numRBs;
                        macTxStep(ueIdx, linkIdx) = macTxStep(ueIdx, linkIdx) + throughputBytes(ueIdx);
                        macNewTxStep(ueIdx, linkIdx) = macNewTxStep(ueIdx, linkIdx) + goodputBytes(ueIdx);
                        bufferStatus(ueIdx, linkIdx) = ueBufferStatus(ueIdx);
                    end
                end
            end
            
            % Extract required metrics of the UEs specified in rntiList
            for idx = 1:numel(obj.PlotIds)
                linkIdx = obj.PlotIds(idx);
                for listIdx = 1:numel(rntiList)
                    ueIdx = find(rntiList(listIdx) == obj.UEIdList);
                    outputStruct(ueIdx, linkIdx).RNTI = rntiList(listIdx);
                    outputStruct(ueIdx, linkIdx).AssignedRBCount = assignedRBsStep(ueIdx, linkIdx);
                    outputStruct(ueIdx, linkIdx).TxBytes = macTxStep(ueIdx, linkIdx);
                    outputStruct(ueIdx, linkIdx).NewTxBytes = macNewTxStep(ueIdx, linkIdx);
                    outputStruct(ueIdx, linkIdx).BufferStatus = bufferStatus(ueIdx, linkIdx);
                end
            end

            dlMetrics =  outputStruct(:, obj.DownlinkIdx); % Downlink Info
            ulMetrics = outputStruct(:, obj.UplinkIdx); % Uplink Info
            % Cell wide metrics
            cellMetrics.DLTxBytes = sum(macTxStep(:, obj.DownlinkIdx));
            cellMetrics.DLNewTxBytes = sum(macNewTxStep(:, obj.DownlinkIdx));
            cellMetrics.ULTxBytes = sum(macTxStep(:, obj.UplinkIdx));
            cellMetrics.ULNewTxBytes = sum(macNewTxStep(:, obj.UplinkIdx));
            if (obj.DuplexMode == obj.TDDDuplexMode)
                RBsScheduledLastStep = obj.NumRBs(obj.DownlinkIdx) * ((stepLogEndIdx - stepLogStartIdx + 1)/obj.StepSize);
                cellMetrics.ULRBsScheduled = (numULSyms/((stepLogEndIdx - stepLogStartIdx +1)/obj.StepSize))*RBsScheduledLastStep;
                cellMetrics.DLRBsScheduled = (numDLSyms/((stepLogEndIdx - stepLogStartIdx +1)/obj.StepSize))*RBsScheduledLastStep;
            else
                cellMetrics.ULRBsScheduled = obj.NumRBs(obj.UplinkIdx) * ((stepLogEndIdx - stepLogStartIdx + 1)/obj.StepSize);
                cellMetrics.DLRBsScheduled =  obj.NumRBs(obj.DownlinkIdx) * ((stepLogEndIdx - stepLogStartIdx + 1)/obj.StepSize);
            end
        end

        function plotRBGrids(obj)
            %plotRBGrids Read the most recent logs to update resource grid
            %visualization
            % plotRBGrids(OBJ) Update the resource grid

            % Check if the figure handle is valid
            if isempty(obj.RGVisualizationFigHandle) || ~ishghandle(obj.RGVisualizationFigHandle)
                return;
            end

            obj.RGTxtHandle.String = ['Frame Number : ' num2str(obj.CurrFrame)];
            if obj.SchedulingType % Symbol based scheduling
                obj.RGSlotTxtHandle.String = ['Slot Number : ' num2str(obj.CurrSlot)];
                frameLogStartIdx = (obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval) + (obj.CurrSlot * obj.LogInterval);
                frameLogEndIdx = frameLogStartIdx + obj.LogInterval;
            else % Slot based scheduling
                frameLogStartIdx = obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval;
                frameLogEndIdx = frameLogStartIdx + (obj.NumSlotsFrame * obj.LogInterval);
            end

            for idx = 1:min(obj.NumLogs, numel(obj.PlotIds))
                plotId = obj.PlotIds(idx);
                if obj.DuplexMode == obj.FDDDuplexMode
                    logIdx = obj.PlotIds(idx);
                else
                    logIdx = 1;
                end

                % Reset the resource grid status
                if obj.SchedulingType % Symbol based scheduling
                    obj.ResourceGrid{logIdx} = zeros(length(obj.SymbolInfo), obj.NumRBs(logIdx));
                else % Slot based scheduling
                    obj.ResourceGrid{logIdx} = zeros(obj.NumSlotsFrame, obj.NumRBs(logIdx));
                end

                slIdx = 0; % Counter to keep track of the number of symbols/slots to be plotted
                for i = frameLogStartIdx+1:obj.StepSize:frameLogEndIdx % For each symbol in the slot or each slot in the frame
                    slIdx = slIdx + 1;
                    slotLog = obj.SchedulingLog{logIdx}(i, :);
                    rbgAssignment = slotLog{obj.ColumnIndexMap('RBG Allocation Map')};
                    harqIds = slotLog{obj.ColumnIndexMap('HARQ Process ID')};
                    txType = slotLog{obj.ColumnIndexMap('Transmission')};
                    if obj.SchedulingType % Symbol based scheduling
                        if obj.DuplexMode == obj.TDDDuplexMode
                            obj.SymbolInfo{slIdx} = slotLog{obj.ColumnIndexMap('Type')};
                        end
                        % Plot the selected visualization of interest
                        if numel(obj.PlotIds) == 1
                            if obj.PlotIds == obj.UplinkIdx && strcmp(obj.SymbolInfo{slIdx}, 'DL') % Skip downlink symbols
                                continue;
                            elseif obj.PlotIds == obj.DownlinkIdx && strcmp(obj.SymbolInfo{slIdx}, 'UL') % Skip uplink symbols
                                continue;
                            end
                        end
                    else
                        if obj.DuplexMode == obj.TDDDuplexMode
                            obj.SlotInfo{slIdx} = slotLog{obj.ColumnIndexMap('Type')};
                            % Plot the selected visualization of interest
                            if numel(obj.PlotIds) == 1
                                if obj.PlotIds == obj.UplinkIdx && strcmp(obj.SlotInfo{slIdx}, 'DL') % Skip downlink slots
                                    continue;
                                elseif obj.PlotIds == obj.DownlinkIdx && strcmp(obj.SlotInfo{slIdx}, 'UL') % Skip uplink slots
                                    continue;
                                end
                            end
                        end
                    end
                    for j = 1 : obj.NumUEs % For each UE
                        if (strcmp(txType(j), 'newTx') || strcmp(txType(j), 'newTx-Start') || strcmp(txType(j), 'newTx-InProgress') || strcmp(txType(j), 'newTx-End'))
                            type = 1; % New transmission
                        else
                            type = 2; % Retransmission
                        end

                        % Updating the resource grid status and related
                        % information
                        RBGAllocationBitmap = rbgAssignment(j, :);
                        for k=1:length(RBGAllocationBitmap) % For all RBGs
                            if(RBGAllocationBitmap(k) == 1)
                                startRBIndex = (k - 1) * obj.RBGSize(plotId) + 1;
                                endRBIndex =  k * obj.RBGSize(plotId);
                                if(k == length(RBGAllocationBitmap) && (mod(obj.NumRBs(plotId), obj.RBGSize(plotId)) ~=0))
                                    % If it is last RBG and it does not
                                    % have same number of RBs as other RBGs
                                    endRBIndex = (k - 1) * obj.RBGSize(plotId) + mod(obj.NumRBs(plotId), obj.RBGSize(plotId));
                                end
                                obj.ResourceGrid{logIdx}(slIdx, startRBIndex : endRBIndex) = j;
                                obj.ResourceGridReTxInfo{logIdx}(slIdx, startRBIndex : endRBIndex) = type;
                                obj.ResourceGridHarqInfo{logIdx}(slIdx, startRBIndex : endRBIndex) = harqIds(j);
                            end
                        end
                    end
                end

                for p = 1:slIdx
                    for q = 1 : obj.NumRBs(plotId)
                        if(obj.ResourceGrid{logIdx}(p, q) == 0)
                            % Clear the previously plotted text in the resource grid
                            obj.ResourceGridInfo{logIdx}(p, q)  = '';
                        else
                            % Create the text to be plotted in the resource
                            % grid
                            obj.ResourceGridInfo{logIdx}(p, q) = strcat('UE-', num2str(obj.ResourceGrid{logIdx}(p, q)), ...
                                '( ', num2str(obj.ResourceGridHarqInfo{logIdx}(p, q)), ')');
                        end
                    end
                end
            end

            % Update the resource grid visualization
            updateResourceGridVisualization(obj);
            drawnow;
        end

        function plotCQIRBGrids(obj)
            %plotCQIRBGrids Update channel quality visualization
            % plotCQIRBGrids(OBJ) Update the channel quality visualization

            % Check if the figure handle is valid
            if isempty(obj.CQIVisualizationFigHandle) || ~ishghandle(obj.CQIVisualizationFigHandle)
                return;
            end

            % Update frame number in the figure
            obj.CVTxtHandle.String = ['Frame Number : '  num2str(obj.CurrFrame)];
            if obj.SchedulingType % Symbol based scheduling
                lwRowIndex = (obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval) + 1;
                upRowIndex = (obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval) + (obj.CurrSlot + 1) * obj.LogInterval;
            else % Slot based scheduling
                lwRowIndex = (obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval) + 1;
                upRowIndex = (obj.CurrFrame * obj.NumSlotsFrame * obj.LogInterval) + (obj.CurrSlot * obj.LogInterval) + 1;
           end

            if (obj.DuplexMode == obj.TDDDuplexMode) % TDD
                % Get the symbols types in the current frame
                symbolTypeInFrame = {obj.SchedulingLog{1}(lwRowIndex:upRowIndex, obj.ColumnIndexMap('Type'))};
                % Get the UL symbol indices
                ulIdx = find(strcmp(symbolTypeInFrame{1}, 'UL'));
                % Get the DL symbol indices
                dlIdx = find(strcmp(symbolTypeInFrame{1}, 'DL'));
                if ~isempty(dlIdx)
                    % Update downlink channel quality based on latest DL
                    % symbol/slot
                    obj.CQIInfo{obj.DownlinkIdx} = obj.SchedulingLog{1}{lwRowIndex + dlIdx(end) - 1, obj.ColumnIndexMap('Channel Quality')};
                end
                if ~isempty(ulIdx)
                    % Update uplink channel quality based on latest UL
                    % symbol/slot
                    obj.CQIInfo{obj.UplinkIdx} = obj.SchedulingLog{1}{lwRowIndex + ulIdx(end) - 1, obj.ColumnIndexMap('Channel Quality')};
                end
            else
                for idx=1:numel(obj.PlotIds)
                    plotId = obj.PlotIds(idx);
                    obj.CQIInfo{plotId} = obj.SchedulingLog{plotId}{upRowIndex, obj.ColumnIndexMap('Channel Quality')};
                end
            end

            % Update the CQI grid visualization
            updateCQIVisualization(obj);
            drawnow;
        end
        
        function logCellSchedulingStats(obj, symbolNum, gNB, UEs, varargin)
            %logCellSchedulingStats Log the MAC layer statistics
            %
            % LOGCELLSCHEDULINGSTATS(OBJ, SYMBOLNUM, GNB, UES, LINKDIR) Logs 
            % the scheduling stats for all the nodes in the cell
            % 
            % SYMBOLNUM - Symbol number in the simulation
            % 
            % GNB - It is an object of type hNRGNB and contains information
            % about the gNB
            %
            % UEs - It is a cell array of length equal to the number of UEs
            % in the cell. Each element of the array is an object of type
            % hNRUE.
            %
            % LINKDIR - Indicates the downlink/uplink direction. 0 and 1
            % denotes downlink and uplink, respectively.
            
            if ~isempty(varargin)
                linkDir = varargin{1};
            else
                linkDir = 2; % For both UL & DL
            end

            % Read UL and DL assignments by gNB MAC scheduler
            % at current time. Resource assignments returned by a scheduler (either
            % UL or DL) are empty if either the scheduler was not scheduled to run at
            % the current time or the scheduler did not schedule any resources
            [resourceAssignmentsUL, resourceAssignmentsDL] = getCurrentSchedulingAssignments(gNB.MACEntity);
            % Read throughput and goodput bytes sent for each UE
            [obj.UEMetricsDL(:, 1), obj.UEMetricsDL(:, 2)] = getTTIBytes(gNB);
            obj.UEMetricsDL(:, 3) = getBufferStatus(gNB); % Read pending buffer (in bytes) on gNB, for all the UEs
            for ueIdx = 1:obj.NumUEs
                obj.HARQProcessStatusUL(ueIdx, :) = getLastNDIFlagHarq(UEs{ueIdx}.MACEntity, 1); % 1 for UL
                obj.HARQProcessStatusDL(ueIdx, :) = getLastNDIFlagHarq(UEs{ueIdx}.MACEntity, 0); % 0 for DL
                % Read the UL channel quality at gNB for each of the UEs for logging
                obj.UplinkChannelQuality(ueIdx,:) = getChannelQuality(gNB, 1, ueIdx); % 1 for UL
                % Read the DL channel quality at gNB for each of the UEs for logging
                obj.DownlinkChannelQuality(ueIdx,:) = getChannelQuality(gNB, 0, ueIdx); % 0 for DL
                % Read throughput and goodput bytes transmitted for the UE in the
                % current TTI for logging
                [obj.UEMetricsUL(ueIdx, 1), obj.UEMetricsUL(ueIdx, 2)] = getTTIBytes(UEs{ueIdx});
                obj.UEMetricsUL(ueIdx, 3) = getBufferStatus(UEs{ueIdx}); % Read pending buffer (in bytes) on UE
            end
            if obj.DuplexMode == 1 % TDD
                symbolType = currentSymbolType(gNB); % Get current symbol type: DL/UL/Guard
                if(symbolType == 0 && linkDir ~= 1) %DL
                    logScheduling(obj, symbolNum, [resourceAssignmentsUL resourceAssignmentsDL], obj.UEMetricsDL, obj.DownlinkChannelQuality, obj.HARQProcessStatusDL, symbolType);
                elseif(symbolType == 1 && linkDir ~= 0) %UL
                    logScheduling(obj, symbolNum, [resourceAssignmentsUL resourceAssignmentsDL], obj.UEMetricsUL, obj.UplinkChannelQuality, obj.HARQProcessStatusUL, symbolType);
                else % Guard
                    logScheduling(obj, symbolNum, [resourceAssignmentsUL resourceAssignmentsDL], zeros(obj.NumUEs, 3), zeros(obj.NumUEs, obj.NumRBs(1)), zeros(obj.NumUEs, 16), symbolType); % UL
                end
            else
                % Store the scheduling logs
                if linkDir ~= 1 %  DL
                    logScheduling(obj, symbolNum, resourceAssignmentsDL, obj.UEMetricsDL, obj.DownlinkChannelQuality, obj.HARQProcessStatusDL, 0); % DL
                end
                if linkDir ~= 0 % UL
                    logScheduling(obj, symbolNum, resourceAssignmentsUL, obj.UEMetricsUL, obj.UplinkChannelQuality, obj.HARQProcessStatusUL, 1); % UL
                end
            end           
        end
        
        function logScheduling(obj, symbolNumSimulation, resourceAssignments, UEMetrics, UECQIs, HarqProcessStatus, type)
            %logScheduling Log the scheduling operations
            %
            % logScheduling(OBJ, SYMBOLNUMSIMULATION, RESOURCEASSIGNMENTS,
            % UEMETRICS, UECQIS, HARQPROCESSSTATUS, RXRESULTUES, TYPE) Logs
            % the scheduling operations based on the input arguments
            %
            % SYMBOLNUMSIMULATION - Cumulative symbol number in the
            % simulation
            %
            % RESOURCEASSIGNMENTS - Resource assignment information.
            %
            % UEMETRICS - N-by-P matrix where N represents the number of
            % UEs and P represents the number of metrics collected.
            %
            % UECQIs - N-by-P matrix where N represents the number of
            % UEs and P represents the number of RBs.
            %
            % HARQPROCESSSTATUS - N-by-P matrix where N represents the number of
            % UEs and P represents the number of HARQ process.
            %
            % TYPE - Type will be based on scheduling type.
            %        - In slot based scheduling type takes two values.
            %          type = 0, represents the downlink and type = 1,
            %          represents uplink.
            %
            %        - In symbol based scheduling type takes three values.
            %          type = 0, represents the downlink, type = 1,
            %          represents uplink and type = 2 represents guard.

            % Determine the log index based on link type and duplex mode
            if obj.DuplexMode == obj.FDDDuplexMode
                if  type == 0
                    linkIndex = obj.DownlinkIdx; % Downlink log
                else
                    linkIndex = obj.UplinkIdx; % Uplink log
                end
            else
                % TDD
                linkIndex = 1;
            end

            % Calculate symbol number in slot (0-13), slot number in frame
            % (0-obj.NumSlotsFrame), and frame number in the simulation.
            slotDuration = 10/obj.NumSlotsFrame;
            obj.CurrSymbol = mod(symbolNumSimulation - 1, obj.NumSym);
            obj.CurrSlot = mod(floor((symbolNumSimulation - 1)/obj.NumSym), obj.NumSlotsFrame);
            obj.CurrFrame = floor((symbolNumSimulation-1)/(obj.NumSym * obj.NumSlotsFrame));
            timestamp = obj.CurrFrame * 10 + (obj.CurrSlot * slotDuration) + (obj.CurrSymbol * (slotDuration / 14));

            columnMap = obj.ColumnIndexMap;
            grantLogsColumnIndexMap = obj.GrantLogsColumnIndexMap;
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, columnMap('Timestamp')} = timestamp;
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, columnMap('Frame Number')} = obj.CurrFrame;
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, columnMap('Slot Number')} = obj.CurrSlot;
            if obj.SchedulingType % Symbol based scheduling
                obj.SchedulingLog{linkIndex}{symbolNumSimulation, columnMap('Symbol Number')} = obj.CurrSymbol;
            end

            if(obj.DuplexMode == obj.TDDDuplexMode) % TDD
                % Log the type: DL/UL/Guard
                switch(type)
                case 0
                    symbolTypeDesc = 'DL';
                case 1
                    symbolTypeDesc = 'UL';
                case 2
                    symbolTypeDesc = 'Guard';
                end
                obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('Type')} = symbolTypeDesc;
            end

            for j = 1:length(resourceAssignments)
                % Fill logs w.r.t. each assignment
                assignment = resourceAssignments{j};
                % Calculate row number in the logs, for the Tx start
                % symbol
                logIndex = (obj.CurrFrame * obj.NumSlotsFrame * obj.NumSym) +  ...
                    ((obj.CurrSlot + assignment.SlotOffset) * obj.NumSym) + assignment.StartSymbol + 1;

                allottedUE = assignment.RNTI;

                % Fill the start Tx symbol logs
                obj.SchedulingLog{linkIndex}{logIndex, columnMap('RBG Allocation Map')}(allottedUE, :) = assignment.RBGAllocationBitmap;
                obj.SchedulingLog{linkIndex}{logIndex, columnMap('UEs MCS')}(allottedUE) = assignment.MCS;
                obj.SchedulingLog{linkIndex}{logIndex, columnMap('HARQ Process ID')}(allottedUE) = assignment.HARQID;
                obj.SchedulingLog{linkIndex}{logIndex, columnMap('Grant NDI Flag')}(allottedUE) = assignment.NDI;
                if obj.SchedulingType % Symbol based scheduling
                    obj.SchedulingLog{linkIndex}{logIndex, columnMap('Transmission')}(allottedUE) = {strcat(assignment.Type, '-Start')};
                    % Fill the logs from the symbol after Tx start, up to
                    % the symbol before Tx end
                    for k = 1:assignment.NumSymbols-2
                        obj.SchedulingLog{linkIndex}{logIndex + k, columnMap('RBG Allocation Map')}(allottedUE, :) = assignment.RBGAllocationBitmap;
                        obj.SchedulingLog{linkIndex}{logIndex + k, columnMap('UEs MCS')}(allottedUE) = assignment.MCS;
                        obj.SchedulingLog{linkIndex}{logIndex + k, columnMap('HARQ Process ID')}(allottedUE) = assignment.HARQID;
                        obj.SchedulingLog{linkIndex}{logIndex + k, columnMap('Grant NDI Flag')}(allottedUE) = assignment.NDI;
                        obj.SchedulingLog{linkIndex}{logIndex + k, columnMap('Transmission')}(allottedUE) = {strcat(assignment.Type, '-InProgress')};
                    end

                    % Fill the last Tx symbol logs
                    obj.SchedulingLog{linkIndex}{logIndex + assignment.NumSymbols -1, columnMap('RBG Allocation Map')}(allottedUE, :) = assignment.RBGAllocationBitmap;
                    obj.SchedulingLog{linkIndex}{logIndex + assignment.NumSymbols -1, columnMap('UEs MCS')}(allottedUE) = assignment.MCS;
                    obj.SchedulingLog{linkIndex}{logIndex + assignment.NumSymbols -1, columnMap('HARQ Process ID')}(allottedUE) = assignment.HARQID;
                    obj.SchedulingLog{linkIndex}{logIndex + assignment.NumSymbols -1, columnMap('Grant NDI Flag')}(allottedUE) = assignment.NDI;
                    obj.SchedulingLog{linkIndex}{logIndex + assignment.NumSymbols -1, columnMap('Transmission')}(allottedUE) = {strcat(assignment.Type, '-End')};
                else % Slot based scheduling
                    obj.SchedulingLog{linkIndex}{logIndex, columnMap('Transmission')}(allottedUE) = {assignment.Type};
                end
                obj.GrantCount  = obj.GrantCount + 1;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('RNTI')} = assignment.RNTI;
                slotNumGrant = mod(obj.CurrSlot + assignment.SlotOffset, obj.NumSlotsFrame);
                if(obj.CurrSlot + assignment.SlotOffset >= obj.NumSlotsFrame)
                    frameNumGrant = obj.CurrFrame + 1; % Assignment is for a slot in next frame
                else
                    frameNumGrant = obj.CurrFrame;
                end
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Frame')} = frameNumGrant;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Slot')} = slotNumGrant;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('RBG Allocation Map')} = mat2str(assignment.RBGAllocationBitmap);
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Start Sym')} = assignment.StartSymbol;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Num Sym')} = assignment.NumSymbols;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('MCS')} = assignment.MCS;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('NumLayers')} = assignment.NumLayers;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('HARQ ID')} = assignment.HARQID;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('NDI Flag')} = assignment.NDI;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('RV')} = assignment.RV;
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Tx Type')} = assignment.Type;
                if(isfield(assignment, 'FeedbackSlotOffset'))
                    % DL grant
                    obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Feedback Slot Offset (DL grants only)')} = assignment.FeedbackSlotOffset;
                    obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Grant type')} = 'DL';
                else
                    % UL Grant
                    obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('Grant type')} = 'UL';
                end
                obj.GrantLog{obj.GrantCount, grantLogsColumnIndexMap('CQI on RBs')} = mat2str(UECQIs(assignment.RNTI, :));
            end

            obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('Channel Quality')} = UECQIs;
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('HARQ process NDI status (at symbol start)')} = HarqProcessStatus;
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('Throughput Bytes')} = UEMetrics(:, 1); % Throughput bytes sent by UEs
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('Goodput Bytes')} = UEMetrics(:, 2); % Goodput bytes sent by UEs
            obj.SchedulingLog{linkIndex}{symbolNumSimulation, obj.ColumnIndexMap('Buffer status of UEs (In bytes)')} = UEMetrics(:, 3); % Current buffer status of UEs in bytes
        end

        function varargout = getSchedulingLogs(obj)
            %getSchedulingLogs Get the per-symbol logs of the whole simulation

            % Get keys of columns (i.e. column names) in sorted order of values (i.e. column indices)
            [~, idx] = sort(cell2mat(values(obj.ColumnIndexMap)));
            columnTitles = keys(obj.ColumnIndexMap);
            columnTitles = columnTitles(idx);
            varargout = cell(obj.NumLogs, 1);

            for logIdx = 1:obj.NumLogs
                if isempty(obj.SchedulingLog{logIdx})
                    continue;
                end
                if obj.SchedulingType
                    % Symbol based scheduling
                    finalLogIndex = (obj.CurrFrame)*obj.NumSlotsFrame*obj.NumSym + (obj.CurrSlot)*obj.NumSym + obj.CurrSymbol + 1;
                    obj.SchedulingLog{logIdx} = obj.SchedulingLog{logIdx}(1:finalLogIndex, :);
                    % For symbol based scheduling, keep 1 row per symbol
                    varargout{logIdx} = [columnTitles; obj.SchedulingLog{logIdx}(1:finalLogIndex, :)];
                else
                    % Slot based scheduling
                    finalLogIndex = (obj.CurrFrame)*obj.NumSlotsFrame*obj.NumSym + (obj.CurrSlot+1)*obj.NumSym;
                    obj.SchedulingLog{logIdx} = obj.SchedulingLog{logIdx}(1:finalLogIndex, :);
                    % For slot based scheduling: keep 1 row per slot and eliminate symbol number as a column title
                    varargout{logIdx} = [columnTitles; obj.SchedulingLog{logIdx}(1:obj.NumSym:finalLogIndex, :)];
                end
            end
        end

        function logs = getGrantLogs(obj)
            %getGrantLogs Get the scheduling assignment logs of the whole simulation

            % Get keys of columns (i.e. column names) in sorted order of values (i.e. column indices)
            [~, idx] = sort(cell2mat(values(obj.GrantLogsColumnIndexMap)));
            columnTitles = keys(obj.GrantLogsColumnIndexMap);
            columnTitles = columnTitles(idx);
            % Read valid rows
            obj.GrantLog = obj.GrantLog(1:obj.GrantCount, :);
            logs = [columnTitles; obj.GrantLog];
        end

        function plotPostSimRBGrids(obj, simSlotNum)
            %plotPostSimRBGrids Post simulation log visualization
            %
            % plotPostSimRBGrids(OBJ, SIMSLOTNUM) To update the resource
            % grid and CQI visualization based on the post simulation logs.
            %
            % SIMSLOTNUM - Cumulative slot number in the simulation

            % Update slot number
            if obj.SchedulingType % Symbol based scheduling
                obj.CurrSlot = mod(simSlotNum-1, obj.NumSlotsFrame);
                if obj.CurrSlot == 0
                    obj.CurrFrame = floor(simSlotNum/obj.NumSlotsFrame);
                end
            else % Slot based scheduling
                obj.CurrSlot = obj.NumSlotsFrame - 1;
                obj.CurrFrame = floor(simSlotNum/obj.NumSlotsFrame) - 1;
            end

            % Update grid information at slot boundary (for symbol based
            % scheduling) and frame boundary (for slot based scheduling)
            % Update resource grid visualization
            plotRBGrids(obj);
            % Update CQI visualization
            plotCQIRBGrids(obj);
        end
        
        function [dlStats, ulStats] = getPerformanceIndicators(obj)
            %getPerformanceIndicators Outputs the data rate, spectral
            % efficiency values 
            %
            % DLSTATS - 5-by-1 array containing the following statistics in
            %           the downlink direction: Theoretical peak data rate,
            %           achieved data rate, theoretical peak spectral
            %           efficiency, achieved spectral efficiency, achieved
            %           goodput
            % ULSTATS - 5-by-1 array containing the following statistics in
            %           the uplink direction: Theoretical peak data rate,
            %           achieved data rate, theoretical peak spectral
            %           efficiency, achieved spectral efficiency, achieved
            %           goodput

            if obj.DuplexMode == obj.FDDDuplexMode
                if ismember(obj.DownlinkIdx, obj.PlotIds)
                    totalDLTxBytes = sum(cell2mat(obj.SchedulingLog{obj.DownlinkIdx}(:,  obj.ColumnIndexMap('Throughput Bytes'))));
                    totalDLNewTxBytes = sum(cell2mat(obj.SchedulingLog{obj.DownlinkIdx}(:,  obj.ColumnIndexMap('Goodput Bytes'))));
                end
                if ismember(obj.UplinkIdx, obj.PlotIds)
                    totalULTxBytes = sum(cell2mat(obj.SchedulingLog{obj.UplinkIdx}(:,  obj.ColumnIndexMap('Throughput Bytes'))));
                    totalULNewTxBytes = sum(cell2mat(obj.SchedulingLog{obj.UplinkIdx}(:,  obj.ColumnIndexMap('Goodput Bytes'))));
                end
            else
                dlIdx = strcmp(obj.SchedulingLog{1}(:, obj.ColumnIndexMap('Type')), 'DL');
                totalDLTxBytes = sum(cell2mat(obj.SchedulingLog{1}(dlIdx,  obj.ColumnIndexMap('Throughput Bytes'))));
                totalDLNewTxBytes = sum(cell2mat(obj.SchedulingLog{1}(dlIdx,  obj.ColumnIndexMap('Goodput Bytes'))));
                ulIdx = strcmp(obj.SchedulingLog{1}(:, obj.ColumnIndexMap('Type')), 'UL');
                totalULTxBytes = sum(cell2mat(obj.SchedulingLog{1}(ulIdx,  obj.ColumnIndexMap('Throughput Bytes'))));
                totalULNewTxBytes = sum(cell2mat(obj.SchedulingLog{1}(ulIdx,  obj.ColumnIndexMap('Goodput Bytes'))));
            end
            dlStats = zeros(5, 1);
            ulStats = zeros(5, 1);

            % Downlink stats
            if ismember(obj.DownlinkIdx, obj.PlotIds)
                dlStats(1, 1) = obj.PeakDataRateDL;
                dlStats(2, 1) = totalDLTxBytes * 8 ./ (obj.NumFrames * 0.01 * 1000 * 1000); % Mbps
                dlStats(3, 1) = obj.PeakDLSpectralEfficiency;
                dlStats(4, 1) = 1e6*dlStats(2, 1)/obj.Bandwidth(obj.DownlinkIdx);
                dlStats(5, 1) = totalDLNewTxBytes * 8 ./ (obj.NumFrames * 0.01 * 1000 * 1000); % Mbps
            end
            % Uplink stats
            if ismember(obj.UplinkIdx, obj.PlotIds)
                ulStats(1, 1) = obj.PeakDataRateUL;
                ulStats(2, 1) = totalULTxBytes * 8 ./ (obj.NumFrames * 0.01 * 1000 * 1000); % Mbps
                ulStats(3, 1) = obj.PeakULSpectralEfficiency;
                ulStats(4, 1) = 1e6*ulStats(2, 1)/obj.Bandwidth(obj.UplinkIdx);
                ulStats(5, 1) = totalULNewTxBytes * 8 ./ (obj.NumFrames * 0.01 * 1000 * 1000); % Mbps
            end
        end
    end

    methods( Access = private)
        function itemList = constructRBItemList(obj, numRBs)
            %constructRBItemList Create the items for the drop-down component

            % Create the items for the drop-down component
            numItems = floor(numRBs / obj.RGMaxRBsToDisplay);
            itemList = cell(numItems, 1);
            for i = 1 : numItems
                k = (i - 1) * obj.RGMaxRBsToDisplay;
                itemList{i} = ['RB ', num2str(k) '-' num2str(k + obj.RGMaxRBsToDisplay - 1)];
            end
            if (mod(numRBs,obj.RGMaxRBsToDisplay) > 0)
                k = i * obj.RGMaxRBsToDisplay;
                itemList{i+1} = ['RB ', num2str(k) '-' num2str(numRBs - 1)];
            end
        end

        function rgSelectedRBRange(obj, dd, hAx)
            %rgSelectedRBRange Handle the event when user selects RB range in resource grid visualization

                obj.RGLowerRBIndex = obj.RGMaxRBsToDisplay * (dd.Value - 1);
                obj.RGUpperRBIndex = obj.RGLowerRBIndex + obj.RGMaxRBsToDisplay;
                if obj.RGUpperRBIndex > obj.NumRBs(obj.RVCurrView)
                    obj.RGUpperRBIndex =  obj.NumRBs(obj.RVCurrView);
                end
            % Update the Y-Axis of the resource grid visualization with
            % selected RB range
            replotResourceGrid(obj, hAx, 'YAxis');
        end

        function rgDropDownForSlotRange(obj, position, hAx)
            %rgDropDownForSlotRange Construct drop-down component for selecting slot range

            % Create the items for the drop-down component
            numItems = floor(obj.NumSlotsFrame / obj.RGMaxSlotsToDisplay);
            itemList = cell(numItems, 1);
            for i = 1 : numItems
                k = (i-1) * obj.RGMaxSlotsToDisplay ;
                itemList{i} = ['Slot ', num2str(k) '-' num2str(k + obj.RGMaxSlotsToDisplay - 1)];
            end
            if (mod(obj.NumSlotsFrame, obj.RGMaxSlotsToDisplay) > 0)
                k = i * obj.RGMaxSlotsToDisplay + 1;
                itemList{i+1} = ['Slot ', num2str(k - 1) '-' num2str(obj.NumSlotsFrame - 1)];
            end

            % Create drop-down component for slot display range selection
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', position, 'String', 'Select Slot', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [position(1)+0.090 position(2:4)], 'String', itemList, 'Callback', @(dd, event) rgSelectedSlotRange(obj, dd, hAx));
        end

        function rgSelectedSlotRange(obj, dd, hAx)
            %rgSelectedSlotRange Handle the event when user selects slot range in resource grid visualization

            obj.RGLowerSlotIndex = obj.RGMaxSlotsToDisplay * (dd.Value - 1);
            obj.RGUpperSlotIndex = obj.RGLowerSlotIndex + obj.RGMaxSlotsToDisplay;
            if obj.RGUpperSlotIndex > obj.NumSlotsFrame
                obj.RGUpperSlotIndex =  obj.NumSlotsFrame;
            end
            % Update the X-Axis of the resource grid visualization with
            % selected slot range
            replotResourceGrid(obj, hAx, 'XAxis');
        end

        function rbSelectedLinkType(obj, dd)
            %rbSelectedLinkType Handle the event when user selects link type in resource grid visualization

            % Update the resource grid visualization with selected link type
            if numel(obj.PlotIds) == 2
                obj.RVCurrView = dd.Value;
            end
            updateResourceGridVisualization(obj);
            drawnow;
        end

        function replotResourceGrid(obj, hAx, coordinate)
            %replotResourceGrid Update the resource grid along X-axis or Y-axis w.r.t to the given input parameters.

            cla(hAx);
            numRBsToDisplay = obj.RGUpperRBIndex - obj.RGLowerRBIndex;

            if obj.SchedulingType % For symbol based scheduling
                lowLogIdx = 0;
                numUnitsToDisplay = obj.NumSym; % Display information of 14 symbols in a slot
            else % For slot based scheduling
                lowLogIdx = obj.RGLowerSlotIndex;
                numUnitsToDisplay = obj.RGUpperSlotIndex - obj.RGLowerSlotIndex;
            end

            [X1, Y1] = meshgrid(0:numUnitsToDisplay, 0 : numRBsToDisplay);
            [X2, Y2] = meshgrid(0:numRBsToDisplay, 0 : numUnitsToDisplay);
            x = linspace(1, numUnitsToDisplay, numUnitsToDisplay);
            y = linspace(1, numRBsToDisplay, numRBsToDisplay);

            for n=1:numUnitsToDisplay
                i = lowLogIdx + n;
                for p = 1 : numRBsToDisplay
                    j = obj.RGLowerRBIndex + p;
                    obj.ResourceGridTextHandles(i, j) = text(hAx, x(n) - .5, y(p) - .5, ' ','FontName', 'FixedWidth', 'FontWeight', 'bold', 'FontUnits', 'normalized', 'HorizontalAlignment', 'center', 'Clipping', 'on');
                end
            end

            hold(hAx, 'on');
            plot(hAx, X1, Y1, 'k', 'Color', 'black', 'LineWidth', 0.1);
            plot(hAx, Y2, X2, 'k', 'Color', 'black', 'LineWidth', 0.1);

            if strcmpi('XAxis', coordinate) == 1
                % Updates X-Axis
                xticks(hAx, (1 : numUnitsToDisplay) - 0.5);
                if obj.SchedulingType % Symbol based scheduling
                    xticklabels(hAx, obj.SymbolInfo);
                else
                    xticklabels(hAx, obj.SlotInfo(obj.RGLowerSlotIndex+1 : obj.RGUpperSlotIndex));
                end
            else
                % Updates Y-Axis
                yticks(hAx, (1 : numRBsToDisplay) - 0.5);
                yTicksLabel = cell(1, 0);
                for i = 1 : numRBsToDisplay
                    yTicksLabel{i} = strcat('    RB- ', num2str(obj.RGLowerRBIndex + i - 1));
                end
                yticklabels(hAx, yTicksLabel);
            end

            % Update the resource grid visualization
            updateResourceGridVisualization(obj);
        end

        function updateResourceGridVisualization(obj)
            %updateResourceGridVisualization Update the resource grid visualization

            if obj.SchedulingType % For symbol based scheduling
                lowLogIdx =  0;
                uppLogIdx = obj.NumSym;
                % Update the axis
                obj.RGVisualizationFigHandle.CurrentAxes.XTickLabel = obj.SymbolInfo;
            else % For slot based scheduling
                lowLogIdx = obj.RGLowerSlotIndex;
                uppLogIdx = obj.RGUpperSlotIndex;
                % Update the axis
                obj.RGVisualizationFigHandle.CurrentAxes.XTickLabel = obj.SlotInfo(obj.RGLowerSlotIndex+1 : obj.RGUpperSlotIndex);
            end

            for n = lowLogIdx+1 : uppLogIdx
                for p = obj.RGLowerRBIndex + 1 : obj.RGUpperRBIndex
                    obj.ResourceGridTextHandles(n, p).String = obj.ResourceGridInfo{obj.RVCurrView}(n, p);
                    if(obj.ResourceGridReTxInfo{obj.RVCurrView}(n, p) == 2) % Re-Tx
                        obj.ResourceGridTextHandles(n, p).Color = 'blue';
                    else
                        obj.ResourceGridTextHandles(n, p).Color = 'black';
                    end
                end
            end

        end

        function cvDropDownForRBRange(obj, position, hAx)
            %cvDropDownForRBRange Construct drop-down component for selecting RB range

            numRBs = obj.NumRBs(obj.CVCurrView);
            % Create the items for the drop-down component
            numItems = floor(numRBs / obj.CVMaxRBsToDisplay);
            itemList = cell(numItems, 1);
            for i=1:numItems
                k = (i - 1) * obj.CVMaxRBsToDisplay;
                itemList{i} = ['RB ', num2str(k) '-' num2str(k + obj.CVMaxRBsToDisplay - 1)];
            end
            if (mod(obj.NumRBs(obj.CVCurrView),obj.CVMaxRBsToDisplay) > 0)
                k = i * obj.CVMaxRBsToDisplay;
                itemList{i+1} = ['RB ', num2str(k) '-' num2str(k+mod(numRBs,obj.CVMaxRBsToDisplay) - 1)];
            end

            obj.CQIVisualizationGridHandles = gobjects(obj.NumUEs, numRBs);

            % Create drop-down component for RB display range selection
            uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', position, 'String', 'Select RB range', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [position(1)+0.090 position(2:4)], 'String', itemList, 'Callback', @(dd, event) cbSelectedRBRange(obj, dd, hAx));
        end

        function cbSelectedRBRange(obj, dd, hAx)
            %cbSelectedRBRange Handle the event when user selects RB range in CQI grid visualization

                obj.CVLowerRBIndex = obj.CVMaxRBsToDisplay * (dd.Value - 1);
                obj.CVUpperRBIndex = obj.CVLowerRBIndex + obj.CVMaxRBsToDisplay;
                if obj.CVUpperRBIndex > obj.NumRBs(obj.CVCurrView)
                    obj.CVUpperRBIndex =  obj.NumRBs(obj.CVCurrView);
                end
            % Update the Y-Axis limits of the CQI grid visualization with
            % selected RB range
            replotCQIGrid(obj, hAx, 'YAxis');
        end

        function cvDropDownForUERange(obj, position, hAx)
            %cvDropDownForUERange Construct drop-down component for selecting UEs

            % Create the items for the drop-down component
            numItems = floor(obj.NumUEs / obj.CVMaxUEsToDisplay);
            itemList = cell(numItems, 1);
            for i = 1 : numItems
                k = (i - 1) * obj.CVMaxUEsToDisplay;
                itemList{i} = ['UE ', num2str(k + 1) '-' num2str(k + obj.CVMaxUEsToDisplay)];
            end
            if (mod(obj.NumUEs,obj.CVMaxUEsToDisplay) > 0)
                k = i * obj.CVMaxUEsToDisplay + 1;
                itemList{i+1} = ['UE ', num2str(k) '-' num2str(k + mod(obj.NumUEs, obj.CVMaxUEsToDisplay) - 1)];
            end

            obj.CQIVisualizationGridHandles = gobjects(obj.NumUEs, obj.NumRBs(obj.CVCurrView));
            % Create drop-down component for UE display range selection
            uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', position, 'String', 'Select UE', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [position(1)+0.090 position(2:4)], 'String', itemList, 'Callback', @(dd, event) cbSelectedUERange(obj, dd, hAx));
        end

        function cbSelectedUERange(obj, dd, hAx)
            %cbSelectedUERange Handle the event when user selects UE range in CQI grid visualization

                obj.CVLowerUEIndex = obj.CVMaxUEsToDisplay * (dd.Value - 1);
                obj.CVUpperUEIndex = obj.CVLowerUEIndex + obj.CVMaxUEsToDisplay;
                if obj.CVUpperUEIndex  > obj.NumUEs
                    obj.CVUpperUEIndex =  obj.NumUEs;
                end
            % Update the X-Axis limits of the CQI grid visualization with
            % selected UE range
            replotCQIGrid(obj, hAx, 'XAxis');
        end

        function cbSelectedLinkType(obj, dd)
            %cbSelectedLinkType Handle the event when user selects link type in CQI grid visualization

            % Update the CQI grid visualization with selected link type
            if numel(obj.PlotIds) == 2
                obj.CVCurrView = dd.Value;
            end
            updateCQIVisualization(obj);
            drawnow;
        end

        function replotCQIGrid(obj, hAx, coordinate)
            %replotCQIGrid Update the CQI grid along X-axis or Y-axis w.r.t to the given input parameters

            cla(hAx);
            numRBsToDisplay = obj.CVUpperRBIndex - obj.CVLowerRBIndex;
            numUEsToDisplay = obj.CVUpperUEIndex - obj.CVLowerUEIndex;
            for x = 1 : numUEsToDisplay
                for y = 1 : numRBsToDisplay
                    obj.CQIVisualizationGridHandles(x, y)  = rectangle(hAx,'Position',[x - 1 y - 1 1 1], 'FaceColor', 'white');
                end
            end

            if strcmpi('XAxis', coordinate)
                % Update X-Axis
                xticks(hAx,(1:numUEsToDisplay) - 0.5);
                xTicksLabel = cell(numUEsToDisplay, 0);
                for i = 1  : numUEsToDisplay
                    xTicksLabel{i} = strcat('UE- ', num2str(obj.CVLowerUEIndex + i ));
                end
                xticklabels(hAx, xTicksLabel);
            else
                % Update Y-Axis
                yticks(hAx, (1:numRBsToDisplay) - 0.5);
                yTicksLabel = cell(numRBsToDisplay, 0);
                for i = 1 : numRBsToDisplay
                    yTicksLabel{i} = strcat('    RB- ', num2str(obj.CVLowerRBIndex + i - 1));
                end
                yticklabels(hAx, yTicksLabel);
            end
            updateCQIVisualization(obj);
        end

        function updateCQIVisualization(obj)
            %updateCQIVisualization Update the CQI grid visualization

            numUEsToDisplay = obj.CVUpperUEIndex - obj.CVLowerUEIndex;
            numRBsToDisplay = obj.CVUpperRBIndex - obj.CVLowerRBIndex;
            for p = 1 : numRBsToDisplay
                % Update the color of the channel
                for q = 1 : numUEsToDisplay
                    rbIndex = obj.CVLowerRBIndex + p ;
                    obj.CQIVisualizationGridHandles(q, p).FaceColor = obj.ColorCoding{obj.CQIInfo{obj.CVCurrView}(q, rbIndex) + 1};
                end
            end
        end

        function logFormat = constructLogFormat(obj, linkIdx, simParam)
            %constructLogFormat Construct log format

            columnIndex = 1;
            logFormat{1, columnIndex} = 0; % Timestamp (in milliseconds)
            obj.ColumnIndexMap('Timestamp') = columnIndex;
            
            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = 0; % Frame number
            obj.ColumnIndexMap('Frame Number') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} =  0; % Slot number
            obj.ColumnIndexMap('Slot Number') = columnIndex;

            if(obj.SchedulingType == 1)
                % Symbol number column is only for symbol-based
                % scheduling
                columnIndex = columnIndex + 1;
                logFormat{1, columnIndex} =  0; % Symbol number
                obj.ColumnIndexMap('Symbol Number') = columnIndex;
            end
            if(obj.DuplexMode == obj.TDDDuplexMode)
                % Slot/symbol type as DL/UL/guard is only for TDD mode
                columnIndex = columnIndex + 1;
                logFormat{1, columnIndex} = 'Guard'; % Symbol type
                obj.ColumnIndexMap('Type') = columnIndex;
            end
            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = zeros(obj.NumUEs, ceil(obj.NumRBs(linkIdx) / obj.RBGSize(linkIdx))); % RBG allocation for UEs
            obj.ColumnIndexMap('RBG Allocation Map') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1*ones(obj.NumUEs, 1); % MCS for assignments
            obj.ColumnIndexMap('UEs MCS') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1*ones(obj.NumUEs, 1); % HARQ IDs for assignments
            obj.ColumnIndexMap('HARQ Process ID') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1*ones(obj.NumUEs, 1); % NDI flag for assignments
            obj.ColumnIndexMap('Grant NDI Flag') = columnIndex;

            % Tx type of the assignments ('newTx' or 'reTx'), 'noTx' if there is no assignment
            txTypeUEs =  cell(obj.NumUEs, 1);
            txTypeUEs(:) = {'noTx'};
            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = txTypeUEs;
            obj.ColumnIndexMap('Transmission') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = obj.CQIInfo{linkIdx}; % Uplink channel quality
            obj.ColumnIndexMap('Channel Quality') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = zeros(obj.NumUEs, obj.NumHARQ); % HARQ process status
            obj.ColumnIndexMap('HARQ process NDI status (at symbol start)') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = zeros(obj.NumUEs, 1); % MAC bytes transmitted
            obj.ColumnIndexMap('Throughput Bytes') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = zeros(obj.NumUEs, 1); % MAC bytes corresponding to a new transmission
            obj.ColumnIndexMap('Goodput Bytes') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = zeros(obj.NumUEs, 1); % UEs' buffer status
            obj.ColumnIndexMap('Buffer status of UEs (In bytes)') = columnIndex;

            % To store reception result from UEs for this slot ('rxSuccess'
            % or 'rxFailure'). Result is 'noRx', if UE was not scheduled to transmit
            rxResultUEs =  cell(obj.NumUEs, 1);
            rxResultUEs(:) = {'noRx'};

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = rxResultUEs;
            obj.ColumnIndexMap('Rx Success/Failure') = columnIndex;

            % Initialize scheduling log for all the symbols in the
            % simulation time. The last time scheduler runs in the
            % simulation, it might assign resources for future slots which
            % are outside of simulation time. Storing those decisions too
            numSlotsSim = simParam.NumFramesSim * obj.NumSlotsFrame; % Simulation time in units of slot duration
            logFormat = repmat(logFormat(1,:), (numSlotsSim + obj.NumSlotsFrame)*obj.NumSym , 1);
        end

        function logFormat = constructGrantLogFormat(obj, simParam)
            %constructGrantLogFormat Construct grant log format

            columnIndex = 1;
            logFormat{1, columnIndex} = -1; % UE's RNTI
            obj.GrantLogsColumnIndexMap('RNTI') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % Frame number
            obj.GrantLogsColumnIndexMap('Frame') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % Slot number
            obj.GrantLogsColumnIndexMap('Slot') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = {''}; % Type: UL or DL
            obj.GrantLogsColumnIndexMap('Grant type') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = {''}; % RBG allocation for UEs
            obj.GrantLogsColumnIndexMap('RBG Allocation Map') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % Start sym
            obj.GrantLogsColumnIndexMap('Start Sym') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % Num sym
            obj.GrantLogsColumnIndexMap('Num Sym') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % MCS Value
            obj.GrantLogsColumnIndexMap('MCS') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % Number of layers
            obj.GrantLogsColumnIndexMap('NumLayers') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % HARQ IDs for assignments
            obj.GrantLogsColumnIndexMap('HARQ ID') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % NDI flag for assignments
            obj.GrantLogsColumnIndexMap('NDI Flag') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = -1; % RV for assignments
            obj.GrantLogsColumnIndexMap('RV') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = {''}; % Tx type: new-Tx or re-Tx
            obj.GrantLogsColumnIndexMap('Tx Type') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = {'NA'}; % PDSCH feedback slot offset (Only applicable for DL grants)
            obj.GrantLogsColumnIndexMap('Feedback Slot Offset (DL grants only)') = columnIndex;

            columnIndex = columnIndex + 1;
            logFormat{1, columnIndex} = {''}; % CQI values
            obj.GrantLogsColumnIndexMap('CQI on RBs') = columnIndex;

            % Initialize scheduling log for all the symbols in the
            % simulation time. The last time scheduler runs in the
            % simulation, it might assign resources for future slots which
            % are outside of simulation time. Storing those decisions too
            if(obj.SchedulingType == 1 && isfield(simParam, 'TTIGranularity'))
                maxRows = obj.NumFrames*obj.NumSlotsFrame*obj.NumUEs*(ceil(obj.NumSym/simParam.TTIGranularity));
            else
                maxRows = obj.NumFrames*obj.NumSlotsFrame*obj.NumUEs;
            end
            logFormat = repmat(logFormat(1,:), maxRows , 1);
        end

        function constructCQIGridVisualization(obj, screenWidth, screenHeight, figureWidth, figureHeight)
            %constructCQIGridVisualization Construct CQI grid visualization

            % Number of RBs to be displayed in the default view of CQI visualization
            maxRBs = max(obj.NumRBs(obj.PlotIds));
            if obj.CVMaxRBsToDisplay <= maxRBs
                obj.CVUpperRBIndex = obj.CVMaxRBsToDisplay;
            else
                obj.CVUpperRBIndex = maxRBs;
            end

            % Number of UEs to be displayed in the default view of CQI visualization
            if obj.NumUEs >= obj.CVMaxUEsToDisplay
                obj.CVUpperUEIndex = obj.CVMaxUEsToDisplay;
            else
                obj.CVUpperUEIndex = obj.NumUEs;
            end

            % Create the figure for CQI Visualization
            obj.CQIVisualizationFigHandle = figure('Name', 'Channel Quality Visualization', 'Position', [screenWidth * 0.05 screenHeight * 0.05 figureWidth figureHeight]);
            % Coordinates of UI control elements displayed in the figure
            dropDownXCoordinate = 0.015;
            dropDownYCoordinate = 0.75;

            % Create label for frame number
            obj.CVTxtHandle = uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [dropDownXCoordinate dropDownYCoordinate 0.13 0.03], 'String', 'Frame Number: ', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');

            % Set axis properties
            hAx = gca;
            set(hAx, 'Units', 'Pixels', 'Position', [figureWidth * 0.28 figureHeight * 0.1 figureWidth * 0.67 figureHeight * 0.8], 'Units', 'normalized');
            obj.CQIVisualizationGridHandles = gobjects(obj.NumSlotsFrame, max(obj.NumRBs(obj.PlotIds)));
            hold on

            % Create title and legend
            title(strcat("Channel Quality Visualization for Cell ID - ", num2str(obj.NCellID)), 'FontSize', 15,'FontUnits', 'normalized', 'Position', [0.5 1.08 0], 'Units', 'normalized');
            legendStr = {'CQI 0-1', 'CQI 2-3', 'CQI 4-5', ...
                'CQI 6-7', 'CQI 8-9', 'CQI 10-11', ...
                'CQI 12-13', 'CQI 14-15'};
            t = gobjects(1, 8);
            t(1) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{1});
            t(2) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{3});
            t(3) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{5});
            t(4) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{7});
            t(5) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{9});
            t(6) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{11});
            t(7) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{13});
            t(8) = plot(NaN, 'square', 'Color', [0 0 0], 'MarkerSize', 10, 'MarkerFaceColor', obj.ColorCoding{15});
            legend(t, legendStr, 'NumColumns', 8, 'Box', 'off', 'Location', 'northoutside', 'FontSize', 10, 'Orientation', 'horizontal', 'AutoUpdate', 'off','Units', 'normalized');
            drawnow;
            % Start coordinates for UI control elements displayed in the figure
            dropDownXCoordinate = 0.015;
            dropDownYCoordinate = dropDownYCoordinate - 2 * 0.03;

            % Create drop down for link type
            if numel(obj.PlotIds) == 2
                uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [dropDownXCoordinate dropDownYCoordinate 0.06 0.025], 'String', 'Select Link', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                uicontrol(obj.CQIVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [dropDownXCoordinate+0.060 dropDownYCoordinate 0.06 0.025], 'String', {'Downlink', 'Uplink'}, 'Callback', @(dd, event) cbSelectedLinkType(obj, dd));
                dropDownYCoordinate = dropDownYCoordinate - 2 * 0.03;
            end

            % Create drop-down for RB range
            vector = [dropDownXCoordinate dropDownYCoordinate 0.09 0.025];
            if obj.CVMaxRBsToDisplay < maxRBs
                cvDropDownForRBRange(obj, vector, hAx);
                dropDownYCoordinate = dropDownYCoordinate - 2 * 0.03;
            end

            % Create drop-down for UE range
            vector = [dropDownXCoordinate dropDownYCoordinate 0.09 0.025];
            if obj.NumUEs > obj.CVMaxUEsToDisplay
                cvDropDownForUERange(obj, vector, hAx);
            end

            % Set CQI-visualization axis label
            replotCQIGrid(obj, hAx, 'XAxis');
            xlabel(hAx, 'UEs', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'FontUnits', 'normalized');
            replotCQIGrid(obj, hAx, 'YAxis');
            ylabel(hAx, 'Resource Blocks', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'FontUnits', 'normalized');
            hAx.TickDir = 'out';
        end

        function constructResourceGridVisualization(obj, screenWidth, screenHeight, figureWidth, figureHeight)
            %constructResourceGridVisualization Construct resource grid
            %visualization

            % Number of RBs to be displayed in the default view of resource grid visualization
            maxRBs = max(obj.NumRBs(obj.PlotIds));
            if obj.RGMaxRBsToDisplay <= maxRBs
                obj.RGUpperRBIndex = obj.RGMaxRBsToDisplay;
            else
                obj.RGUpperRBIndex = maxRBs;
            end

            % Number of slots to be displayed in the default view of resource grid visualization
            if obj.NumSlotsFrame >= obj.RGMaxSlotsToDisplay
                obj.RGUpperSlotIndex = obj.RGMaxSlotsToDisplay;
            else
                obj.RGUpperSlotIndex = obj.NumSlotsFrame;
            end

            % Construct the drop-down item list
            for idx = 1:min(obj.NumLogs, numel(obj.PlotIds))
                if obj.DuplexMode == obj.FDDDuplexMode
                    plotId = obj.PlotIds(idx);
                else
                    plotId = idx;
                end
                % Construct the drop down based on number of RBs
                if obj.RGMaxRBsToDisplay < obj.NumRBs(plotId)
                    obj.RBItemsList{plotId} = constructRBItemList(obj, obj.NumRBs(plotId));
                end
            end

            % Create the figure for resource grid visualization
            obj.RGVisualizationFigHandle = figure('Name', 'Resource Grid Allocation', 'Position', [screenWidth * 0.05 screenHeight * 0.05 figureWidth figureHeight]);

            % Starting coordinates for UI control elements displayed in the figure
            xCoordinate = 0.015;
            yCoordinate = 0.85;
            % Gap between each UI element vertically
            gap = 0.03;

            % Create legend
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'UE-x(n) : Transmission', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            yCoordinate = yCoordinate - gap;
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'UE-x(n) : Retransmission', 'FontSize', 10, 'ForegroundColor', 'blue', 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            yCoordinate = yCoordinate - gap;
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'x : UE RNTI', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            yCoordinate = yCoordinate - gap;
            uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'n : HARQ Process ID', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
            yCoordinate = yCoordinate - gap;

            % Set axis properties
            hAx = gca;
            set(hAx, 'Units', 'Pixels', 'Position', [figureWidth * 0.28 figureHeight * 0.1 figureWidth * 0.67 figureHeight * 0.8], 'Units', 'normalized');
            obj.ResourceGridTextHandles  = gobjects(obj.NumSlotsFrame, maxRBs);

            % Create title
            title(strcat("Resource Grid Allocation for Cell ID - ", num2str(obj.NCellID)), 'FontSize', 15,'FontUnits', 'normalized', 'Position', [0.5 1.02 0], 'Units', 'normalized');

            % If post simulation log analysis enabled
            if isempty(obj.IsLogReplay) || obj.IsLogReplay
                % Create label for frame number
                obj.RGTxtHandle = uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'Frame Number: ', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                if obj.SchedulingType % Symbol based scheduling
                    % Create label for slot number
                    yCoordinate = yCoordinate - gap;
                    obj.RGSlotTxtHandle = uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', 'Slot Number: ', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                end
                yCoordinate = yCoordinate - gap;
            else
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.13 0.03], 'String', ['Total Frames: ' num2str(obj.NumFrames)], 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                yCoordinate = yCoordinate - gap;
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.09 0.02], 'String', 'Enter Frame Number', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'edit', 'String', ' ', 'Units', 'normalized', 'Position', [xCoordinate+0.09 yCoordinate 0.03 0.02], 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'Callback',@(H, E) showFrame(obj, H));
                yCoordinate = yCoordinate - gap;
                if obj.SchedulingType % Symbol based scheduling
                    uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.09 0.02], 'String', 'Enter Slot Number', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                    uicontrol(obj.RGVisualizationFigHandle, 'Style', 'edit', 'String', ' ', 'Units', 'normalized', 'Position', [xCoordinate+0.09 yCoordinate 0.03 0.02], 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'Callback',@(H, E) showSlot(obj, H));
                    yCoordinate = yCoordinate - gap;
                    obj.CurrFrame  = 0;
                    obj.CurrSlot = 0;
                end
                drawnow;
            end

            if obj.SchedulingType
                % Initialize the symbol pattern in a slot
                for sidx =1:obj.NumSym
                    obj.SymbolInfo{sidx} = strcat("Symbol-", num2str(sidx-1));
                end
            else
                % Initialize the slot pattern in a frame
                for sidx =1:obj.NumSlotsFrame
                    obj.SlotInfo{sidx} = strcat("Slot-", num2str(sidx-1));
                end
            end

            % Set resource-grid visualization axis label
            replotResourceGrid(obj, hAx, 'XAxis');
            if obj.SchedulingType
                xlabel(hAx, 'Symbols in Slot', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'FontUnits', 'normalized');
            else
                xlabel(hAx, 'Slots in 10 ms Frame', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'FontUnits', 'normalized');
            end
            replotResourceGrid(obj, hAx, 'YAxis');
            ylabel(hAx, 'Resource Blocks', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'FontUnits', 'normalized');
            hAx.TickDir = 'out';

            % Create drop-down for link type
            if min(obj.NumLogs, numel(obj.PlotIds))== 2
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.09 0.025], 'String', 'Select Link', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [xCoordinate+0.090 yCoordinate 0.06 0.025], 'String', {'Downlink', 'Uplink'}, 'Callback', @(dd, event) rbSelectedLinkType(obj, dd));
                yCoordinate = yCoordinate - gap;
            end

            if obj.RGMaxRBsToDisplay < maxRBs
                % Create drop-down component for RB display range selection
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'text', 'Units', 'normalized', 'Position', [xCoordinate yCoordinate 0.09 0.025], 'String', 'Select RB range', 'FontSize', 10, 'HorizontalAlignment', 'left', 'FontUnits', 'normalized');
                uicontrol(obj.RGVisualizationFigHandle, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [xCoordinate+0.090 yCoordinate 0.09 0.025], 'String', obj.RBItemsList{obj.RVCurrView}, 'Callback', @(dd, event) rgSelectedRBRange(obj, dd, hAx));
                yCoordinate = yCoordinate - gap;
            end

            if obj.SchedulingType == obj.SlotBased
                % Create drop-down for Slot range
                vector = [xCoordinate yCoordinate  0.09 0.025];
                if obj.NumSlotsFrame > obj.RGMaxSlotsToDisplay
                    rgDropDownForSlotRange(obj, vector, hAx);
                end
            end
        end

        function showFrame(obj, h)
            %showFrame Handle the event when user enters a
            % number to visualize a particular frame number in the
            % simulation

            frameNumber = str2double(get(h,'string'));
            if isnan(frameNumber)
                warndlg('Frame number must be an integer','Warning');
                return;
            end
            if frameNumber >= obj.NumFrames || frameNumber < 0
                msg = strcat("Frame number must be in between 0 and ", num2str(obj.NumFrames - 1));
                warndlg(msg,'Warning');
                return;
            end

            if obj.SchedulingType == obj.SlotBased % Slot based scheduling
                % Set the current slot as last slot of the frame
                obj.CurrSlot = obj.NumSlotsFrame - 1;
            end
            obj.CurrFrame = frameNumber;

            % Update the resource grid and CQI grid visualization
            plotRBGrids(obj);
            plotCQIRBGrids(obj);
        end

        function showSlot(obj, h)
            %showSlot Handle the event when user enters a
            % number to visualize a particular slot number in the
            % simulation

            slotNumber = str2double(get(h,'string'));
            if isnan(slotNumber)
                warndlg('Slot number must be an integer','Warning');
                return;
            end
            if slotNumber >= obj.NumSlotsFrame || slotNumber < 0
                msg = strcat("Slot number must be in between 0 and ", num2str(obj.NumSlotsFrame - 1));
                warndlg(msg,'Warning');
                return;
            end

            obj.CurrSlot = slotNumber;
            % Update the resource grid and CQI grid visualization
            plotRBGrids(obj);
            plotCQIRBGrids(obj);
        end
    end
end