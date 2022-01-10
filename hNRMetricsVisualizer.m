classdef hNRMetricsVisualizer < handle
%hNRMetricsVisualizer Creates metrics visualization object
%   The class implements visualization of the metrics. The following three types of
%   visualizations are shown:
%       (i) Display of RLC metrics
%      (ii) Display of MAC Scheduler performance metrics
%     (iii) Display of Phy metrics
        
%   Copyright 2020-2021 The MathWorks, Inc.

    properties
        % SimTime Simulation time (in seconds)
        SimTime
        
        %CellOfInterest Cell id to which the visualization belongs
        CellOfInterest

        % UELegend Legend for the UE
        UELegend

        % NumMetricsSteps Number of times metrics plots are updated
        NumMetricsSteps

        % RLCVisualization Timescope to display UE RLC layer's logical channel throughput
        RLCVisualization
        
        % MACVisualization Timescope to display the downlink and uplink scheduler performance metrics
        MACVisualization = cell(2, 1);

        % PhyVisualization Timescope to display the downlink and uplink block error rates
        PhyVisualization

        % LCHInfo Logical channel information. It is structure and contains
        % following fields.
        %   RNTI - Radio network temporary identifier of a UE
        %   LCID - Specifies the logical channel id of a UE
        %   EntityDir - Specifies the logical channel type
        %   corresponding to the logical channel specified in LCID.
        %      - 0 represents the logical channel is in downlink direction
        %      - 1 represents the logical channel is in uplink direction
        %      - 2 represents the logical channel is in both downlink & uplink direction
        LCHInfo

        % MetricsStepSize Number of slots in one metrics step
        MetricsStepSize

        % MetricsStepDuration Duration of 1 metrics step
        MetricsStepDuration

        %RLCLogger RLC logger handle object
        RLCLogger
        
        %MACLogger MAC logger handle object
        MACLogger

        %PhyLogger Phy logger handle object
        PhyLogger

        % RLCMetricName Name of the RLC metric to plot
        % It can take one of the values 'TxDataPDU', 'TxDataBytes',
        % 'ReTxDataPDU', 'ReTxDataBytes', 'TxControlPDU', 'TxControlBytes',
        % 'TxPacketsDropped', 'TimerPollRetransmitTimedOut', 'RxDataPDU',
        % 'TxBytesDropped', 'RxDataPDUBytes', 'RxDataPDUDropped',
        % 'RxDataBytesDropped', 'RxDataPDUDuplicate', 'RxDataBytesDuplicate',
        % 'RxControlPDU', 'RxControlBytes', 'TimerReassemblyTimedOut',
        % 'TimerStatusProhibitTimedOut'. The default value is 'TxDataBytes'
        RLCMetricName = 'TxDataBytes'
        
        % RLCMetricCustomLabel Custom label name to display in y-axis for the RLC metrics
        RLCMetricCustomLabel = 'Transmitted Bytes';
        
        % PeakDataRateDL Theoretical peak data rate 
        % A vector of two elements. First element represents downlink and
        % second element represents uplink theoretical peak data rate respectively
        PeakDataRate = zeros(2, 1);
        
        %YXC begin
        siteIdx
        
        YUO
        %YXC end
    end

    properties (Access = private)
        % VisualizationFlag  Indicates the plots to visualize
        % It takes the values 0, 1, 2 and represent downlink, uplink, and both
        % respectively. Default value is 2.
        VisualizationFlag = 2

        %PlotIds Represent the IDs of the plots
        PlotIDs = [1 2]

        % UEOfInterestListInfo Information about the list of UEs of interest
        % It is a M-by-3 matrix, where M represents the number of
        % UEs and the column1, column2, and column3 represents the UE id,
        % number of logical channels in downlink, uplink respectively
        UEOfInterestListInfo
    end

    properties (Access = private, Constant, Hidden)
        % NumLogicalChannels Maximum number of logical channels in each UE
        NumLogicalChannels = 32;

        % Constants related to downlink and uplink information. These
        % constants are used for indexing logs and identifying plots
        % DownlinkIdx Index for all downlink information
        DownlinkIdx = 1;
        % UplinkIdx Index for all uplink information
        UplinkIdx = 2;
    end

    methods (Access = public)
        function obj = hNRMetricsVisualizer(param, varargin)
            %hNRMetricsVisualizer Constructs metrics visualization object
            %
            % OBJ = hNRMetricsVisualizer(PARAM) Create metrics visualization
            % object for downlink and uplink plots.
            %
            % OBJ = hNRMetricsVisualizer(PARAM, Name, Value) creates a metrics visualization
            % object, OBJ, with properties specified by one or more name-value
            % pairs. You can specify additional name-value pair arguments in any
            % order as (Name1,Value1,...,NameN,ValueN).
            %
            % PARAM - It is a structure and contain simulation
            % configuration information.
            %
            %    NumFramesSim      - Number of frames in simulation
            %    CellOfInterest    - Cell of interest
            %    SCS               - Subcarrier spacing
            %    UEOfInterest      - List of UEs of interest
            %    NumUEs            - Number of UEs
            %    NumMetricsSteps   - Number of times metrics plots to be
            %                        updated
            %    MetricsStepSize   - Interval at which metrics visualization
            %                        updates in terms of number of slots
            %
            % VARARGIN - Optional arguments as name-value pairs

            % Initialize the properties
            for idx = 1:2:numel(varargin)
                obj.(varargin{idx}) = varargin{idx+1};
            end
            obj.SimTime = (10 * param.NumFramesSim) / 1000; % Simulation time (in seconds)
            if isempty(obj.NumMetricsSteps)
                obj.NumMetricsSteps = param.NumMetricsSteps;
            end
            % Interval at which metrics visualization updates in terms of number of
            % slots. Make sure that MetricsStepSize is an integer
            if isempty(obj.MetricsStepSize)
                obj.MetricsStepSize = param.MetricsStepSize;
            end
            obj.MetricsStepDuration = obj.MetricsStepSize * (15 / param.SCS);
            
            if isempty(obj.CellOfInterest)
                if isfield(param, 'CellOfInterest')
                    obj.CellOfInterest = param.CellOfInterest;
                else
                    obj.CellOfInterest = 1;
                end
            end
            if isfield(param, 'UEOfInterest')
                ueOfInterestList = param.UEOfInterest;
            elseif ~isempty(obj.LCHInfo)
                ueOfInterestList = [obj.LCHInfo.RNTI];
            else
                ueOfInterestList = 1:param.NumUEs;
            end
            % Create legend information for the plots
            numUEs = numel(ueOfInterestList);
            obj.UELegend = cell(1, numUEs);
            obj.UEOfInterestListInfo = zeros(numUEs, 3);
            for idx = 1:numUEs
                obj.UEOfInterestListInfo(idx, 1) = ueOfInterestList(idx); % Update the UE id
                obj.UELegend{idx} = ['UE-' num2str(ueOfInterestList(idx)) ' '];
            end

            % Determine the plots
            if obj.VisualizationFlag == 0
                obj.PlotIDs = obj.DownlinkIdx; % Downlink plot id
            elseif obj.VisualizationFlag == 1
                obj.PlotIDs = obj.UplinkIdx; % Uplink plot id
            end

            % Create RLC visualization
            if ~isempty(obj.RLCLogger) && ~isempty(obj.LCHInfo)
                addRLCVisualization(obj);
            end

            % Create Phy visualization
            if ~isempty(obj.PhyLogger)
                addPhyVisualization(obj);
            end

            % Create MAC visualization
            if ~isempty(obj.MACLogger)
                addMACVisualization(obj);
            end
        end

        function addRLCVisualization(obj, varargin)
            %addRLCVisualization Create RLC visualization
            %
            % addRLCVisualization(OBJ) Create and configure RLC
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink.
            %
            % addRLCVisualization(OBJ, RLCLOGGER) Create and configure RLC
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink and also set the RLCLogger value
            %
            % RLCLOGGER - RLC logger. It is an object of type hNRRLCLogger

            % Create the timescope
            if isempty(obj.RLCVisualization)
                obj.RLCVisualization = timescope('Name', 'RLC Metrics Visualization');
            end
            
            % Set RLCLogger
            if nargin == 2
                obj.RLCLogger = varargin{1};
            end

            lchNames = cell(1,1);
            count = 0;
            % Create the logical channel names for legend
            for ueIdx=1:size(obj.UEOfInterestListInfo, 1)
                ueId = obj.UEOfInterestListInfo(ueIdx, 1);
                idx = find(ueId == [obj.LCHInfo.RNTI], 1);

                % Logical channels in downlink
                dlIdx = sort([find(obj.LCHInfo(idx).EntityDir == 0); find(obj.LCHInfo(idx).EntityDir == 2)]);
                dlLogicalChannels = obj.LCHInfo(idx).LCID(dlIdx);
                obj.UEOfInterestListInfo(ueIdx, obj.DownlinkIdx+1) = numel(dlLogicalChannels); % Number of channels
                for lcIdx=1:numel(dlLogicalChannels)
                    count = count + 1;
                    lchNames{count, 1} = [obj.UELegend{ueIdx} 'LCH-' num2str(dlLogicalChannels(lcIdx))];
                end
            end
            for ueIdx=1:size(obj.UEOfInterestListInfo, 1)
                ueId = obj.UEOfInterestListInfo(ueIdx, 1);
                idx = find(ueId == [obj.LCHInfo.RNTI], 1);
                % Logical channels in uplink
                ulIdx = sort([find(obj.LCHInfo(idx).EntityDir == 1); find(obj.LCHInfo(idx).EntityDir == 2)]);
                ulLogicalChannels = obj.LCHInfo(idx).LCID(ulIdx);
                obj.UEOfInterestListInfo(ueIdx, obj.UplinkIdx+1) = numel(ulLogicalChannels); % Number of channels
                for lcIdx=1:numel(ulLogicalChannels)
                    count = count + 1;
                    lchNames{count, 1} = [obj.UELegend{ueIdx} 'LCH-' num2str(ulLogicalChannels(lcIdx))];
                end
            end

            % Update the timescope properties
            release(obj.RLCVisualization);
            set(obj.RLCVisualization, 'LayoutDimensions', [numel(obj.PlotIDs) 1], 'ShowLegend', true, ...
                'SampleRate', obj.NumMetricsSteps/obj.SimTime, 'TimeSpanSource', 'property','ChannelNames', lchNames, 'TimeSpan', obj.SimTime);

            % Initialize the plots
            if numel(obj.PlotIDs) == 1
                obj.RLCVisualization(zeros(1, sum(obj.UEOfInterestListInfo(:, obj.DownlinkIdx+1) + obj.UEOfInterestListInfo(:, obj.UplinkIdx+1))));
            else
                obj.RLCVisualization(zeros(1, sum(obj.UEOfInterestListInfo(:, obj.DownlinkIdx+1))), zeros(1, sum(obj.UEOfInterestListInfo(:, obj.UplinkIdx+1))));
            end

            % Add the titles and legends
            for idx=1:numel(obj.PlotIDs)
                obj.RLCVisualization.ActiveDisplay = idx;
                obj.RLCVisualization.YLabel = ['Cell-' num2str(obj.CellOfInterest) ' ' obj.RLCMetricCustomLabel];
                obj.RLCVisualization.AxesScaling = 'Updates';
                obj.RLCVisualization.AxesScalingNumUpdates = 1;

                if obj.PlotIDs(idx) == obj.DownlinkIdx
                    obj.RLCVisualization.Title = 'Downlink Logical Channels (LCH)';
                else
                    obj.RLCVisualization.Title = 'Uplink Logical Channels (LCH)';
                end
            end

        end
        
        function addMACVisualization(obj, varargin)
            %addMACVisualization Create MAC visualization
            %
            % addMACVisualization(OBJ) Create and configure MAC
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink.
            %
            % addMACVisualization(OBJ, MACLOGGER) Create and configure MAC
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink and also set the MACLogger value
            %
            % MACLOGGER - MAC logger. It is an object of type hNRSchedulingLogger
            
            % Set MACLogger
            if nargin == 2
                obj.MACLogger = varargin{1};
            end

            numUEs = size(obj.UEOfInterestListInfo, 1);
            nodeMetrics = zeros(1, numUEs);
            % Plot titles and Y-axis label prefix
            title = {'Downlink Scheduler Performance Metrics', ...
                'Uplink Scheduler Performance Metrics'};
            tag = {['Cell-' num2str(obj.CellOfInterest) ' DL '], ...
                ['Cell-' num2str(obj.CellOfInterest) ' UL ']};
            channelNames = [obj.UELegend 'Cell' 'Peak Data Rate' obj.UELegend obj.UELegend 'Cell' 'Peak Data Rate' obj.UELegend];
            
            % Set peak data rate
            if isempty(obj.MACLogger.PeakDataRateDL)
                obj.PeakDataRate(obj.DownlinkIdx) = 0;
            else
                obj.PeakDataRate(obj.DownlinkIdx) = obj.MACLogger.PeakDataRateDL;
            end
            if isempty(obj.MACLogger.PeakDataRateUL)
                obj.PeakDataRate(obj.UplinkIdx) = 0;
            else
                obj.PeakDataRate(obj.UplinkIdx) = obj.MACLogger.PeakDataRateUL;
            end
            % Create time scope and add labels
            for idx=1:numel(obj.PlotIDs)
                windowId = obj.PlotIDs(idx);

                if isempty(obj.MACVisualization{windowId})
                    obj.MACVisualization{windowId} = timescope('Name', title{windowId});
                end

                release(obj.MACVisualization{windowId});
                set(obj.MACVisualization{windowId}, 'LayoutDimensions',[2 2], 'ChannelNames', channelNames,...
                    'ActiveDisplay',1, 'YLabel',[tag{windowId} 'Throughput (Mbps)'], 'ShowLegend',true,'AxesScaling', 'Updates', ...
                     'AxesScalingNumUpdates', 1, 'TimeSpanSource', 'property', 'TimeSpan', obj.SimTime, ...
                    'ActiveDisplay',2, 'YLabel',[tag{windowId} 'Resource Share (%)'], ...
                    'ShowLegend',true, 'YLimits',[1 100],'AxesScaling', 'Updates','AxesScalingNumUpdates', 1, ...
                    'SampleRate', obj.NumMetricsSteps/obj.SimTime, 'TimeSpanSource', 'property', 'TimeSpan', obj.SimTime, ...
                    'ActiveDisplay',3, 'YLabel',[tag{windowId} 'Goodput (Mbps)'], 'ShowLegend',true,'AxesScaling', 'Updates', 'AxesScalingNumUpdates', 1, ...
                    'SampleRate', obj.NumMetricsSteps/obj.SimTime, 'TimeSpanSource', 'property', 'TimeSpan', obj.SimTime, ...
                    'ActiveDisplay',4, 'YLabel',[tag{windowId} 'Buffer Status (KB)'], 'ShowLegend',true,'AxesScaling', 'Updates', 'AxesScalingNumUpdates', 1, ...
                    'SampleRate', obj.NumMetricsSteps/obj.SimTime, 'TimeSpanSource', 'property', 'TimeSpan', obj.SimTime);
                obj.MACVisualization{windowId}([nodeMetrics 0 obj.PeakDataRate(windowId)], nodeMetrics, [nodeMetrics 0 obj.PeakDataRate(windowId)], nodeMetrics);
            end
        end

        function addPhyVisualization(obj, varargin)
            %addPhyVisualization Create Phy visualization
            %
            % addPhyVisualization(OBJ) Create and configure Phy
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink.
            %
            % addPhyVisualization(OBJ, PHYLOGGER) Create and configure Phy
            % visualization. It creates figures for visualizing metrics
            % in both downlink and uplink and also set the PhyLogger value
            %
            % PHYLOGGER - Phy logger. It is an object of type hNRPhyLogger

            % Set PhyLogger
            if nargin == 2
                obj.PhyLogger = varargin{1};
            end

            % Create and configure the timescope
            if isempty(obj.PhyVisualization)
                obj.PhyVisualization = timescope('Name', 'Block Error Rate (BLER) Visualization');
            end
            release(obj.PhyVisualization);
            set(obj.PhyVisualization, 'LayoutDimensions', [numel(obj.PlotIDs) 1], 'ShowLegend', true, ...
                'SampleRate', obj.NumMetricsSteps/obj.SimTime,'TimeSpanSource', 'property','ChannelNames', repmat(obj.UELegend, [1 numel(obj.PlotIDs)]), 'TimeSpan', obj.SimTime);

            numUEs = size(obj.UEOfInterestListInfo, 1);
            blerData = zeros(1, numUEs);
            titles = {'Downlink BLER', 'Uplink BLER'};
            % Initialize the plots
            if numel(obj.PlotIDs) == 1
                obj.PhyVisualization(blerData);
            else
                obj.PhyVisualization(blerData, blerData);
            end

            % Add the titles and legends
            for idx=1:numel(obj.PlotIDs)
                obj.PhyVisualization.ActiveDisplay = idx;
                obj.PhyVisualization.YLimits = [0 1];
                obj.PhyVisualization.YLabel = ['Cell-' num2str(obj.CellOfInterest) ' BLER'];
                obj.PhyVisualization.Title = titles{obj.PlotIDs(idx)};
            end

        end

        function plotMetrics(obj, slotNum)
            %plotMetrics Updates the metric plots
            %
            % plotMetrics(OBJ, SLOTNUM) Updates the metrics plots
            %
            % SLOTNUM - Slot number in simulation. It is used to calculate
            % the time interval corresponding to which metrics plots have
            % to be updated

            % RLC metrics visualization
            if ~isempty(obj.RLCVisualization)
                plotRLCMetrics(obj, slotNum);
            end

            % MAC metrics visualization
            if ~isempty(obj.MACVisualization{1}) || ~isempty(obj.MACVisualization{2})
                plotMACMetrics(obj, slotNum);
            end

            % PHY metrics visualization
            if ~isempty(obj.PhyVisualization)
                plotPhyMetrics(obj, slotNum);
            end
        end
    end

    methods(Access = private)
        function plotRLCMetrics(obj, slotNum)
            %plotRLCMetrics Plots the RLC metrics
            %
            % plotRLCMetrics(OBJ, SLOTNUM) Plots the metrics of each logical
            % channel of each UE
            
            dlIdx = obj.DownlinkIdx+1; % DL index
            ulIdx = obj.UplinkIdx+1; % UL index
            numDLChannels = sum(obj.UEOfInterestListInfo(:, dlIdx));
            numULChannels = sum(obj.UEOfInterestListInfo(:, ulIdx));
            metricInfo = zeros(2, max(numULChannels, numDLChannels));
            dlCount = 1;
            ulCount = 1;
            [dlMetrics, ulMetrics] = getRLCMetrics(obj.RLCLogger, ...
                slotNum-obj.MetricsStepSize+1, slotNum, obj.UEOfInterestListInfo(:, 1), obj.RLCMetricName);

            for ueIdx=1:size(obj.UEOfInterestListInfo, 1)
                if ~isempty(dlMetrics)
                    metricInfo(obj.DownlinkIdx, dlCount:dlCount+obj.UEOfInterestListInfo(ueIdx, dlIdx)-1) = dlMetrics(ueIdx).MetricValue;
                    dlCount = dlCount+obj.UEOfInterestListInfo(ueIdx, dlIdx);
                end

                if ~isempty(ulMetrics)
                    metricInfo(obj.UplinkIdx, ulCount:ulCount+obj.UEOfInterestListInfo(ueIdx, ulIdx)-1) = ulMetrics(ueIdx).MetricValue;
                    ulCount = ulCount+obj.UEOfInterestListInfo(ueIdx, ulIdx);
                end
            end

            % Update the plots
            if numel(obj.PlotIDs) == 1
                obj.RLCVisualization(metricInfo(obj.PlotIDs(1), :));
            else
                obj.RLCVisualization(metricInfo(obj.DownlinkIdx, 1:numDLChannels), metricInfo(obj.UplinkIdx, 1:numULChannels));
            end
        end

        function plotMACMetrics(obj, slotNum)
            %plotMACMetrics Plot the MAC metrics
            %
            % plotMACMetrics(OBJ, SLOTNUM) Plots the metrics of each UE

            numUEs = numel(obj.UELegend);
            throughputServed = zeros(2, numUEs+2);
            goodput = zeros(2, numUEs+2);
            bufferstatus = zeros(2, numUEs);
            resourceshare = zeros(2, numUEs);

            [dlMetrics, ulMetrics, cellMetrics] = getMACMetrics(obj.MACLogger, slotNum-obj.MetricsStepSize+1, slotNum, obj.UEOfInterestListInfo(:, 1));
            if ~isempty(dlMetrics)
                throughputServed(obj.DownlinkIdx, 1:numUEs) = [dlMetrics.TxBytes] .* 8 ./ (obj.MetricsStepDuration * 1000);
                throughputServed(obj.DownlinkIdx, numUEs+1) = cellMetrics.DLTxBytes .* 8 ./ (obj.MetricsStepDuration * 1000); % Cell throughput
                throughputServed(obj.DownlinkIdx, numUEs+2) = obj.PeakDataRate(obj.DownlinkIdx); % Peak datarate
                goodput(obj.DownlinkIdx, 1:numUEs) = [dlMetrics.NewTxBytes] .* 8 ./ (obj.MetricsStepDuration * 1000);
                goodput(obj.DownlinkIdx, numUEs+1) = cellMetrics.DLNewTxBytes .* 8 ./ (obj.MetricsStepDuration * 1000); % Cell goodput
                goodput(obj.DownlinkIdx, numUEs+2) = obj.PeakDataRate(obj.DownlinkIdx); % Peak datarate
                bufferstatus(obj.DownlinkIdx, 1:numUEs) = [dlMetrics.BufferStatus] ./ 1000; % In KB
                resourceshare(obj.DownlinkIdx, 1:numUEs) = ([dlMetrics.AssignedRBCount] * 100) ./ cellMetrics.DLRBsScheduled; % Percent share
            end

            if ~isempty(ulMetrics)
                throughputServed(obj.UplinkIdx, 1:numUEs) = [ulMetrics.TxBytes] .* 8 ./ (obj.MetricsStepDuration * 1000);
                throughputServed(obj.UplinkIdx, numUEs+1) = cellMetrics.ULTxBytes .* 8 ./ (obj.MetricsStepDuration * 1000); % Cell throughput
                throughputServed(obj.UplinkIdx, numUEs+2) = obj.PeakDataRate(obj.UplinkIdx); % Peak datarate
                goodput(obj.UplinkIdx, 1:numUEs) = [ulMetrics.NewTxBytes] .* 8 ./ (obj.MetricsStepDuration * 1000);
                goodput(obj.UplinkIdx, numUEs+1) = cellMetrics.ULNewTxBytes .* 8 ./ (obj.MetricsStepDuration * 1000); % Cell goodput
                goodput(obj.UplinkIdx, numUEs+2) = obj.PeakDataRate(obj.UplinkIdx); % Peak datarate
                bufferstatus(obj.UplinkIdx, 1:numUEs) = [ulMetrics.BufferStatus] ./ 1000; % In KB
                resourceshare(obj.UplinkIdx, 1:numUEs) = ([ulMetrics.AssignedRBCount] * 100) ./ cellMetrics.ULRBsScheduled; % Percent share
            end

            % Update the plots
            for plotIdx = 1:numel(obj.PlotIDs)
                plotId = obj.PlotIDs(plotIdx);
                obj.MACVisualization{plotId}(throughputServed(plotId, :), resourceshare(plotId, :), goodput(plotId, :), bufferstatus(plotId, :));
            end
            
            %YXC begin
            % Store throughput in YUO
            storeThroughput(obj.YUO, throughputServed, obj.siteIdx, slotNum);
            %YXC end
        end

        function plotPhyMetrics(obj, slotNum)
            %plotPhyMetrics Plot the Phy metrics
            %
            % plotPhyMetrics(OBJ, SLOTNUM) Plots the metrics of each UE

            numUEs = size(obj.UEOfInterestListInfo, 1);
            blerData = zeros(2, numUEs);

            % Get Phy metrics
            [dlMetrics, ulMetrics] = getPhyMetrics(obj.PhyLogger, slotNum-obj.MetricsStepSize+1, slotNum, obj.UEOfInterestListInfo(:, 1));

            if ~isempty(dlMetrics)
                blerData(obj.DownlinkIdx, :) = [dlMetrics.ErroneousPackets] ./ [dlMetrics.TotalPackets];
            end
            if ~isempty(ulMetrics)
                blerData(obj.UplinkIdx, :) = [ulMetrics.ErroneousPackets] ./ [ulMetrics.TotalPackets];
            end

            blerData(isnan(blerData)) = 0; % To handle NaN
            % Update the plots
            if numel(obj.PlotIDs) == 1
                obj.PhyVisualization(blerData(obj.PlotIDs(1), :));
            else
                obj.PhyVisualization(blerData(obj.DownlinkIdx, :), blerData(obj.UplinkIdx, :));
            end
        end
    end
end