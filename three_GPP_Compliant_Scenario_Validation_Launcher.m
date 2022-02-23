% this launcher file is based on the modified version of NRInterferenceModelingWithToroidalWrapAroundExample
% Changes have been made to make the simulation compliant with 3GPP specifications
% This launcher is based on 3GPP TR37.910 V16.1.0, Annex 
% detailed calibration parameters and assumptions are found in section 4 of RP-180524
%


%% clear workspace, clear command window, close all figures
% clear
% clc
% close all

%start timer
tic

%MXC_2
disp('Simulation starting...');
disp('Loading simulation parameter configurations...');

%% Simulation parameter configuration
rng('default'); % Reset the random number generator
simParameters = []; % Clear the simParameters variable

% simulation configuration
simParameters.NumFramesSim = 2; % Simulation time, in number of 10 ms frames
simParameters.EnableWrapAround = true; % Enable wrap-around modeling
simParameters.Scenario='RMa'; %UMi, UMa or RMa
simParameters.ScenarioConfiguration = 'A'; % A, B
simParameters.ChannelModelType = 'CDL';
simParameters.SchedulingType = 0; % Slot-based scheduling
simParameters.EnableAllVisualization = false;%true;

%MXC_2
switch simParameters.Scenario
    case 'UMa'
        switch simParameters.ScenarioConfiguration
            case 'A'
                % Carrier frequency, channel bandwidth, and subcarrier spacing (SCS)
                % See 3GPP TS 38.104 section 5.3.2
                simParameters.SCS = 15; % kHz
                simParameters.NumRBs = 52; 
                % 15KHz SCS: 25->5MHz, 52->10MHz, 106->20MHz, 216->40MHz
                % 30Khz SCS: 11->5MHz, 24->10Mhz, 51->20MHz, 106->40MHz
                % see TR37.910 Table 8.1.1-2 for detail
                % Assume that the UL and DL carriers have the same channel bandwidth
                simParameters.DLBandwidth = 10e6; % Hz
                simParameters.ULBandwidth = 10e6; % Hz

                simParameters.DLCarrierFreq = 4e9; % Hz
                simParameters.ULCarrierFreq = 4e9; % Hz
                %see 3GPP TR37.910 Table 8.2.1-1 for detailed NR bands
                %also found in Table 8.2.1-1 and 8.2.1-2 in TS38.104
                simParameters.DuplexMode = 1; %0 - FDD; 1 - TDD
                
                simParameters.InterSiteDistance = 200; % Distance between adjacent gNBs in meters
                simParameters.minUEgNBDistance = 10;
                
                % gNB configuration
                simParameters.gNBHeight = 25; % meters
                simParameters.AntennaDowntilt = 6; % degrees
                simParameters.AntennaSlant = 0; % degrees
                simParameters.GNBTxPower = 41; % Tx power for gNBs in dBm
                simParameters.GNBRxGain = 8; % Receiver antenna gain at gNB in dBi
                simParameters.GNBTxAnts = 128;
                simParameters.GNBRxAnts = 1;
                simParameters.GNBTxAntPanelSize = [8 8 2 1 1]; %[M N P Mg Ng]
                simParameters.GNBRxAntPanelSize = [1 1 1 1 1]; %[M N P Mg Ng]
                simParameters.GNBTxAntElementSpacing = [0.5 0.8 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.GNBRxAntElementSpacing = [0.5 0.8 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.GNBTxAntPolarizationAngles = [45 -45];
                simParameters.GNBRxAntPolarizationAngles = 45;
                simParameters.GNBAntElement = '38.901';
                simParameters.GNBAntPolarizationModel = 'Model-2';
                
                % UE configuration
                simParameters.NumUEsCell = 2; % Number of UEs in each cell
                simParameters.UEHeight = 1.5; % meters
                simParameters.UETxPower = 23; % Tx power for all the UEs in dBm
                simParameters.UETxAnts = 1;
                simParameters.UERxAnts = 4;
                simParameters.UETxAntPanelSize = [1 1 1 1 1];
                simParameters.UERxAntPanelSize = [1 2 2 1 1]; %[M N P Mg Ng]
                simParameters.UETxAntElementSpacing = [0.5 0.5 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.UERxAntElementSpacing = [0.5 0.5 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.UETxAntPolarizationAngles = 0;
                simParameters.UERxAntPolarizationAngles = [0 90];
                simParameters.UEAntElement = 'isotropic';
                simParameters.UEAntPolarizationModel = 'Model-2';
                
                
                
                
                
            case 'B'
                disp('UMa_B');
                error('NOT YET!');
            otherwise
            error('NO!');
        end
    case 'RMa'
        switch simParameters.ScenarioConfiguration
            case 'A'
                % Carrier frequency, channel bandwidth, and subcarrier spacing (SCS)
                % See 3GPP TS 38.104 section 5.3.2
                simParameters.SCS = 15; % kHz
                simParameters.NumRBs = 52; 
                % 15KHz SCS: 25->5MHz, 52->10MHz, 106->20MHz, 216->40MHz
                % 30Khz SCS: 11->5MHz, 24->10Mhz, 51->20MHz, 106->40MHz
                % see TR37.910 Table 8.1.1-2 for detail
                % Assume that the UL and DL carriers have the same channel bandwidth
                simParameters.DLBandwidth = 10e6; % Hz
                simParameters.ULBandwidth = 10e6; % Hz

                simParameters.DLCarrierFreq = 710e6; % Hz
                simParameters.ULCarrierFreq = 690e6; % Hz
                %see 3GPP TR37.910 Table 8.2.1-1 for detailed NR bands
                %also found in Table 8.2.1-1 and 8.2.1-2 in TS38.104
                simParameters.DuplexMode = 0; %0 - FDD; 1 - TDD
                %should be FDD, but TDD for now
                
                simParameters.InterSiteDistance = 1732; % Distance between adjacent gNBs in meters
                simParameters.minUEgNBDistance = 10;
                
                % gNB configuration
                simParameters.gNBHeight = 35; % meters
                simParameters.AntennaDowntilt = 10; % degrees
                simParameters.AntennaSlant = 0; % degrees
                simParameters.GNBTxPower = 46; % Tx power for gNBs in dBm
                simParameters.GNBRxGain = 8; % Receiver antenna gain at gNB in dBi
                simParameters.GNBTxAnts = 64;
                simParameters.GNBRxAnts = 1;
                %YXC begin
                simParameters.GNBTxAntPanelSize = [8 4 2 1 1]; %[M N P Mg Ng]
                %YXC end
                simParameters.GNBRxAntPanelSize = [1 1 1 1 1]; %[M N P Mg Ng]
                simParameters.GNBTxAntElementSpacing = [0.5 0.8 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.GNBRxAntElementSpacing = [0.5 0.8 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.GNBTxAntPolarizationAngles = [45 -45];
                simParameters.GNBRxAntPolarizationAngles = 45;
                simParameters.GNBAntElement = '38.901';
                simParameters.GNBAntPolarizationModel = 'Model-2';
                
                
                % UE configuration
                simParameters.NumUEsCell = 10; % Number of UEs in each cell
                simParameters.UEHeight = 1.5; % meters
                simParameters.UETxPower = 23; % Tx power for all the UEs in dBm
                simParameters.UETxAnts = 1;
                simParameters.UERxAnts = 2;
                simParameters.UETxAntPanelSize = [1 1 1 1 1];
                simParameters.UERxAntPanelSize = [1 1 2 1 1]; %[M N P Mg Ng]
                simParameters.UETxAntElementSpacing = [0.5 0.5 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.UERxAntElementSpacing = [0.5 0.5 1 1]; % [dH dV dgv dgh] vertical and horzontal element spacing and panel spacing
                simParameters.UETxAntPolarizationAngles = 0;
                simParameters.UERxAntPolarizationAngles = [0 90];
                simParameters.UEAntElement = 'isotropic';
                simParameters.UEAntPolarizationModel = 'Model-2';
            case 'B'
                disp('RMa_B');
                error('NOT YET!');
            otherwise
            error('NO!');
        end
    otherwise
    error('NO!');
end

% Fast fading tables
load FastFadingTabs.mat; % Load 3GPP TR 38.901 Table 7.5-6, 7.5-7, 7.5-8, 7.5-9, 7.5-10, 7.5-11 to workspace
simParameters.FastFadingTabs = FastFadingTabs; % Functions cannot use the variables from workspace directly




% Specify the signal-to-interference-plus-noise ratio (SINR) to CQI mapping table for a block error rate (BLER) of 0.1.
simParameters.SINR90pc = [-5.46 -0.46 4.54 9.05 11.54 14.04 15.54 18.04 ...
    20.04 22.43 24.93 25.43 27.43 30.43 33.43];

% Specify the scheduling strategy and the maximum limit on the RBs allotted for PUSCH and PDSCH.
% The transmission limit applies only to new PUSCH and PDSCH assignments, and not to the retransmissions.
simParameters.SchedulerStrategy = 'RR'; % Supported scheduling strategies: 'PF', 'RR', and 'BestCQI'
simParameters.RBAllocationLimitUL = 20; % For PUSCH
simParameters.RBAllocationLimitDL = 20; % For PDSCH 

%% Logging and Visualization Configuration
% Specify the ID of cell of interest as an integer between 1 and 56, inclusive. 
% The example shows visualizations and metrics for this cell and the central cell, cell-0.
simParameters.CellOfInterest = 2; %Set a value from 1 to 56. Set the value to 0 to visualize cell-0 only.
% column vector to be consistant  
% i dont understand why this needs to be a scalar, comment validateattributes out for now
validateattributes(simParameters.CellOfInterest, {'numeric'}, {'integer', 'scalar', '>=', 0, '<=', 56}, 'simParameters.CellOfInterest', 'CellOfInterest')
% The CQIVisualization and RBVisualization parameters control the display of the CQI visualization of RBs and the RB assignment visualization. 
% To enable the RB visualization plot, set the RBVisualization field to true.
%MXC_3
if simParameters.EnableAllVisualization
    simParameters.CQIVisualization = true;
    simParameters.RBVisualization = true;
else
    simParameters.CQIVisualization = false;
    simParameters.RBVisualization = false;
end
%MXC_3
%Update the output metrics plots periodically, specifying NumMetricsSteps updates within the simulation.
simParameters.NumMetricsSteps = 10;
%Write the logs to MAT-files. The example uses these logs for post-simulation analysis and visualization.
simParameters.ParametersLogFile = 'simParameters'; % For logging the simulation parameters
simParameters.SimulationLogFile = 'simulationLogs'; % For logging the simulation logs

%% Application traffic configuration
% Automatically generates the parameters in accordance to the number of UEs
ulPacketPeriodicity = 10 * ones(simParameters.NumUEsCell,1); % Periodicity, in ms, at which the UEs generate UL packets
ulPacketSize = 2e4 * ones(simParameters.NumUEsCell,1); % Size of generated UL packets, in bytes
dlPacketPeriodicity = 10 * ones(simParameters.NumUEsCell,1); % Periodicity, in ms, at which gNB generates DL packets
dlPacketSize = 2e4 * ones(simParameters.NumUEsCell,1); % Size of generated DL packets, in bytes

%% Derived Parameters
% Compute the derived parameters based on the primary configuration parameters specified in the previous section, and set some example-specific constants.


simParameters.NumClusters = 1;
simParameters.NumSitesPerCluster = 57; % Number of gNBs per cluster

% Set the UE and gNB positions.
[simParameters.gNBBearing, simParameters.GNBPositions, simParameters.UEPositions, simParameters.CellPositions, simParameters.UEStates] = hMacrocellTopology(simParameters);
simParameters.NCellIDList = 0:simParameters.NumSitesPerCluster-1; % List of physical cell IDs

%MXC_2 begin
%{
% Specify the CSI-RS resource configuration, assuming that all UEs measure channel quality on the same CSI-RS resource.
%simParameters.CSIRSRowNumber = 2; % Possible row numbers for single transmit antenna case are 1 and 2
simParameters.CSIRSRowNumber = 12; 
simParameters.SubbandSize = 8; % Size of sub-band for CQI reporting in terms of number of RBs
%}

%copied from mimo example
%YXC begin
% Temporarily make the number of antenna ports and the number of antenna
% elements equal, as 5G Toolbox does not support flexible layer mapping
simParameters.CSIRSRowNumber = 11; % CSI-RS row number as per 3GPP TS 38.211 Table 7.4.1.5.3-1
simParameters.CSIRSSubcarrierLocation = [1 3 5 7];
simParameters.CSIRSSymbolLocation = 0;
simParameters.CSIRSPeriod = [4 0];  % Slot periodicity and offset
simParameters.PMIMode = 'Subband';  % 'Wideband' or 'Subband'
simParameters.CQIMode = 'Subband';  % 'Wideband' or 'Subband'
simParameters.PanelDimensions = [4 2]; % [N1 N2] as per 3GPP TS 38.214 Table 5.2.2.2.12
simParameters.SubbandSize = 4; % Refer TS 38.214 Table 5.2.1.4-2 for valid subband sizes
simParameters.CodebookMode = 2; % 1 or 2
simParameters.RankIndicator = 2; 
% Copied from https://www.sharetechnote.com/html/lte_toolbox/Matlab_LteToolbox_5G_CSI_RS.html#Example_p8_Ex01
% simParameters.CSIRSRowNumber = 6; 
% simParameters.CSIRSSubcarrierLocation = [2 4 6 8];
% simParameters.CSIRSSymbolLocation = 3;
% simParameters.CSIRSPeriod = [40 1];
% simParameters.PMIMode = 'Subband';
% simParameters.CQIMode = 'Subband';
% simParameters.PanelDimensions = [2 2]; % This is different from the one shown in the URL
% simParameters.SubbandSize = 4;
% simParameters.CodebookMode = 1;
% simParameters.RankIndicator = 2; 
%{
% Change DL to SISO for testing
% Configuration copied from the original example
% simParameters.CSIRSRowNumber = 2; % Possible row numbers for single transmit antenna case are 1 and 2
% simParameters.SubbandSize = 8; % Size of sub-band for CQI reporting in terms of number of RBs
% simParameters.PanelDimensions = [1,1];
%}
%YXC end
%MXC_2 end

% Calculate the slot duration for the selected SCS and the number of slots in a 10 ms frame.
slotDuration = 1/(simParameters.SCS/15); % In ms
numSlotsFrame = 10/slotDuration; % Number of slots in a 10 ms frame
numSlotsSim = simParameters.NumFramesSim*numSlotsFrame; % Number of slots in the simulation

% Set the interval at which the example updates metrics visualization, number of slots. 
% Because this example uses a time granularity of one slot, the MetricsStepSize field must be an integer.
simParameters.MetricsStepSize = ceil(numSlotsSim / simParameters.NumMetricsSteps);
if mod(numSlotsSim, simParameters.NumMetricsSteps) ~= 0
    % Update the NumMetricsSteps parameter if it does not exactly divide NumSlotsSim
    simParameters.NumMetricsSteps = floor(numSlotsSim / simParameters.MetricsStepSize);
end

% Specify one logical channel in each UE, and set the logical channel configuration for all nodes (UEs and gNBs) in the example.
numLogicalChannels = 1;
simParameters.LCHConfig.LCID = 4; % Logical channel ID (logical channel ID of data radio bearers starts from 4)

% Specify the RLC entity direction as 0 to specify DL only, 1 to specify UL only, or 2 to specify both UL and DL.
simParameters.RLCConfig.EntityDir = 2;

% Create RLC channel configuration structure.
rlcChannelConfigStruct.LCGID = 1; % Mapping between logical channel and logical channel group ID
rlcChannelConfigStruct.Priority = 1; % Priority of each logical channel
rlcChannelConfigStruct.PBR = 8; % Prioritized bitrate (PBR), in kilobytes per second, of each logical channel
rlcChannelConfigStruct.BSD = 10; % Bucket size duration (BSD), in ms, of each logical channel
rlcChannelConfigStruct.EntityType = simParameters.RLCConfig.EntityDir;
rlcChannelConfigStruct.LogicalChannelID = simParameters.LCHConfig.LCID;

% Set the maximum RLC service data unit (SDU) length, in bytes, as specified in 3GPP TS 38.323, and the total number of nodes in the simulation.
simParameters.maxRLCSDULength = 9000;
simParameters.MaxReceivers = simParameters.NumSitesPerCluster*(simParameters.NumUEsCell + 1);

if ~isfield(simParameters, 'SchedulingType') || simParameters.SchedulingType == 0 % No scheduling type or slot-based scheduling
    rbAssignmentPlotPeriodicity = numSlotsFrame; % Update RB assignment visualization every frame (10 ms)
    tickGranularity = 14;
    simParameters.PUSCHMappingType = 'A';
    simParameters.PDSCHMappingType = 'A';
else % Symbol-based scheduling
    rbAssignmentPlotPeriodicity = 1; % Update RB assignment visualization every slot
    tickGranularity = 1;
    simParameters.PUSCHMappingType = 'B';
    simParameters.PDSCHMappingType = 'B';
end

%% Multicell Setup
% Set up the 57-cell topology.
% For each cell, create the gNB and UE objects, initialize the channel quality information for UEs, and set up the logical channel at gNB and UEs. 
% The hNRGNB and hNRUE helper classes create gNB and UE nodes respectively, containing the RLC, MAC, and PHY layers.

%MXC_2
disp('Generating topology...');
gNB = cell(simParameters.NumSitesPerCluster, 1);
UEs = cell(simParameters.NumSitesPerCluster, simParameters.NumUEsCell);

%Create DL and UL packet distribution objects, and initialize a wrap-around distance calculator callback.
dlPacketDistributionObj = hNRPacketDistribution(simParameters, 0); % 0 for DL
ulPacketDistributionObj = hNRPacketDistribution(simParameters, 1); % 1 for UL
if simParameters.EnableWrapAround
    distCalc = @(TxPos,RxPos) wDistCalc(simParameters.InterSiteDistance/3,TxPos,RxPos);
else
    distCalc = @(TxPos,RxPos) wDistCalc(simParameters.InterSiteDistance/3,TxPos,RxPos,'NoWrapAround');
end

%Initialize YusUtilityObj
YUO = YusUtilityObj(simParameters, numSlotsSim);

for siteIdx = 1:simParameters.NumSitesPerCluster
    simParameters.NCellID = simParameters.NCellIDList(siteIdx); % Cell ID
    simParameters.Position = simParameters.GNBPositions(siteIdx, :);
    simParameters.NumUEs = simParameters.NumUEsCell; % Number of UEs in a cell
    % Create scheduler
    switch(simParameters.SchedulerStrategy)
        case 'RR' % Round robin scheduler
            scheduler = hNRSchedulerRoundRobin(simParameters);
        case 'PF' % Proportional fair scheduler
            scheduler = hNRSchedulerProportionalFair(simParameters);
        case 'BestCQI' % Best CQI scheduler
            scheduler = hNRSchedulerBestCQI(simParameters);
    end
    
    % Create gNB
    gNB{siteIdx} = hNRGNB(simParameters);
    addScheduler(gNB{siteIdx}, scheduler); % Add scheduler to gNB
    gNB{siteIdx}.PhyEntity = hNRGNBPhy(simParameters, siteIdx); % Create PHY layer instance
    configurePhy(gNB{siteIdx}, simParameters); % Configure PHY layer
    % Register distance calculator for wrap-around distance computations
    gNB{siteIdx}.DistanceCalculatorFcn = distCalc;
    setPhyInterface(gNB{siteIdx}); % Set up the interface to PHY layer

    for ueIdx = 1:simParameters.NumUEsCell
        simParameters.Position = simParameters.UEPositions{siteIdx}(ueIdx, :); % Position of UE in (x,y,z) coordinates
        UEs{siteIdx, ueIdx} = hNRUE(simParameters, ueIdx);
        UEs{siteIdx, ueIdx}.PhyEntity = hNRUEPhy(simParameters, siteIdx, ueIdx); % Create PHY layer instance
        UEs{siteIdx, ueIdx}.PhyEntity.StoreCQIInfo = @YUO.storeCQIInfo; % Register the callback to YusUtilityObj
        UEs{siteIdx, ueIdx}.PhyEntity.YusUtilityParameter.YUO = YUO;    % Store handle of YUO in UE PHY
        configurePhy(UEs{siteIdx, ueIdx}, simParameters); % Configure PHY layer
        % Register distance calculator for wrap-around distance computations
        UEs{siteIdx, ueIdx}.DistanceCalculatorFcn =  distCalc;
        setPhyInterface(UEs{siteIdx, ueIdx}); % Set up the interface to PHY

        % Set up logical channel at gNB for the UE
        configureLogicalChannel(gNB{siteIdx}, ueIdx, rlcChannelConfigStruct);
        % Set up logical channel at UE
        configureLogicalChannel(UEs{siteIdx, ueIdx}, ueIdx, rlcChannelConfigStruct);

        % Add data traffic pattern generators to gNB and UE nodes
        % Calculate the data rate (in kbps) of On-Off traffic pattern using
        % packet size (in bytes) and packet interval (in ms)
        ulDataRate = ceil(1000/ulPacketPeriodicity(ueIdx)) * ulPacketSize(ueIdx) * 8e-3;
        % Limit the size of the generated application packet to the maximum RLC
        % SDU size. The maximum supported RLC SDU size is 9000 bytes
        if ulPacketSize(ueIdx) > simParameters.maxRLCSDULength
            ulPacketSize(ueIdx) = simParameters.maxRLCSDULength;
        end
        % Create an object for On-Off network traffic pattern and add it to the
        % specified UE. This object generates the uplink data traffic on the UE
        ulApp = networkTrafficOnOff('PacketSize', ulPacketSize(ueIdx), 'GeneratePacket', true, ...
            'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, 'DataRate', ulDataRate);
        UEs{siteIdx, ueIdx}.addApplication(ueIdx, simParameters.LCHConfig.LCID, ulApp);
        
        dlDataRate = ceil(1000/dlPacketPeriodicity(ueIdx)) * dlPacketSize(ueIdx) * 8e-3;
        if dlPacketSize(ueIdx) > simParameters.maxRLCSDULength
            dlPacketSize(ueIdx) = simParameters.maxRLCSDULength;
        end
        % Create an object for On-Off network traffic pattern for the specified
        % UE and add it to the gNB. This object generates the downlink data
        % traffic on the gNB for the UE
        dlApp = networkTrafficOnOff('PacketSize', dlPacketSize(ueIdx), 'GeneratePacket', true, ...
            'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, 'DataRate', dlDataRate);
        gNB{siteIdx}.addApplication(ueIdx, simParameters.LCHConfig.LCID, dlApp);
    end
    
    % Setup the UL and DL packet distribution mechanism
    hNRSetUpPacketDistribution(simParameters, gNB{siteIdx}, UEs(siteIdx, :), dlPacketDistributionObj, ulPacketDistributionObj);
end

% Set up logging and visualization, specifying the central cell (cell 0) and the cell of interest.
%MXC_2
%cellsOfInterest = unique([0; simParameters.CellOfInterest]);
cellsOfInterest = (0:56)';%[0; 1; 2];
numCellsOfInterest = length(cellsOfInterest); % Number of cells that the example logs and visualizes
%YXC begin
YUO.cellOfInterestIdx = numCellsOfInterest;
%YXC end

% Visualize the network topology
hTopologyVisualizer(simParameters);   

% Log and visualize PHY and MAC metrics in cell arrays.
simSchedulingLogger = cell(numCellsOfInterest, 1);
% Cell array for PHY metrics logging and visualization
simPhyLogger = cell(numCellsOfInterest, 1);
visualizer = cell(numCellsOfInterest, 1);

for siteIdx = 1:numCellsOfInterest
    simParameters.NCellID = cellsOfInterest(siteIdx);
    simParameters.CellOfInterest = simParameters.NCellID;
    % Create an object for MAC scheduling information visualization and logging
    simSchedulingLogger{siteIdx} = hNRSchedulingLogger(simParameters); 
    % Create an object for PHY layer metrics logging
    simPhyLogger{siteIdx} = hNRPhyLogger(simParameters);
    % Create visualization object for MAC and PHY metrics 
    %comment out visualization code to speed up simulation
    
    if simParameters.EnableAllVisualization
        %YXC begin
        visualizer{siteIdx} = hNRMetricsVisualizer(simParameters, 'MACLogger', simSchedulingLogger{siteIdx}, 'PhyLogger', simPhyLogger{siteIdx},'siteIdx',siteIdx,'YUO',YUO);
        %YXC end
    end
end
%YXC begin
ConfTime = toc;
ConfTime = seconds(ConfTime);
disp(['Time elapsed: ',datestr(ConfTime,'MM:SS')])
disp('Configuration finished. Start processing loop.')
%YXC end
%% Processing Loop
% Simulation is run slot by slot. For each cell, in each slot, these operations are executed:
% *Run the MAC and PHY layers of gNB
% *Run the MAC and PHY layers of UEs
% *Log and visualize metrics for each layer
% *Advance the timer for the nodes and send a trigger to application and RLC layers every millisecond. 
% *The application layer and RLC layer execute their scheduled operations based on a 1 ms timer trigger.


slotNum = 0;
numSymbolsSim = numSlotsSim * 14; % Simulation time in units of symbol duration
% Run processing loop
for symbolNum = 1 : tickGranularity : numSymbolsSim

    %MXC_2
    %show which symbol is currently being run
    %fprintf('running slot # %d out of %d \n',(((symbolNum-1)/14)+1), numSlotsSim);
    %MXC_2

    if mod(symbolNum - 1, 14) == 0
        slotNum = slotNum + 1;
    end
    
    % Because all the cells operate on the same SCS, slot durations do not vary
    for siteIdx = 1:simParameters.NumSitesPerCluster % For each site
        run(gNB{siteIdx});
        
        % Run MAC and PHY layers of UEs
        for ueIdx = 1:simParameters.NumUEsCell

            %Push indices into YusUtilityObj
            pushTriple(YUO,slotNum,siteIdx,ueIdx);
            
            run(UEs{siteIdx, ueIdx});
            %YXC begin
            % Log simulation progress
            Time = toc;
            Time = seconds(Time);
            LoopTime = Time - ConfTime;
            Progress = ((symbolNum-1)/tickGranularity*simParameters.NumSitesPerCluster*...
                simParameters.NumUEsCell+(siteIdx-1)*simParameters.NumUEsCell+ueIdx)/...
                (numSlotsSim*simParameters.NumSitesPerCluster*simParameters.NumUEsCell);
            TimeLeft = LoopTime/Progress + ConfTime - LoopTime;
            disp(['Time elapsed: ',datestr(Time,'dd:HH:MM:SS'), '. Remaining time: ',datestr(TimeLeft,'dd:HH:MM:SS')])
            disp(['Progress: ',num2str(floor(100*Progress)),'%. ',...
                'Processed: UE(',num2str(ueIdx),'/',num2str(simParameters.NumUEsCell),') in ',...
                'cell(',num2str(siteIdx),'/',num2str(simParameters.NumSitesPerCluster),') of ',...
                'slot(',num2str(floor((symbolNum-1)/tickGranularity)+1),'/',num2str(numSlotsSim),').'])
            %YXC end
        end

        cellIdx = find((siteIdx-1) == cellsOfInterest, 1);
        if ~isempty(cellIdx)
            % MAC logging
            logCellSchedulingStats(simSchedulingLogger{cellIdx}, symbolNum, gNB{siteIdx}, UEs(siteIdx, :));

            % PHY logging
            logCellPhyStats(simPhyLogger{cellIdx}, symbolNum, gNB{siteIdx}, UEs(siteIdx, :));
        end
    end
    
    for siteIdx = 1:simParameters.NumSitesPerCluster
        % Advance timer ticks for gNB and UEs by the number of symbols per slot
        advanceTimer(gNB{siteIdx}, tickGranularity);
        for ueIdx = 1:simParameters.NumUEsCell

            %Push indices into YusUtilityObj
            pushTriple(YUO,slotNum,siteIdx,ueIdx);

            advanceTimer(UEs{siteIdx, ueIdx}, tickGranularity);
        end
    end
    %MXC_3
    if simParameters.EnableAllVisualization
        for idx = 1:numCellsOfInterest
            % Visualization
            % Check slot boundary
            if symbolNum > 1 && ((simParameters.SchedulingType == 1 && mod(symbolNum, 14) == 0) || (simParameters.SchedulingType == 0 && mod(symbolNum-1, 14) == 0))
                % RB assignment visualization (if enabled)
                if simParameters.RBVisualization
                    if mod(slotNum, rbAssignmentPlotPeriodicity) == 0
                        % Plot at slot boundary, if the update periodicity is reached
                        plotRBGrids(simSchedulingLogger{idx});
                    end
                end
                % CQI grid visualization (if enabled)
                if simParameters.CQIVisualization
                    if mod(slotNum, numSlotsFrame) == 0 % Plot at frame boundary
                        plotCQIRBGrids(simSchedulingLogger{idx});
                    end
                end
                % If the update periodicity is reached, plot scheduler metrics and PHY metrics visualization
                % at slot boundary
                if mod(slotNum, simParameters.MetricsStepSize) == 0
                    plotMetrics(visualizer{idx}, slotNum);
                end
            end
        end
    end
    
    %MXC_3
    
    
    %yuxtime = toc;
    %lineToPrint = ['Elapsed time: ', datestr(seconds(yuxtime),'dd:HH:MM:SS')];
    %disp(lineToPrint);

end

for idx = 1:numCellsOfInterest
    [dlStats, ulStats] = getPerformanceIndicators(simSchedulingLogger{idx});
    [logInfo.BLERLogs, logInfo.AvgBLERLogs] = getBLERLogs(simPhyLogger{idx}); % Block Error rate logs
    fprintf('\n\nMetrics for cell %d :\n\n', cellsOfInterest(idx));
    fprintf('Peak UL throughput: %0.2f Mbps. Achieved average UL Throughput: %0.2f Mbps', ulStats(1, 1), ulStats(2, 1));
    fprintf('\nPeak DL throughput: %0.2f Mbps. Achieved average DL Throughput: %0.2f Mbps', dlStats(1, 1), dlStats(2, 1));
    fprintf('\nAchieved average UL Goodput: %0.2f Mbps. Achieved average DL Goodput: %0.2f Mbps', ulStats(5, 1), dlStats(5, 1));
    fprintf('\nPeak UL spectral efficiency: %0.2f bits/s/Hz. Achieved average UL spectral efficiency: %0.2f bits/s/Hz ', ulStats(3, 1), ulStats(4, 1));
    fprintf('\nPeak DL spectral efficiency: %0.2f bits/s/Hz. Achieved average DL spectral efficiency: %0.2f bits/s/Hz', dlStats(3, 1), dlStats(4, 1));
    disp(['Block error rate for each UE in the uplink direction: [' num2str(round(logInfo.AvgBLERLogs(:, 2)', 2)) ']']);
    disp(['Block error rate for each UE in the downlink direction: [' num2str(round(logInfo.AvgBLERLogs(:, 1)', 2)) ']']);
end

%% Simulation Logs
% The parameters used for simulation and the simulation logs are saved in MAT-files for post simulation analysis and visualization.
% The simulation parameters are saved in a MAT-file with the file name as the value of configuration parameter simParameters.ParametersLogFile. 
% For more details, see the NR Intercell Interference Modeling example.

% Get the logs
if(simParameters.DuplexMode == 0) % FDD
    logInfo = struct('NCellID', [], 'DLTimeStepLogs', [], 'ULTimeStepLogs', [], 'SchedulingAssignmentLogs', [], 'BLERLogs', [], 'AvgBLERLogs', []);
else
    logInfo = struct('TimeStepLogs', [], 'SchedulingAssignmentLogs', [], 'BLERLogs', [], 'AvgBLERLogs', []);
end
simulationLogs = cell(simParameters.NumSitesPerCluster, 1);
for siteIdx = 1:numCellsOfInterest
    logInfo.NCellID = cellsOfInterest(siteIdx);
    if(simParameters.DuplexMode == 0) % FDD
        [logInfo.DLTimeStepLogs, logInfo.ULTimeStepLogs] = getSchedulingLogs(simSchedulingLogger{siteIdx});
    else % TDD
        logInfo.TimeStepLogs = getSchedulingLogs(simSchedulingLogger{siteIdx});
    end
    logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger{siteIdx}); % Scheduling assignments log
    [logInfo.BLERLogs, logInfo.AvgBLERLogs] = getBLERLogs(simPhyLogger{siteIdx}); % Block error rate logs
    simulationLogs{siteIdx, 1} = logInfo;
end
save(simParameters.ParametersLogFile, 'simParameters'); % Save simulation parameters in a MAT-file
save(simParameters.SimulationLogFile, 'simulationLogs'); % Save simulation logs in a MAT-file
SaveFile(YUO); % Save data in YUO in a MAT-file
%MXC_2
SINR_plotting_trail(YUO);
%YXC begim
plotThroughputCDF(YUO,'DL');
%YXC end
%MXC_2
% disp('simulation complete');
% elapsed_time=toc;
% lineToPrint = ['Elapsed time: ', datestr(datenum(0,0,0,0,0,elapsed_time),'dd:HH:MM:SS')];
% disp(lineToPrint);
%MXC_2

%YXC begin
EndTime = toc;
EndTime = seconds(EndTime);
disp(['Simulation finished. Total time: ',datestr(EndTime,'dd:HH:MM:SS')])
%YXC end




