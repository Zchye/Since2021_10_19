classdef hNRUEPhy < hNRPhyInterface
    %hNRUEPhy 5G NR Phy Tx and Rx processing chains at UE
    %   The class implements the Phy Tx and Rx processing chains of 5G NR
    %   at UE. It also implements the interfaces for information exchange
    %   between Phy and higher layers. It only supports transmission of
    %   physical uplink shared channel (PUSCH) along with its demodulation
    %   reference signals (DM-RS). It only supports reception of physical
    %   downlink shared channel (PDSCH) along with its DM-RS. A single
    %   bandwidth part is assumed to cover the entire carrier bandwidth.
    %   Note that setCarrierInformation and setCellConfig need to be called
    %   on the created class object before using it.
    %
    %   hNRUEPhy methods:
    %       hNRUEPhy                - Construct a UE Phy object 
    %       run                     - Run the UE Phy layer operations
    %       setCarrierInformation   - Set the carrier configuration
    %       enablePacketLogging     - Enable packet logging
    %       registerMACInterfaceFcn - Register MAC interface functions at 
    %                                 Phy, for sending information to MAC
    %       registerInBandTxFcn     - Register callback for transmission on
    %                                 PUSCH
    %       txDataRequest           - Tx request from MAC to Phy for 
    %                                 starting PUSCH transmission
    %       dlControlRequest        - Downlink control (non-data) reception request 
    %                                 from MAC to Phy
    %       rxDataRequest           - Rx request from MAC to Phy for 
    %                                 starting PDSCH reception
    %       phyTx                   - Physical layer processing and
    %                                 transmission
    %       storeReception          - Receive the waveform and add it to the
    %                                 reception buffer
    %       phyRx                   - Physical layer reception and sending
    %                                 of decoded information to MAC layer
    %       getDLBLER               - Get the block error statistics of the
    %                                 current slot
    %
    %   Example: 
    %   % Generate a hNRUEPhy object. Configure the carrier and cell 
    %   % properties using setCarrierInformation and setCellConfig methods,
    %   % respectively
    %
    %   phyParam = struct();
    %   phyParam.SCS = 15;
    %   phyParam.NumRBs = 52;
    %   phyParam.UETxPower = 23;
    %   phyParam.DLCarrierFreq = 2.1e9;
    %   phyParam.UERxBufferSize = 1;
    %   phyParam.ChannelModelType = 'CDL';
    %   rnti = 1;
    %   phy = hNRUEPhy(phyParam, rnti);
    %
    %   carrierParam = struct();
    %   carrierParam.SubcarrierSpacing = phyParam.SCS;
    %   carrierParam.NRBsUL = phyParam.NumRBs;
    %   carrierParam.NRBsDL = phyParam.NumRBs;
    %   carrierParam.DLFreq = phyParam.DLCarrierFreq;
    %   setCarrierInformation(phy, carrierParam);
    %
    %   cellParam.NCellID = 1;
    %   cellParam.DuplexMode = 0;
    %   setCellConfig(phy, cellParam);
    %
    %   See also hNRPhyInterface
    
    %   Copyright 2020-2021 The MathWorks, Inc.
    
    properties (Access = private)
        %RNTI RNTI of the UE
        RNTI (1, 1){mustBeInRange(RNTI, 1, 65519)} = 1;
        
        %ULSCHEncoder Uplink shared channel (UL-SCH) encoder system object
        % It is an object of type nrULSCH
        ULSCHEncoder
        
        %DLSCHDecoder Downlink shared channel (DL-SCH) decoder system object
        % It is an object of type nrDLSCHDecoder
        DLSCHDecoder
        
        %WaveformInfoDL Downlink waveform information
        WaveformInfoDL
        
        %WaveformInfoUL Uplink waveform information
        WaveformInfoUL
        
        %TxAntPanel Tx antenna panel geometry
        % It is an array of the form [M, N, P, Mg, Ng] where M and N are 
        % the number of rows and columns in the antenna array, P is the 
        % number of polarizations (1 or 2), Mg and Ng are the number of row
        % and column array panels, respectively
        TxAntPanel

        %RxAntPanel Rx antenna panel geometry
        % It is an array of the form [M, N, P, Mg, Ng] where M and N are 
        % the number of rows and columns in the antenna array, P is the 
        % number of polarizations (1 or 2), Mg and Ng are the number of row
        % and column array panels, respectively
        RxAntPanel

        %TxPower Tx power in dBm
        TxPower(1, 1) {mustBeFinite, mustBeNonnegative, mustBeNonNan} = 23;
        
        %RxGain Rx antenna gain in dBi
        RxGain(1, 1) {mustBeFinite, mustBeNonnegative, mustBeNonNan} = 0;
        
        %PUSCHPDU Physical uplink shared channel (PUSCH) information sent by MAC for the current slot
        % PUSCH PDU is an object of type hNRPUSCHInfo. It has the
        % information required by Phy to transmit the MAC PDU stored in
        % object property MacPDU
        PUSCHPDU = {}
        
        %MacPDU PDU sent by MAC which is scheduled to be transmitted in the currrent slot
        % The information required to transmit this PDU is stored in object
        % property PUSCHPDU
        MacPDU = {}
        
        %CSIRSContext Rx context for the channel state information reference signals (CSI-RS)
        % This information is populated by MAC and is used by Phy to
        % receive UE's scheduled CSI-RS. It is a cell array of size 'N'
        % where N is the number of slots in a 10 ms frame. The cell
        % elements are populated with objects of type nrCSIRSConfig. An
        % element at index 'i' contains the context of CSI-RS which is sent
        % at slot index 'i-1'. Cell element at index 'i' is empty if no
        % CSI-RS reception was scheduled for the UE in the slot index 'i-1'
        CSIRSContext
        
        %CSIRSIndicationFcn Function handle to send the measured DL channel quality to MAC
        CSIRSIndicationFcn
        
        %CSIReportConfig CSI report configuration
        % The detailed explanation of this structure and its fields is
        % present as ReportConfig in <a href="matlab:help('hCQISelect')">hCQISelect</a> function
        CSIReportConfig

        %SINRvsCQI SINR to CQI mapping
        % SINRTable is a vector of 15 SINR values in dB, each corresponding to a
        % CQI value that have been computed according to the BLER condition as
        % mentioned in TS 38.214 Section 5.2.2.1
        SINRTable
        
        %DLBlkErr Downlink block error information
        % It is an array of two elements containing the number of
        % erroneously received packets and total received packets,
        % respectively
        DLBlkErr
        
        %PrevCumulativeDLBlkErr Cumulative downlink block error information returned in the last query
        % It is an array of two elements containing the number of
        % erroneously received packets and total received packets,
        % respectively
        PrevCumulativeDLBlkErr
        
        %NoiseFigure Noise figure at the receiver
        NoiseFigure = 1;
        
        %Temperature at node in Kelvin
        % It is used for thermal noise calculation
        Temperature = 300
        
        %ChannelModel Information about the propagation channel model
        % This property is an object of type nrCDLChannel if the
        % ChannelModelType is specified as 'CDL', otherwise empty
        ChannelModel
        
        %MaxChannelDelay Maximum delay introduced due to multipath components and implementation delays
        MaxChannelDelay = 0;
        
        %TimingOffset Receiver timing offset 
        TimingOffset
        
        %RxBuffer Reception buffer object to store received waveforms
        RxBuffer
        
        %SendWaveformFcn Function handle to send the waveform
        SendWaveformFcn
        
        %RankIndicator Rank advised by the UE in channel state information (CSI) report
        RankIndicator = 1;

        %PacketLogger Contains handle of the PCAP object
        PacketLogger
        
        %PacketMetaData Contains the information required for logging MAC
        %packets into PCAP
        PacketMetaData
        
        %MXC_1
        siteIdx
        
        %YXC begin
        % Yuxiao has moved this property to public properties
%         YusUtilityParameter
        %MXC_1
        %YXC end
    end
    
    %YXC begin
    properties (Access = public)
        %StoreCQIInfo stores the callback to the YusUtilityObj for storing
        %CQIInfo returned by the function hCQISelect in runtime
        StoreCQIInfo
        
        YusUtilityParameter
    end
    %YXC end
    
    methods
        function obj = hNRUEPhy(param, siteIdx, rnti)
            %MXC_1 extra constructor input
            %hNRUEPhy Construct a UE Phy object
            %   OBJ = hNRUEPHY(PARAM, RNTI) constructs a UE Phy object. It
            %   also creates the context of UL-SCH encoder system object
            %   and DL-SCH decoder system object.
            %
            %   PARAM is structure with the fields:
            %       SCS              - Subcarrier spacing
            %       UETxPower        - UE Tx Power in dBm
            %       UERxGain         - UE Rx antenna gain in dBi
            %       SINR90pc         - SINR to CQI look up table. An array of
            %                          16 SINR values correspond to 16 CQI
            %                          values (0 to 15). The look up table
            %                          contains the CQI resulting in a
            %                          maximum of 0.1 BLER for the
            %                          corresponding SINR.
            %       NumRBs           - Number of resource blocks
            %       DLCarrierFreq    - Downlink carrier frequency in Hz
            %       SubbandSize      - Size of CQI measurement sub-band in RBs
            %       UERxBufferSize   - Maximum number of waveforms to be
            %                          stored
            %       ChannelModelType - Propagation channel model type
            %       GNBTxAnts        - Number of GNB Tx antennas (Required
            %                          for CDL Tx antenna configuration)
            %       UETxAnts         - Number of Tx antennas on UEs
            %       UERxAnts         - Number of Rx antennas on UEs
            %       RankIndicator    - DL Rank to calculate precoding matrix
            %                          indicator (PMI) and channel quality 
            %                          indicator (CQI)
            %       CSIRSRowNumber   - Row number of CSI-RS resource as per TS 38.211 - Table 7.4.1.5.3-1
            %       PanelDimensions  - Antenna panel configuration as a
            %                          two-element vector in the form of [N1 N2].
            %                          N1 represents the number of antenna
            %                          elements in horizontal direction and N2
            %                          represents the number of antenna elements
            %                          in vertical direction. Valid combinations
            %                          of [N1 N2] are defined in TS 38.214 Table
            %                          5.2.2.2.12
            %       PMIMode          - PMI reporting mode. Value as 'Subband' or 'Wideband'
            %       CQIMode          - CQI reporting mode. Value as 'Subband' or 'Wideband'
            %       SubbandSize      - Subband size for CQI or PMI reporting as per TS 38.214 Table 5.2.1.4-2
            %       CodebookMode     - Codebook mode. Value as 1 or 2
            %
            %   RNTI - RNTI of the UE
            
            % Validate the subcarrier spacing
            if ~ismember(param.SCS, [15 30 60 120 240])
                error('nr5g:hNRUEPhy:InvalidSCS', 'The subcarrier spacing ( %d ) must be one of the set (15, 30, 60, 120, 240).', param.SCS);
            end
            
            obj.RNTI = rnti;
            
            %MXC_1
            obj.siteIdx = siteIdx;
            
            obj.YusUtilityParameter.Scenario = param.Scenario;
            %MXC_1
            
            obj.YusUtilityParameter.UEStat = param.UEStates{siteIdx, rnti};
            
            % Create UL-SCH encoder system object
            ulschEncoder = nrULSCH;
            ulschEncoder.MultipleHARQProcesses = true;
            obj.ULSCHEncoder = ulschEncoder;
            
            % Create DL-SCH decoder system object
            dlschDecoder = nrDLSCHDecoder;
            dlschDecoder.MultipleHARQProcesses = true;
            dlschDecoder.LDPCDecodingAlgorithm = 'Normalized min-sum';
            dlschDecoder.MaximumLDPCIterationCount = 6;
            obj.DLSCHDecoder = dlschDecoder;
            
            % Set the number of erroneous packets and total number of
            % packets received by the UE to zero
            obj.DLBlkErr = zeros(1, 2);
            obj.PrevCumulativeDLBlkErr = zeros(1, 2);
            
            obj.CSIRSContext = cell(10*(param.SCS/15)*14, 1); % Create the context for all the symbols in the frame
            % Set SINR vs CQI lookup table
            if isfield(param, 'SINR90pc')
                obj.SINRTable = param.SINR90pc;
            else
                obj.SINRTable = [-5.46 -0.46 4.54 9.05 11.54 14.04 15.54 18.04 ...
                    20.04 22.43 24.93 25.43 27.43 30.43 33.43];
            end
            
            if ~isfield(param, 'GNBTxAnts')
                param.GNBTxAnts = 1;
            % Validate the number of transmitter antennas on gNB
            elseif ~ismember(param.GNBTxAnts, [1,2,4,8,16,32,64,128,256,512,1024])
                error('nr5g:hNRUEPhy:InvalidAntennaSize',...
                    'Number of gNB Tx antennas (%d) must be a member of [1,2,4,8,16,32,64,128,256,512,1024].', param.GNBTxAnts);                                
            end            
            if ~isfield(param, 'GNBRxAnts')
                param.GNBRxAnts = 1;
            % Validate the number of receiver antennas on gNB
            elseif ~ismember(param.GNBRxAnts, [1,2,4,8,16,32,64,128,256,512,1024])
                error('nr5g:hNRUEPhy:InvalidAntennaSize',...
                    'Number of gNB Rx antennas (%d) must be a member of [1,2,4,8,16,32,64,128,256,512,1024].', param.GNBRxAnts);
            end
            if ~isfield(param, 'UETxAnts')
                param.UETxAnts = 1;
            % Validate the number of transmitter antennas on UEs
            elseif ~ismember(param.UETxAnts, [1,2,4,8,16])
                error('nr5g:hNRUEPhy:InvalidAntennaSize',...
                    'Number of UE Tx antennas (%d) must be a member of [1,2,4,8,16].', param.UETxAnts(rnti));                 
            end            
            if ~isfield(param, 'UERxAnts')
                param.UERxAnts = 1;
            % Validate the number of receiver antennas on UEs
            elseif ~ismember(param.UERxAnts, [1,2,4,8,16])              
                error('nr5g:hNRUEPhy:InvalidAntennaSize',...
                    'Number of UE Rx antennas (%d) must be a member of [1,2,4,8,16].', param.UERxAnts(rnti));
            end
            % Validate the rank indicator
            if isfield(param, 'RankIndicator')
                validateattributes(param.RankIndicator, {'numeric'}, ...
                    {'nonempty', 'integer', 'scalar', '>=', 1, '<=', min(param.GNBTxAnts, param.UERxAnts)},...
                    'param.RankIndicator', 'RankIndicator');
                obj.RankIndicator = param.RankIndicator;
            end
            
            
            %MXC_2
            %{
            % Downlink direction i.e. gNB's Tx antenna panel
            % configuration and UE's Rx panel configuration
            [gNBTxAntPanel, obj.RxAntPanel] = ...
                hArrayGeometry(param.GNBTxAnts, param.UERxAnts, 'downlink');
            % Uplink direction i.e. UE's Tx antenna panel
            % configuration and gNB's Rx panel configuration
            [obj.TxAntPanel, ~] = ...
                hArrayGeometry( param.UETxAnts, param.GNBRxAnts, 'uplink');
            %}
            %MXC_2
            
            

            % Set CSI report config
            obj.CSIReportConfig.NStartBWP = 0;
            obj.CSIReportConfig.NSizeBWP = param.NumRBs;
            % Set panel dimensions:[N1, N2]
            csirsConfig = nrCSIRSConfig;
            if isfield(param, 'CSIRSRowNumber')
                csirsConfig.RowNumber = param.CSIRSRowNumber;
            else
                % Possible CSI-RS resource row numbers for single transmit antenna case are 1 and 2
                csirsConfig.RowNumber = 2;
            end
            % PanelDimension property is applicable only to 2 or more GNBTxAnts antennas
            if param.GNBTxAnts > 1
                    if csirsConfig.NumCSIRSPorts ~= 2*prod(param.PanelDimensions)
                        error('nr5g:hNRUEPhy:InvalidPanelDimension',...
                            'Number of CSI-RS ports (%d) must be equal to the number of CSI-RS antenna elements (%d)',...
                            csirsConfig.NumCSIRSPorts, 2*prod(param.PanelDimensions));
                    end
                    obj.CSIReportConfig.PanelDimensions = param.PanelDimensions;
            end
            % Set CQI mode
            if isfield(param, 'CQIMode')
                if ~ismember(param.CQIMode, {'Wideband', 'Subband'})
                    error('nr5g:hNRUEPhy:InvalidCQIMode', "CQIMode must be either 'Wideband' or 'Subband'.");
                end
                obj.CSIReportConfig.CQIMode =  param.CQIMode;
            else
                obj.CSIReportConfig.CQIMode = 'Subband';
            end
            % Set precoding matrix indicator (PMI) mode
            if isfield(param, 'PMIMode')
                if ~ismember(param.PMIMode, {'Wideband', 'Subband'})
                    error('nr5g:hNRUEPhy:InvalidPMIMode', "PMIMode must be either 'Wideband' or 'Subband'.");
                end
                obj.CSIReportConfig.PMIMode =  param.PMIMode;
            else
                obj.CSIReportConfig.PMIMode = 'Subband';
            end
            % Set subband size
            if isfield(param, 'SubbandSize')
                obj.CSIReportConfig.SubbandSize =  param.SubbandSize;
            else
                obj.CSIReportConfig.SubbandSize = 4;
            end
            % Set precoding resource block group (PRG) size
            if isfield(param, 'PRGSize')
                obj.CSIReportConfig.PRGSize =  param.PRGSize;
            else
                obj.CSIReportConfig.PRGSize = [];
            end
            % Set codebook mode
            if isfield(param, 'CodebookMode')
                validateattributes(param.CodebookMode, 'numeric', {'scalar', 'integer', '>=', 1, '<=', 2}, 'param.CodebookMode', 'CodebookMode')
                obj.CSIReportConfig.CodebookMode =  param.CodebookMode;
            else
                obj.CSIReportConfig.CodebookMode = 1;
            end
            
            %Initialize DMRS and CSI-RS timing offsets
            obj.TimingOffset = 0;
            
            % Set Tx power in dBm
            if isfield(param, 'UETxPower')
                obj.TxPower = param.UETxPower;
            end
            
            % Set Rx antenna gain in dBi
            if isfield(param, 'UERxGain')
                obj.RxGain = param.UERxGain;
            end
            
            waveformInfo = nrOFDMInfo(param.NumRBs, param.SCS);

            % Create channel model object
            if isfield(param, 'ChannelModelType')
                if strcmpi(param.ChannelModelType, 'CDL')
                    %MXC_1
                    %{
                    obj.ChannelModel = nrCDLChannel; % CDL channel object
                    obj.ChannelModel.DelayProfile = 'CDL-C';
                    obj.ChannelModel.DelaySpread = 300e-9;
                    % Set the carrier frequency to downlink
                    obj.ChannelModel.CarrierFrequency = param.DLCarrierFreq;
                    % Size of antenna array [M N P Mg Ng], where:
                    %    - M and N are the number of rows and columns in
                    %      the antenna array, respectively.
                    %    - P is the number of polarizations (1 or 2).
                    %    - Mg and Ng are the number of row and column
                    %      array panels, respectively.
                    obj.ChannelModel.TransmitAntennaArray.Size = gNBTxAntPanel;
                    obj.ChannelModel.ReceiveAntennaArray.Size = obj.RxAntPanel;
                    obj.ChannelModel.SampleRate = waveformInfo.SampleRate;
                    chInfo = info(obj.ChannelModel);
                    % Update the maximum delay caused due to CDL channel model
                    obj.MaxChannelDelay = ceil(max(chInfo.PathDelays*obj.ChannelModel.SampleRate)) + chInfo.ChannelFilterDelay;
                    %}
                    obj.ChannelModel = nrCDLChannel; % CDL channel object
                    obj.ChannelModel.DelayProfile = 'Custom';
                    
                    LinkDir = 0; % DL
                    CHParam = ConfigChannel(param, LinkDir, obj.siteIdx, obj.RNTI);
                    
                    
                    obj.ChannelModel.CarrierFrequency = param.DLCarrierFreq;
                    obj.ChannelModel.HasLOSCluster = logical(CHParam.LOS);
                    obj.ChannelModel.AngleSpreads = [CHParam.LSPs.ASD, CHParam.LSPs.ASA, CHParam.LSPs.ZSD, CHParam.LSPs.ZSA];
                    if CHParam.LOS
                        obj.ChannelModel.KFactorFirstCluster = CHParam.LSPs.K;
                    end
                    obj.ChannelModel.PathDelays = CHParam.tau_n';
                    obj.ChannelModel.AveragePathGains = CHParam.P_n';
                    obj.ChannelModel.AnglesAoA = CHParam.ADAngles.phi_n_AOA';
                    obj.ChannelModel.AnglesAoD = CHParam.ADAngles.phi_n_AOD';
                    obj.ChannelModel.AnglesZoA = CHParam.ADAngles.theta_n_ZOA';
                    obj.ChannelModel.AnglesZoD = CHParam.ADAngles.theta_n_ZOD';
                    obj.ChannelModel.XPR = CHParam.XPR;
                    %MXC_1 quick fix
                    %causes problems with 100M element array, change to zero for now
                    obj.ChannelModel.NumStrongestClusters = 2;
%                     obj.ChannelModel.NumStrongestClusters = 0;
                    %Mxc_1
                    %not relevant anymore due to quick fix
                    obj.ChannelModel.ClusterDelaySpread = CHParam.c_DS;
                    %MXC_1
                    
                    %MXC_2
                    obj.TxAntPanel = param.UETxAntPanelSize;
                    obj.RxAntPanel = param.UERxAntPanelSize;
                    
                    %gNB Antenna Config
                    obj.ChannelModel.TransmitAntennaArray.Size = param.GNBTxAntPanelSize;
                    obj.ChannelModel.TransmitAntennaArray.ElementSpacing = param.GNBTxAntElementSpacing;
                    obj.ChannelModel.TransmitAntennaArray.PolarizationAngles = param.GNBTxAntPolarizationAngles;
                    obj.ChannelModel.TransmitAntennaArray.Element = param.GNBAntElement;
                    obj.ChannelModel.TransmitAntennaArray.PolarizationModel = param.GNBAntPolarizationModel;
                    obj.ChannelModel.TransmitArrayOrientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
                    
                    %UE Antenna Config
                    obj.ChannelModel.ReceiveAntennaArray.Size = obj.RxAntPanel;
                    obj.ChannelModel.ReceiveAntennaArray.ElementSpacing = param.UERxAntElementSpacing;
                    obj.ChannelModel.ReceiveAntennaArray.PolarizationAngles = param.UERxAntPolarizationAngles;
                    obj.ChannelModel.ReceiveAntennaArray.Element = param.UEAntElement;
                    obj.ChannelModel.ReceiveAntennaArray.PolarizationModel = param.UEAntPolarizationModel;
                    obj.ChannelModel.ReceiveArrayOrientation = [0; 0; 0];
                    
                    % UE Mobility
                    c = physconst('lightspeed'); % speed of light in m/s
                    fd = param.UEStates{obj.siteIdx,obj.RNTI}.Speed/c*param.DLCarrierFreq; % UE max Doppler frequency in Hz
                    obj.ChannelModel.MaximumDopplerShift = fd;
                    obj.ChannelModel.UTDirectionOfTravel = [param.UEStates{obj.siteIdx,obj.RNTI}.DirectionOfTravel; 90];
                    %MXC_2
                    
                    obj.ChannelModel.SampleRate = waveformInfo.SampleRate;
                    chInfo = info(obj.ChannelModel);
                    % Update the maximum delay caused due to CDL channel model
                    obj.MaxChannelDelay = ceil(max(chInfo.PathDelays*obj.ChannelModel.SampleRate)) + chInfo.ChannelFilterDelay;
                    %}
                end
            end
            
            % Set receiver noise figure
            if isfield(param, 'NoiseFigure')
                obj.NoiseFigure = param.NoiseFigure;
            end
            
            % Create reception buffer object
            if isfield(param, 'UERxBufferSize')
                obj.RxBuffer = hNRPhyRxBuffer('BufferSize', param.UERxBufferSize, 'NRxAnts', param.UERxAnts);
            else
                obj.RxBuffer = hNRPhyRxBuffer('NRxAnts', param.UERxAnts);
            end
        end
        
        function run(obj)
            %run Run the UE Phy layer operations
            
            % Phy processing and transmission of PUSCH (along with its DM-RS).
            % It assumes that MAC has already loaded the Phy Tx
            % context for anything scheduled to be transmitted at the
            % current symbol
            phyTx(obj);
            
            % Phy reception of PDSCH (along with its DM-RS) and CSI-RS, and then sending the decoded information to MAC.
            % PDSCH Rx is done in the symbol after the last symbol in PDSCH
            % duration (till then the packets are queued in Rx buffer). Phy
            % calculates the last symbol of PDSCH duration based on
            % 'rxDataRequest' call from MAC (which comes at the first
            % symbol of PDSCH Rx time) and the PDSCH duration. CSI-RS
            % reception is done at the start of slot which is after the
            % scheduled CSI-RS reception slot
            phyRx(obj);
        end
        
        function setCarrierInformation(obj, carrierInformation)
            %setCarrierInformation Set the carrier configuration
            %   setCarrierInformation(OBJ, CARRIERINFORMATION) sets the
            %   carrier configuration, CARRIERINFORMATION.
            %   CARRIERINFORMATION is a structure including the following
            %   fields:
            %       SubcarrierSpacing  - Sub carrier spacing used. Assuming
            %                            single bandwidth part in the whole
            %                            carrier
            %       NRBsDL             - Downlink bandwidth in terms of
            %                            number of resource blocks
            %       NRBsUL             - Uplink bandwidth in terms of
            %                            number of resource blocks
            %       DLBandwidth        - Downlink bandwidth in Hz
            %       ULBandwidth        - Uplink bandwidth in Hz
            %       DLFreq             - Downlink carrier frequency in Hz
            %       ULFreq             - Uplink carrier frequency in Hz
            
            setCarrierInformation@hNRPhyInterface(obj, carrierInformation);
            
            % Initialize data Rx context
            obj.DataRxContext = cell(obj.CarrierInformation.SymbolsPerFrame, 1);
            % Set waveform properties
            setWaveformProperties(obj, obj.CarrierInformation);
        end
        
        function timestamp = getCurrentTime(obj)
            %getCurrentTime Return the current timestamp of node in microseconds
            
            % Calculate number of samples from the start of the current
            % frame to the current symbol
            numSubFrames = floor(obj.CurrSlot / obj.WaveformInfoDL.SlotsPerSubframe);
            numSlotSubFrame = mod(obj.CurrSlot, obj.WaveformInfoDL.SlotsPerSubframe);
            symbolNumSubFrame = numSlotSubFrame*obj.WaveformInfoDL.SymbolsPerSlot + obj.CurrSymbol;
            numSamples = (numSubFrames * sum(obj.WaveformInfoDL.SymbolLengths))...
                + sum(obj.WaveformInfoDL.SymbolLengths(1:symbolNumSubFrame));
            
            % Timestamp in microseconds
            timestamp = (obj.AFN * 0.01) + (numSamples *  1 / obj.WaveformInfoDL.SampleRate);
            timestamp = (1e6 * timestamp);
        end
        
        function enablePacketLogging(obj, fileName)
            %enablePacketLogging Enable packet logging
            %
            % FILENAME - Name of the PCAP file
            
            % Create packet logging object
            obj.PacketLogger = hNRPacketWriter('FileName', fileName, 'FileExtension', 'pcap');
            obj.PacketMetaData = hNRPacketInfo;
            if obj.CellConfig.DuplexMode % Radio type
                obj.PacketMetaData.RadioType = obj.PacketLogger.RadioTDD;
            else
                obj.PacketMetaData.RadioType = obj.PacketLogger.RadioFDD;
            end
            obj.PacketMetaData.RNTIType = obj.PacketLogger.CellRNTI;
            obj.PacketMetaData.RNTI = obj.RNTI;
        end
        
        function registerMACInterfaceFcn(obj, sendMACPDUFcn, sendDLChannelQualityFcn)
            %registerMACInterfaceFcn Register MAC interface functions at Phy, for sending information to MAC
            %   registerMACInterfaceFcn(OBJ, SENDMACPDUFCN,
            %   SENDDLCHANNELQUALITYFCN) registers the callback function to
            %   send decoded MAC PDUs and measured DL channel quality to MAC.
            %
            %   SENDMACPDUFCN Function handle provided by MAC to Phy for
            %   sending PDUs to MAC.
            %
            %   SENDDLCHANNELQUALITYFCN Function handle provided by MAC to Phy for
            %   sending the measured DL channel quality (measured on CSI-RS).
            
            obj.RxIndicationFcn = sendMACPDUFcn;
            obj.CSIRSIndicationFcn = sendDLChannelQualityFcn;
        end
        
        function registerInBandTxFcn(obj, sendWaveformFcn)
            %registerInBandTxFcn Register callback for transmission on PUSCH
            %
            %   SENDWAVEFORMFCN Function handle provided by packet
            %   distribution object for sending packets to nodes.
            
            obj.SendWaveformFcn = sendWaveformFcn;
        end
        
        function txDataRequest(obj, PUSCHInfo, macPDU)
            %txDataRequest Tx request from MAC to Phy for starting PUSCH transmission
            %  txDataRequest(OBJ, PUSCHINFO, MACPDU) sets the Tx context to
            %  indicate PUSCH transmission in the current symbol
            %
            %  PUSCHInfo is an object of type hNRPUSCHInfo sent by MAC. It
            %  contains the information required by the Phy for the
            %  transmission.
            %
            %  MACPDU is the uplink MAC PDU sent by MAC for transmission.
            
            obj.PUSCHPDU = PUSCHInfo;
            obj.MacPDU = macPDU;
        end
        
        
        function dlControlRequest(obj, pduType, dlControlPDU)
            %dlControlRequest Downlink control (non-data) reception request from MAC to Phy
            %   dlControlRequest(OBJ, PDUTYPES, DLCONTROLPDUS) is a request
            %   from MAC for downlink receptions. MAC sends it at the start
            %   of a DL slot for all the scheduled non-data DL receptions
            %   in the slot (Data i.e. PDSCH reception information is sent
            %   by MAC using rxDataRequest interface of this class).
            %
            %   PDUTYPE is an array of packet types. Currently, only
            %   packet type 0 (CSI-RS) is supported.
            %
            %   DLCONTROLPDU is an array of DL control information PDUs,
            %   corresponding to packet types in PDUTYPE. Currently
            %   supported information CSI-RS PDU is an object of type
            %   nrCSIRSConfig.
            
            % Update the Rx context for DL receptions
            for i = 1:length(pduType)
                switch(pduType(i))
                    case obj.CSIRSPDUType
                        % CSI-RS would be read at the start of next slot
                        nextSlot = mod(obj.CurrSlot+1, obj.CarrierInformation.SlotsPerSubframe*10);
                        obj.CSIRSContext{nextSlot+1} = dlControlPDU{i};
                end
            end
        end
        
        function rxDataRequest(obj, pdschInfo)
            %rxDataRequest Rx request from MAC to Phy for starting PDSCH reception
            %   rxDataRequest(OBJ, PDSCHINFO) is a request to start PDSCH
            %   reception. It starts a timer for PDSCH end time (which on
            %   triggering receives the complete PDSCH). The Phy expects
            %   the MAC to send this request at the start of reception
            %   time.
            %
            %   PDSCHInfo is an object of type hNRPDSCHInfo. It contains the
            %   information required by the Phy for the reception.
            
            symbolNumFrame = obj.CurrSlot*14 + obj.CurrSymbol; % Current symbol number w.r.t start of 10 ms frame
            
            % PDSCH to be read in the symbol after the last symbol in
            % PDSCH reception
            numPDSCHSym =  pdschInfo.PDSCHConfig.SymbolAllocation(2);
            pdschRxSymbolFrame = mod(symbolNumFrame + numPDSCHSym, obj.CarrierInformation.SymbolsPerFrame);
            
            % Add the PDSCH RX information at the index corresponding to
            % the symbol just after PDSCH end time
            obj.DataRxContext{pdschRxSymbolFrame+1} = pdschInfo;
        end
        
        function phyTx(obj)
            %phyTx Physical layer processing and transmission
            
            if ~isempty(obj.PUSCHPDU) % If any PUSCH is scheduled for current symbol
                % Initialize Tx grid
                numTxAnts = prod(obj.TxAntPanel);
                txSlotGrid = zeros(obj.CarrierInformation.NRBsUL*12, obj.WaveformInfoUL.SymbolsPerSlot, numTxAnts);
                
                % Fill PUSCH symbols in the grid
                txSlotGrid = populatePUSCH(obj, obj.PUSCHPDU, obj.MacPDU, txSlotGrid);
                
                % OFDM modulation
                carrier = nrCarrierConfig;
                carrier.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
                carrier.NSizeGrid = obj.CarrierInformation.NRBsDL;
                carrier.NSlot = obj.CurrSlot;
                txWaveform = nrOFDMModulate(carrier, txSlotGrid);
                
                % Trim txWaveform to span only the transmission symbols
                numTxSymbols = obj.PUSCHPDU.PUSCHConfig.SymbolAllocation(2);
                slotNumSubFrame = mod(obj.CurrSlot, obj.WaveformInfoUL.SlotsPerSubframe);
                startSymSubframe = slotNumSubFrame*obj.WaveformInfoUL.SymbolsPerSlot + 1; % Start symbol of current slot in the subframe
                lastSymSubframe = startSymSubframe + obj.WaveformInfoUL.SymbolsPerSlot -1; % Last symbol of current slot in the subframe
                symbolLengths = obj.WaveformInfoUL.SymbolLengths(startSymSubframe : lastSymSubframe); % Length of symbols of current slot
                startSample = sum(symbolLengths(1:obj.CurrSymbol)) + 1;
                endSample = sum(symbolLengths(1:obj.CurrSymbol+numTxSymbols));
                txWaveform = txWaveform(startSample:endSample);
                
                % Apply Tx power and gain
                gain = 0;
                txWaveform = applyTxPowerLevelAndGain(obj, txWaveform, gain);
                
                % Construct packet information
                packetInfo.Waveform = txWaveform;
                packetInfo.RNTI = obj.RNTI;
                packetInfo.Position = obj.Node.NodePosition;
                packetInfo.CarrierFreq = obj.CarrierInformation.ULFreq;
                packetInfo.TxPower = obj.TxPower;
                packetInfo.NTxAnts = numTxAnts;
                packetInfo.SampleRate = obj.WaveformInfoUL.SampleRate;
                
                % Waveform transmission by sending it to packet
                % distribution entity
                obj.SendWaveformFcn(packetInfo);
            end
            
            % Clear the Tx contexts
            obj.PUSCHPDU = {};
            obj.MacPDU = {};
        end
        
        function storeReception(obj, waveformInfo)
            %storeReception Receive the waveform and add it to the reception
            % buffer
            
            % Apply channel model
            rxWaveform = applyChannelModel(obj, waveformInfo);
            currTime = getCurrentTime(obj);
            rxWaveformInfo = struct('Waveform', rxWaveform, ...
                'NumSamples', size(rxWaveform, 1), ...
                'SampleRate', waveformInfo.SampleRate, ...
                'StartTime', currTime);
            
            % Store the received waveform in the buffer
            addWaveform(obj.RxBuffer, rxWaveformInfo);
        end
        
        function phyRx(obj)
            %phyRx Physical layer reception and sending of decoded information to MAC layer
            
            symbolNumFrame = obj.CurrSlot*14 + obj.CurrSymbol; % Current symbol number w.r.t start of 10 ms frame
            pdschInfo = obj.DataRxContext{symbolNumFrame + 1};
            csirsInfo = obj.CSIRSContext{obj.CurrSlot + 1};
            
            if isempty(pdschInfo) && isempty(csirsInfo)
                return; % No packet is scheduled to be read at the current symbol
            end
            
            currentTime = getCurrentTime(obj);
            duration = 0;
               
            % Calculate the reception duration
            if ~isempty(pdschInfo)
                startSymPDSCH = pdschInfo.PDSCHConfig.SymbolAllocation(1);
                numSymPDSCH = pdschInfo.PDSCHConfig.SymbolAllocation(2);
                % Calculate the symbol start index w.r.t start of 1 ms sub frame
                slotNumSubFrame = mod(pdschInfo.NSlot, obj.WaveformInfoDL.SlotsPerSubframe);
                pdschSymbolSet = startSymPDSCH : startSymPDSCH+numSymPDSCH-1;
                symbolSetSubFrame = (slotNumSubFrame * 14) + pdschSymbolSet + 1;
                duration = 1e6 * (1/obj.WaveformInfoDL.SampleRate) * ...
                    sum(obj.WaveformInfoDL.SymbolLengths(symbolSetSubFrame)); % In microseconds
            end
            
            if ~isempty(csirsInfo)
                % The CSI-RS which is currently being read was sent in the
                % last slot. Read the complete last slot
                %duration = 1e6 * (1e-3 / obj.WaveformInfoDL.SlotsPerSubframe);  % In microseconds
                
                % Calculate the symbol start index w.r.t start of 1 ms sub frame
                if obj.CurrSlot > 0
                    txSlot = obj.CurrSlot-1;
                else
                    txSlot = obj.WaveformInfoDL.SlotsPerSubframe*10-1;
                end
                slotNumSubFrame = mod(txSlot, obj.WaveformInfoDL.SlotsPerSubframe);
                symbolSet = 0:13;
                symbolSetSubFrame = (slotNumSubFrame * 14) + symbolSet + 1;
                duration = 1e6 * (1/obj.WaveformInfoDL.SampleRate) * ...
                    sum(obj.WaveformInfoDL.SymbolLengths(symbolSetSubFrame)); % In microseconds
                
            end
           
            % Convert channel delay into microseconds
            maxChannelDelay = 1e6 * (1/obj.WaveformInfoDL.SampleRate) * obj.MaxChannelDelay;
            
            % Get the received waveform
            duration = duration + maxChannelDelay;
            rxWaveform = getReceivedWaveform(obj.RxBuffer, currentTime + maxChannelDelay - duration, duration, obj.WaveformInfoDL.SampleRate);
            
            % Process the waveform and send the decoded information to MAC
            phyRxProcessing(obj, rxWaveform, pdschInfo, csirsInfo);
           
            % Clear the Rx contexts
            obj.DataRxContext{symbolNumFrame + 1} = {};
            obj.CSIRSContext{obj.CurrSlot + 1} = {};
        end
 
        function dlBLER = getDLBLER(obj)
            %getDLBLER Get the block error statistics of the current slot
            
            dlBLER = obj.DLBlkErr - obj.PrevCumulativeDLBlkErr;
            obj.PrevCumulativeDLBlkErr = obj.DLBlkErr;
        end
    end
    
    methods (Access = private)
        function setWaveformProperties(obj, carrierInformation)
            %setWaveformProperties Set the UL and DL waveform properties
            %   setWaveformProperties(OBJ, CARRIERINFORMATION) sets the UL
            %   and DL waveform properties ae per the information in
            %   CARRIERINFORMATION. CARRIERINFORMATION is a structure
            %   including the following fields:
            %       SubcarrierSpacing  - Subcarrier spacing used
            %       NRBsDL             - Downlink bandwidth in terms of number of resource blocks
            %       NRBsUL             - Uplink bandwidth in terms of number of resource blocks
            
            % Set the UL waveform properties
            obj.WaveformInfoUL = nrOFDMInfo(carrierInformation.NRBsUL, carrierInformation.SubcarrierSpacing);
            
            % Set the DL waveform properties
            obj.WaveformInfoDL = nrOFDMInfo(carrierInformation.NRBsDL, carrierInformation.SubcarrierSpacing);
        end
        
        function updatedSlotGrid = populatePUSCH(obj, puschInfo, macPDU, txSlotGrid)
            %populatePUSCH Populate PUSCH symbols in the Tx grid and return the updated grid
            
            % Set transport block in the encoder. In case of empty MAC PDU
            % sent from MAC (indicating retransmission), no need to set
            % transport block as it is already buffered in UL-SCH encoder
            % object
            if ~isempty(macPDU)
                % A non-empty MAC PDU is sent by MAC which indicates new
                % transmission
                macPDUBitmap = de2bi(macPDU, 8);
                macPDUBitmap = reshape(macPDUBitmap', [], 1); % Convert to column vector
                setTransportBlock(obj.ULSCHEncoder, macPDUBitmap, puschInfo.HARQID);
            end
            
            if ~isempty(obj.PacketLogger) % Packet capture enabled
                % Log uplink packets
                if isempty(macPDU)
                    tbID = 0; % Transport block id
                    macPDUBitmap = getTransportBlock(obj.ULSCHEncoder, tbID, puschInfo.HARQID);
                    macPDUBitmap = reshape(macPDUBitmap, 8, [])';
                    macPDU = bi2de(macPDUBitmap);
                end
                logPackets(obj, puschInfo, macPDU, 1)
            end
            
            % Calculate PUSCH and DM-RS information
            carrierConfig = nrCarrierConfig;
            carrierConfig.NSizeGrid = obj.CarrierInformation.NRBsUL;
            carrierConfig.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
            carrierConfig.NSlot = puschInfo.NSlot;
            carrierConfig.NCellID = obj.CellConfig.NCellID;
            [puschIndices, puschIndicesInfo] = nrPUSCHIndices(carrierConfig, puschInfo.PUSCHConfig);
            dmrsSymbols = nrPUSCHDMRS(carrierConfig, puschInfo.PUSCHConfig);
            dmrsIndices = nrPUSCHDMRSIndices(carrierConfig, puschInfo.PUSCHConfig);
            
            % UL-SCH encoding
            obj.ULSCHEncoder.TargetCodeRate = puschInfo.TargetCodeRate;
            codedTrBlock = obj.ULSCHEncoder(puschInfo.PUSCHConfig.Modulation, puschInfo.PUSCHConfig.NumLayers, ...
                puschIndicesInfo.G, puschInfo.RV, puschInfo.HARQID);
            
            % PUSCH modulation
            puschSymbols = nrPUSCH(carrierConfig, puschInfo.PUSCHConfig, codedTrBlock);
            
            % PUSCH mapping in the grid
            [~,puschAntIndices] = nrExtractResources(puschIndices,txSlotGrid);
            txSlotGrid(puschAntIndices) = puschSymbols;
            
            % PUSCH DM-RS precoding and mapping
            [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,1),txSlotGrid);
            txSlotGrid(dmrsAntIndices) = txSlotGrid(dmrsAntIndices) + dmrsSymbols(:,1);
            
            updatedSlotGrid = txSlotGrid;
        end
        
        function rxWaveform = applyChannelModel(obj, pktInfo)
            %applyChannelModel Return the waveform after applying channel model
            
            rxWaveform = pktInfo.Waveform;
            % Check if channel model is specified between gNB and UE
            if ~isempty(obj.ChannelModel)
                rxWaveform = [rxWaveform; zeros(obj.MaxChannelDelay, size(rxWaveform,2))];
                rxWaveform = obj.ChannelModel(rxWaveform);
            end
            
            % Apply path loss on the waveform
            [rxWaveform, pathloss] = applyPathLoss(obj, rxWaveform, pktInfo);
            pktInfo.TxPower = pktInfo.TxPower - pathloss;
            
            % Apply receiver antenna gain
            rxWaveform = applyRxGain(obj, rxWaveform);
            pktInfo.TxPower = pktInfo.TxPower + obj.RxGain;
            
            % Add thermal noise to the waveform
            selfInfo.Temperature = obj.Temperature;
            selfInfo.Bandwidth = obj.CarrierInformation.DLBandwidth;
            rxWaveform = applyThermalNoise(obj, rxWaveform, pktInfo, selfInfo);
        end
        
        function phyRxProcessing(obj, rxWaveform, pdschInfo, csirsInfo)
            %phyRxProcessing Process the waveform and send the decoded information to MAC
            
            carrier = nrCarrierConfig;
            carrier.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
            carrier.NSizeGrid = obj.CarrierInformation.NRBsDL;
            
            % Get the Tx slot
            if obj.CurrSymbol ==0 % Current symbol is first in the slot hence transmission was done in the last slot
                if obj.CurrSlot > 0
                    txSlot = obj.CurrSlot-1;
                    txSlotAFN = obj.AFN; % Tx slot was in the current frame
                else
                    txSlot = obj.WaveformInfoDL.SlotsPerSubframe*10-1;
                    % Tx slot was in the previous frame
                    txSlotAFN = obj.AFN - 1;
                end
                lastSym = obj.WaveformInfoDL.SymbolsPerSlot-1; % Last symbol number of the waveform
            else % Transmission was done in the current slot
                txSlot = obj.CurrSlot;
                txSlotAFN = obj.AFN; % Tx slot was in the current frame
                lastSym = obj.CurrSymbol - 1; % Last symbol number of the waveform
            end
            if ~isempty(csirsInfo)
                startSym = 0; % Read full slot
            else
                % Read from PDSCH start symbol
                startSym = lastSym - pdschInfo.PDSCHConfig.SymbolAllocation(2) + 1;
            end
            
            carrier.NSlot = txSlot;
            carrier.NFrame = txSlotAFN;
            
            % Populate the the received waveform at appropriate indices in the slot-length waveform
            slotNumSubFrame = mod(txSlot, obj.WaveformInfoDL.SlotsPerSubframe);
            startSampleIndexSlot = slotNumSubFrame*obj.WaveformInfoDL.SymbolsPerSlot + 1; % Start sample index of tx slot w.r.t start of subframe
            endSampleIndexSlot = startSampleIndexSlot + obj.WaveformInfoDL.SymbolsPerSlot -1; % End sample index of tx slot w.r.t start of subframe
            slotWaveform = zeros(sum(obj.WaveformInfoDL.SymbolLengths(startSampleIndexSlot:endSampleIndexSlot)) + obj.MaxChannelDelay, prod(obj.RxAntPanel));
            startIndex = sum(obj.WaveformInfoDL.SymbolLengths(1 : startSym))+1;  % Calculate the symbol start index w.r.t start of 1 ms subframe
            slotWaveform(startIndex : startIndex+length(rxWaveform)-1, :) = rxWaveform;

            % Calculate timing offset
            if ~isempty(pdschInfo) % If PDSCH is scheduled to be received in the waveform
                % Calculate PDSCH and DM-RS information
                carrier.NCellID = obj.CellConfig.NCellID;
                [pdschIndices, ~] = nrPDSCHIndices(carrier, pdschInfo.PDSCHConfig);
                dmrsSymbols = nrPDSCHDMRS(carrier, pdschInfo.PDSCHConfig);
                dmrsIndices = nrPDSCHDMRSIndices(carrier, pdschInfo.PDSCHConfig);
                % Calculate timing offset
                [t,mag] = nrTimingEstimate(carrier, slotWaveform, dmrsIndices, dmrsSymbols);
                obj.TimingOffset = hSkipWeakTimingOffset(obj.TimingOffset, t, mag);
            else
                % If only CSI-RS is present in the waveform
                csirsInd = nrCSIRSIndices(carrier, csirsInfo);
                csirsSym = nrCSIRS(carrier, csirsInfo);
                % Calculate timing offset
                [t,mag] = nrTimingEstimate(carrier, slotWaveform, csirsInd, csirsSym);
                obj.TimingOffset = hSkipWeakTimingOffset(obj.TimingOffset, t, mag);
            end
            
            if obj.TimingOffset > obj.MaxChannelDelay
                % Ignore the timing offset estimate resulting from weak correlation 
                obj.TimingOffset = 0;
            end
            slotWaveform = slotWaveform(1+obj.TimingOffset:end, :);
            % Perform OFDM demodulation on the received data to recreate the
            % resource grid, including padding in the event that practical
            % synchronization results in an incomplete slot being demodulated
            rxGrid = nrOFDMDemodulate(carrier, slotWaveform);
            [K,L,R] = size(rxGrid);
            if (L < obj.WaveformInfoDL.SymbolsPerSlot)
                rxGrid = cat(2,rxGrid,zeros(K, obj.WaveformInfoDL.SymbolsPerSlot-L, R));
            end
            
            % Decode MAC PDU if PDSCH is present in waveform
            if ~isempty(pdschInfo)
                obj.DLSCHDecoder.TransportBlockLength = pdschInfo.TBS*8;
                obj.DLSCHDecoder.TargetCodeRate = pdschInfo.TargetCodeRate;
                [estChannelGrid,noiseEst] = nrChannelEstimate(rxGrid, dmrsIndices, dmrsSymbols, 'CDMLengths', pdschInfo.PDSCHConfig.DMRS.CDMLengths);
                % Get PDSCH resource elements from the received grid
                [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGrid);
                
                % Equalization
                [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
                
                % PDSCH decoding
                [dlschLLRs,rxSymbols] = nrPDSCHDecode(pdschEq, pdschInfo.PDSCHConfig.Modulation, pdschInfo.PDSCHConfig.NID, ...
                    pdschInfo.PDSCHConfig.RNTI, noiseEst);
                
                % Scale LLRs by CSI
                csi = nrLayerDemap(csi); % CSI layer demapping
                
                cwIdx = 1;
                Qm = length(dlschLLRs{1})/length(rxSymbols{cwIdx}); % bits per symbol
                csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);   % expand by each bit per symbol
                dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale
                
                %YXC begin
                % Store SINR calculated from DMRS in YUO
                L_csi = length(csi);
                mean_csi = zeros(L_csi,1);
                for ii = 1:L_csi
                    mean_csi(ii) = mean(csi{ii}(:));
                end
                DMRSSINR = mean(mean_csi)/noiseEst-1;   % SINR in linear scale
                obj.YusUtilityParameter.YUO.storeDMRSSINR(DMRSSINR);
                %YXC end
                
                [decbits, crcFlag] = obj.DLSCHDecoder(dlschLLRs, pdschInfo.PDSCHConfig.Modulation, ...
                    pdschInfo.PDSCHConfig.NumLayers, pdschInfo.RV, pdschInfo.HARQID);
                
                if pdschInfo.RV == 1
                    % The last redundancy version as per the order [0 3 2
                    % 1] failed. Reset the soft buffer
                    resetSoftBuffer(obj.DLSCHDecoder, 0, pdschInfo.HARQID);
                end
                
                % Convert bit stream to byte stream
                decbits = (reshape(decbits, 8, []))';
                macPDU = bi2de(decbits);
                
                % Rx callback to MAC
                macPDUInfo = hNRRxIndicationInfo;
                macPDUInfo.RNTI = pdschInfo.PDSCHConfig.RNTI;
                macPDUInfo.TBS = pdschInfo.TBS;
                macPDUInfo.HARQID = pdschInfo.HARQID;
                obj.RxIndicationFcn(macPDU, crcFlag, macPDUInfo); % Send PDU to MAC
                
                % Increment the number of erroneous packets
                obj.DLBlkErr(1) = obj.DLBlkErr(1) + crcFlag;
                % Increment the total number of received packets
                obj.DLBlkErr(2) = obj.DLBlkErr(2) + 1;
                
                if ~isempty(obj.PacketLogger) % Packet capture enabled
                    logPackets(obj, pdschInfo, macPDU, 0); % Log DL packets
                end
            end
            
            % If CSI-RS is present in waveform, measure DL channel quality
            if ~isempty(csirsInfo)
                csirsSym = nrCSIRS(carrier, csirsInfo);
                csirsRefInd = nrCSIRSIndices(carrier, csirsInfo);
                if ~isempty(csirsRefInd)
                    cdmType = csirsInfo.CDMType;
                    if ~iscell(csirsInfo.CDMType)
                        cdmType = {csirsInfo.CDMType};
                    end
                    mapping = containers.Map({'noCDM','FD-CDM2','CDM4','CDM8'},{[1 1],[2 1],[2 2],[2 4]});
                    cdmLen = mapping(cdmType{1});
                    % Estimated channel and noise variance
                    [Hest,nVar] = nrChannelEstimate(rxGrid, csirsRefInd, csirsSym, 'CDMLengths', cdmLen);

                    rank = obj.RankIndicator;
                    %YXC begin
%                     [cqi, pmiSet, ~, ~] = hCQISelect(carrier, csirsInfo, obj.CSIReportConfig, rank, Hest, nVar, obj.SINRTable);
                    [cqi, pmiSet, cqiInfo, ~] = hCQISelect(carrier, csirsInfo, obj.CSIReportConfig, rank, Hest, nVar, obj.SINRTable);
                    %Store cqiInfo into YusUtilityObj
                    obj.StoreCQIInfo(cqiInfo);
                    %YXC end
                    
                    % CQI value reported for each slot is stored in a new column
                    % In subband case, a column of CQI values is reported, where each element corresponds to each subband
                    % Convert CQI of sub-bands to per-RB CQI
                    cqiRBs = zeros(obj.CarrierInformation.NRBsDL, 1);
                    if strcmp(obj.CSIReportConfig.CQIMode, 'Subband')
                        subbandSize = obj.CSIReportConfig.SubbandSize;
                        cqiOffsetLevel = [0 1 2 -1]; %  Corresponding to offset values (0, 1, 2 and 3) as per TS 38.214 Table 5.2.2.1-1
                        % Fill same CQI for all the RBs in the sub-band
                        for i = 1:obj.CarrierInformation.NRBsDL/subbandSize
                            cqiRBs((i-1)*subbandSize+1 : i*subbandSize) = cqi(1) + cqiOffsetLevel(cqi(i+1)+1);
                        end
                        if mod(obj.CarrierInformation.NRBsDL, subbandSize)
                            cqiRBs((length(cqi)-2)*subbandSize+1 : end) = cqi(1) + cqiOffsetLevel(cqi(end)+1);
                        end
                    else
                        cqiRBs(:) = cqi(1); % Wideband CQI
                    end
                    cqiRBs(cqiRBs<=1) = 1; % Ensuring minimum CQI as 1
                    
                    % Report the CQI to MAC
                    obj.CSIRSIndicationFcn(rank, pmiSet, cqiRBs);
                end
            end
        end
       
        function waveformOut = applyTxPowerLevelAndGain(obj, waverformIn, gain)
            %applyTxPowerLevel Apply Tx power level to IQ samples
            
            % Apply Tx power to IQ samples
            scale = 10.^((-30 + obj.TxPower + gain)/20);
            waveformOut = waverformIn * scale;
        end
        
        function [waveformOut, pathloss] = applyPathLoss(obj, waveformIn, txInfo)
            %MXC_1
            %{
            %applyPathloss Apply free space path loss to the received waveform
            
            % Calculate the distance between source and destination nodes
            distance = getNodeDistance(obj.Node, txInfo.Position);
            % Wavelength
            lambda = physconst('LightSpeed')/txInfo.CarrierFreq;
            % Calculate the path loss
            pathloss = fspl(distance, lambda);
            %}
            
            gNBPos = txInfo.Position;
            UEPos = obj.Node.NodePosition;
            
%             LOS = obj.ChannelModel.HasLOSCluster;
            
            %Calculate new gNB position in wrap-around mode
            [~,gNBPos] = obj.Node.DistanceCalculatorFcn(gNBPos,UEPos);
            
            % Calculate pathloss
            YUF = YusUtilityFunctions;
            % Construct the structure of the parameters for input for
            % calculating pathloss
            paramForPL.Scenario = obj.YusUtilityParameter.Scenario;
            paramForPL.DLCarrierFreq = obj.ChannelModel.CarrierFrequency;
            LinkDir = 0; % 0 for DL, 1 for UL
            [pathloss, sigma_SF] = calculatePathloss(YUF, paramForPL, gNBPos, obj.YusUtilityParameter.UEStat, LinkDir);
            pathloss = normrnd(pathloss,sigma_SF); % Add shadow fading to pathloss
            %{
            % The procedure for calculating pathloss is moved to
            % YusUtilityFunctions.m
            %get distance and heights
            d_2D=norm(gNBPos(1:2)-UEPos(1:2));
            d_3D=norm(gNBPos(1:3)-UEPos(1:3));
            h_BS = gNBPos(3);
            h_UT=UEPos(3);
    
            %constants
            f_c = obj.ChannelModel.CarrierFrequency / 1e9; % frequency is stored in Hz, need in GHz for pathloss calculations
            
            %select the function to calculate LOS probability
            switch obj.YusUtilityParameter.Scenario
                case 'UMi'
                    %see 3GPP TR 38.901 Table 7.4.1-1
    
                    %generate coefficients
                    c = 300000000;
                    h_e = 1;
                    h_bs_prime = h_BS - h_e;
                    h_ut_prime = h_UT - h_e;
                    d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;
    
                    %calculate LOS pathloss
                    if (d_2D >= 10) && (d_2D <= d_bp_prime)
                        PL_umi_los = 32.4 + 21*log10(d_3D) + 20*log10(f_c);
                    elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
                        PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
                    else
                        %error('d_2D out of bound.'); 
                        PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
                    end 
    
                    if LOS
                        pathloss = PL_umi_los;
                    else
                        %calculate NLOS pathloss
                        PL_umi_nlos_prime = 35.3*log10(d_3D)+22.4+21.3*log10(f_c)-0.3*(h_UT-1.5);
                        pathloss = max(PL_umi_los, PL_umi_nlos_prime);
                    end
                    
                case 'UMa'
                    
                    %see 3GPP TR 38.901 Table 7.4.1-1
    
                    %generate coefficients
                    c = 300000000;
                    h_e = 1;
                    h_bs_prime = h_BS - h_e;
                    h_ut_prime = h_UT - h_e;
                    d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;
    
                    %calculate LOS pathloss
                    if (d_2D >= 10) && (d_2D <= d_bp_prime)
                        PL_uma_los = 28.0 + 22*log10(d_3D) + 20*log10(f_c);
                    elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
                        PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
                    else
                        %error('d_2D out of bound.');
                        PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
                    end 
    
                    if LOS
                        pathloss = PL_uma_los;
                    else
                        %calculate NLOS pathloss
                        PL_uma_nlos_prime = 13.54 + 39.08*log10(d_3D)+20*log10(f_c)-0.6*(h_UT-1.5);
                        pathloss = max(PL_uma_los, PL_uma_nlos_prime);
                    end
                    
                case 'RMa'
                    %generate coefficients
                    c = 300000000;
                    h = 5;
                    w = 20;
                    d_bp = (2 * pi * h_BS * h_UT * f_c * 1000000000) / c;
    
                    %calculate LOS pathloss
                    if (d_2D >= 10) && (d_2D <= d_bp)
                        PL_rma_los = 20*log10(40*pi*d_3D*f_c/3)+min(0.03*(h^1.72),10)*log10(d_3D)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_3D;
                    elseif (d_2D >= d_bp) && (d_2D <=10000)
                        PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
                    else
                        %error('d_2D out of bound.');
                        PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
                    end 
    
                    if LOS
                        pathloss = PL_rma_los;
                    else
                        %calculate NLOS pathloss
                        PL_rma_nlos_prime = 161.04-7.1*log10(w)+7.5*log10(h)-(24.34-((h/h_BS)^2))*log10(h_BS)+(43.42-3.1*log10(h_BS))*(log10(d_3D)-3)+20*log10(f_c)-(3.2*((log10(11.75*h_UT))^2)-4.97);
                        pathloss = max(PL_rma_los, PL_rma_nlos_prime);
                    end
                    
                otherwise
                    error('Scenarios other than RMa, UMi, UMa are yet to be implemented.');
            end
            %}
            
            %MXC_1
            % Apply path loss on IQ samples
            scale = 10.^(-pathloss/20);
            waveformOut = waveformIn * scale;
        end
        
        function waveformOut = applyRxGain(obj, waveformIn)
            %applyRxGain Apply receiver antenna gain
            
            scale = 10.^(obj.RxGain/20);
            waveformOut = waveformIn.* scale;
        end
        
        function waveformOut = applyThermalNoise(obj, waveformIn, pktInfo, selfInfo)
            %applyThermalNoise Apply thermal noise
            
            % Thermal noise(in Watts) = BoltzmannConstant * Temperature (in Kelvin) * bandwidth of the channel.
            Nt = physconst('Boltzmann') * selfInfo.Temperature * selfInfo.Bandwidth;
            thermalNoise = (10^(obj.NoiseFigure/10)) * Nt; % In watts
            totalnoise = thermalNoise;
            % Calculate SNR
            SNR = pktInfo.TxPower - ((10*log10(totalnoise)) + 30);
            % Add noise
            waveformOut = awgn(waveformIn,SNR,pktInfo.TxPower-30,'db');
        end
        
        function logPackets(obj, info, macPDU, linkDir)
            %logPackets Capture the MAC packets to a PCAP file
            %
            % logPackets(OBJ, INFO, MACPDU, LINKDIR)
            %
            % INFO - Contains the PUSCH/PDSCH information
            %
            % MACPDU - MAC PDU
            %
            % LINKDIR - 1 represents UL and 0 represents DL direction
            
            timestamp = round(obj.getCurrentTime()); % Timestamp in microseconds
            obj.PacketMetaData.HARQID = info.HARQID;
            obj.PacketMetaData.SlotNumber = info.NSlot;
            
            if linkDir % Uplink
                obj.PacketMetaData.SystemFrameNumber = mod(obj.AFN, 1024);
                obj.PacketMetaData.LinkDir = obj.PacketLogger.Uplink;
            else % Downlink
                % Get frame number of previous slot i.e the Tx slot. Reception ended at the
                % end of previous slot
                if obj.CurrSlot > 0
                    prevSlotAFN = obj.AFN; % Previous slot was in the current frame
                else
                    % Previous slot was in the previous frame
                    prevSlotAFN = obj.AFN - 1;
                end
                obj.PacketMetaData.SystemFrameNumber = mod(prevSlotAFN, 1024);
                obj.PacketMetaData.LinkDir = obj.PacketLogger.Downlink;
            end
            write(obj.PacketLogger, macPDU, timestamp, 'PacketInfo', obj.PacketMetaData);
        end
    end

    methods (Hidden = true)
        function dlTTIRequest(obj, pduType, dlControlPDU)
            dlControlRequest(obj, pduType, dlControlPDU);
        end
    end
end