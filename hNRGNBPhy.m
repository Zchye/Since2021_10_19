classdef hNRGNBPhy < hNRPhyInterface
    %hNRGNBPhy 5G NR Phy Tx and Rx processing chains at gNB
    %   The class implements the Phy Tx and Rx processing chains of 5G NR
    %   at gNB. It also implements the interfaces for information exchange
    %   between Phy and higher layers. It supports transmission of physical
    %   downlink shared channel (PDSCH) along with its demodulation
    %   reference signals (DM-RS), and channel state information reference
    %   signals (CSI-RS). It only supports reception of physical uplink
    %   shared channel (PUSCH) along with its DM-RS. gNB is assumed to
    %   serve a single cell. A single bandwidth part is assumed to cover
    %   the entire carrier bandwidth. Note that setCarrierInformation and 
    %   setCellConfig need to be called on the created class object before
    %   using it.
    %
    %   hNRGNBPhy methods:
    %       hNRGNBPhy               - Construct a gNB Phy object
    %       run                     - Run the gNB Phy layer operations
    %       setCarrierInformation   - Set the carrier configuration
    %       enablePacketLogging     - Enable packet logging
    %       registerMACInterfaceFcn - Register MAC interface functions at 
    %                                 Phy, for sending information to MAC
    %       registerInBandTxFcn     - Register callback for transmission 
    %                                 on PDSCH
    %       txDataRequest           - Tx request from MAC to Phy for 
    %                                 starting PDSCH transmission
    %       dlControlRequest        - Downlink control (non-data) transmission 
    %                                 request from MAC to Phy
    %       rxDataRequest           - Rx request from MAC to Phy for 
    %                                 starting PUSCH reception
    %       phyTx                   - Physical layer processing and 
    %                                 transmission
    %       storeReception          - Receive the incoming waveform and add
    %                                 it to the reception buffer
    %       phyRx                   - Physical layer reception and sending 
    %                                 of decoded information to MAC layer
    %       getULBLER               - Get block error statistics of the slot
    %                                 for each UE
    % 
    %   Example:
    %   % Generate a hNRGNBPhy object. Configure the carrier and cell 
    %   % properties using setCarrierInformation and setCellConfig methods,
    %   respectively
    %
    %   phyParam = struct();
    %   phyParam.NumUEs = 1;
    %   phyParam.SCS = 15;
    %   phyParam.NumRBs = 52;
    %   phyParam.GNBTxPower = 29;
    %   phyParam.GNBRxGain = 11;
    %   phyParam.GNBRxBufferSize = 1;
    %   phyParam.ULCarrierFreq = 2.1e9;
    %   phyParam.ChannelModelType = 'CDL';
    %   phy = hNRGNBPhy(phyParam);
    %
    %   carrierParam = struct();
    %   carrierParam.SubcarrierSpacing = phyParam.SCS;
    %   carrierParam.NRBsUL = phyParam.NumRBs;
    %   carrierParam.NRBsDL = phyParam.NumRBs;
    %   carrierParam.DLFreq = phyParam.ULCarrierFreq;
    %   setCarrierInformation(phy, carrierParam);
    %
    %   cellParam.NCellID = 1;
    %   cellParam.DuplexMode = 0;
    %   setCellConfig(phy, cellParam);
    %
    %   See also hNRPhyInterface
    
    %   Copyright 2020-2021 The MathWorks, Inc.
    
    properties (Access = private)
        %UEs RNTIs in the cell
        UEs
        
        %DLSCHEncoders Downlink shared channel (DL-SCH) encoder system objects for the UEs
        % Vector of length equal to the number of UEs in the cell. Each
        % element is an object of type nrDLSCH
        DLSCHEncoders
        
        %ULSCHDecoders Uplink shared channel (UL-SCH) decoder system objects for the UEs
        % Vector of length equal to the number of UEs in the cell. Each
        % element is an object of type nrULSCHDecoder
        ULSCHDecoders
        
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
        TxPower(1, 1) {mustBeFinite, mustBeNonnegative, mustBeNonNan} = 29;
        
        %RxGain Rx antenna gain in dBi
        RxGain(1, 1) {mustBeFinite, mustBeNonnegative, mustBeNonNan} = 11;
        
        %CSIRSPDU CSI-RS information PDU sent by MAC for the current slot
        % It is an object of type nrCSIRSConfig containing the
        % configuration of CSI-RS to be sent in current slot. If empty,
        % then CSI-RS is not scheduled for the current slot
        CSIRSPDU = {}
        
        %PDSCHPDU PDSCH information sent by MAC for the current slot
        % It is an array of objects of type hNRPDSCHInfo. An object at
        % index 'i' contains the information required by Phy to transmit a
        % MAC PDU stored at index 'i' of object property 'MacPDU'
        PDSCHPDU = {}
        
        %MacPDU PDUs sent by MAC which are scheduled to be sent in the currrent slot
        % It is an array of downlink MAC PDUs to be sent in the current
        % slot. Each object in the array corresponds to one object in
        % object property PDSCHPDU
        MacPDU = {}
        
        %ULBlkErr Uplink block error information
        % It is an array of size N-by-2 where N is the number of UEs,
        % columns 1 and 2 contains the number of erroneously received
        % packets and total received packets, respectively
        ULBlkErr
        
        %PrevCumulativeULBlkErr Cumulative uplink block error information returned in the last query
        % It is an array of size N-by-2 where N is the number of UEs,
        % columns 1 and 2 contains the number of erroneously received
        % packets and total received packets, respectively
        PrevCumulativeULBlkErr
        
        %ChannelModel Information about the propagation channel model
        % It is a cell array of length equal to the number of UEs. The
        % array contains objects of type nrCDLChannel, if the channel model
        % type is specified as 'CDL', otherwise empty. An object at index
        % 'i' models the channel between the gNB and UE with RNTI 'i'
        ChannelModel
        
        %MaxChannelDelay Maximum delay introduced by multipath components and implementation delays
        % It is an array of length equal to the number of UEs. Each element
        % at index 'i' corresponds to maximum channel delay between the gNB
        % and UE with RNTI 'i'
        MaxChannelDelay
        
        %TimingOffset Receiver timing offset
        % Receiver timing offset used for practical synchronization. It is
        % an array of length equal to the number of UEs. Each element at
        % index 'i' corresponds to the timing offset experienced during
        % reception of waveform from UE with RNTI 'i'
        TimingOffset
        
        %NoiseFigure Noise figure at the receiver
        NoiseFigure = 1;
        
        %RxBuffer Reception buffer object to store received waveforms
        RxBuffer
        
        %SendWaveformFcn Function handle to transmit the waveform
        SendWaveformFcn
        
        %Temperature Temperature at node in Kelvin
        % It is used for thermal noise calculation
        Temperature = 300
        
        %PacketLogger Contains handle of the packet capture (PCAP) object
        PacketLogger
        
        %PacketMetaData Contains the information required for logging MAC packets into PCAP file
        PacketMetaData
        
        %MXC_1
        siteIdx
        
        YusUtilityParameter
        %MXC_1
    end
    
    methods
        function obj = hNRGNBPhy(param, siteIdx)
            %MXC_1 extra constructor input
            %hNRGNBPhy Construct a gNB Phy object
            %   OBJ = hNRGNBPHY(numUEs) constructs a gNB Phy object. It
            %   also creates the context of DL-SCH encoders system objects
            %   and UL-SCH decoders system objects for all the UEs.
            %
            %   PARAM is structure with the fields:
            %
            %       NumUEs           - Number of UEs in the cell
            %       SCS              - Subcarrier spacing
            %       NumRBs           - Number of resource blocks
            %       GNBTxPower       - Tx Power in dBm
            %       GNBRxGain        - Receiver antenna gain at gNB in dBi
            %       GNBRxBufferSize  - Maximum number of waveforms to be
            %                          stored
            %       ULCarrierFreq    - Uplink carrier frequency in Hz
            %       ChannelModelType - Propagation channel model type
            %       GNBTxAnts        - Number of GNB Tx antennas
            %       GNBRxAnts        - Number of GNB Rx antennas
            %       UETxAnts         - Number of Tx antennas on UEs. Vector of length 'N' where N is number of UEs.
            %                          Value at index 'i' contains Tx antennas at UE with RNTI 'i'
            %       UERxAnts         - Number of Rx antennas on UEs. Vector of length 'N' where N is number of UEs.
            %                          Value at index 'i' contains Rx antennas at UE with RNTI 'i'
            
            % Validate the number of UEs
            validateattributes(param.NumUEs, {'numeric'}, {'nonempty', 'integer', 'scalar', '>', 0, '<=', 65519}, 'param.NumUEs', 'NumUEs')
            
            obj.UEs = 1:param.NumUEs;
            
            %MXC_1
            obj.siteIdx = siteIdx;
            
            obj.YusUtilityParameter.Scenario = param.Scenario;
            %MXC_1
            
            % Create DL-SCH encoder system objects for the UEs
            obj.DLSCHEncoders = cell(param.NumUEs, 1);
            for i=1:param.NumUEs
                obj.DLSCHEncoders{i} = nrDLSCH;
                obj.DLSCHEncoders{i}.MultipleHARQProcesses = true;
            end
            
            % Create UL-SCH decoder system objects for the UEs
            obj.ULSCHDecoders = cell(param.NumUEs, 1);
            for i=1:param.NumUEs
                obj.ULSCHDecoders{i} = nrULSCHDecoder;
                obj.ULSCHDecoders{i}.MultipleHARQProcesses = true;
                obj.ULSCHDecoders{i}.LDPCDecodingAlgorithm = 'Normalized min-sum';
                obj.ULSCHDecoders{i}.MaximumLDPCIterationCount = 6;
            end
            
            % Set the number of erroneous packets and the total number of
            % packets received from each UE to zero
            obj.ULBlkErr = zeros(param.NumUEs, 2);
            obj.PrevCumulativeULBlkErr = zeros(param.NumUEs, 2);
            
            % Set Tx power in dBm
            if isfield(param, 'GNBTxPower')
                obj.TxPower = param.GNBTxPower;
            end
            % Set Rx antenna gain in dBi
            if isfield(param, 'GNBRxGain')
                obj.RxGain = param.GNBRxGain;
            end
            
            % Initialize the ChannelModel and MaxChannelDelay properties
            obj.ChannelModel = cell(1, param.NumUEs);
            obj.MaxChannelDelay = zeros(1, param.NumUEs);
            
            % Initialize timing offsets to 0
            obj.TimingOffset = zeros(1, param.NumUEs);
            
            waveformInfo = nrOFDMInfo(param.NumRBs, param.SCS);

            %MXC_2
            
            if ~isfield(param, 'GNBTxAnts')
                param.GNBTxAnts = 1;
            % Validate the number of transmitter antennas on gNB
            elseif ~ismember(param.GNBTxAnts, [1,2,4,8,16,32,64,128,256,512,1024])
                error('nr5g:hNRGNBPhy:InvalidAntennaSize',...
                    'Number of gNB Tx antennas (%d) must be a member of [1,2,4,8,16,32,64,128,256,512,1024].', param.GNBTxAnts);
            end
            if ~isfield(param, 'GNBRxAnts')
                param.GNBRxAnts = 1;
            % Validate the number of receiver antennas on gNB
            elseif ~ismember(param.GNBRxAnts, [1,2,4,8,16,32,64,128,256,512,1024])
                error('nr5g:hNRGNBPhy:InvalidAntennaSize',...
                    'Number of gNB Rx antennas (%d) must be a member of [1,2,4,8,16,32,64,128,256,512,1024].', param.GNBRxAnts);
            end
            if ~isfield(param, 'UETxAnts')
                param.UETxAnts = ones(1, param.NumUEs);
            % Validate the number of transmitter antennas on UEs
            elseif any(~ismember(param.UETxAnts, [1,2,4,8,16]))
                error('nr5g:hNRGNBPhy:InvalidAntennaSize',...
                    'Number of UE Tx antennas must be a member of [1,2,4,8,16].');
            end
            if ~isfield(param, 'UERxAnts')
                param.UERxAnts = ones(1, param.NumUEs);
            % Validate the number of receiver antennas on UEs
            elseif any(~ismember(param.UERxAnts, [1,2,4,8,16]))                
                error('nr5g:hNRGNBPhy:InvalidAntennaSize',...
                    'Number of UE Rx antennas must be a member of [1,2,4,8,16].');
            end
            
            %{
            ueTxAntPanel = ones(param.NumUEs, 5);
            for i = 1:param.NumUEs
                % Downlink direction i.e. gNB's Tx antenna panel
                % configuration and UEs' Rx panel configuration
                [obj.TxAntPanel] = ...
                    hArrayGeometry(param.GNBTxAnts, param.UERxAnts(i), 'downlink');

                % Uplink direction i.e. UEs' Tx antenna panel
                % configuration and gNB's Rx panel configuration
                [ueTxAntPanel(i, :), obj.RxAntPanel] = ...
                    hArrayGeometry( param.UETxAnts(i), param.GNBRxAnts, 'uplink');
            end
            %}
            
            %MXC_2
            if isfield(param, 'ChannelModelType')
                if strcmpi(param.ChannelModelType, 'CDL')
                    for ueIdx = 1:param.NumUEs
                        %MXC_1
                        %{
                        channel = nrCDLChannel; % CDL channel object
                        channel.DelayProfile = 'CDL-C';
                        channel.DelaySpread = 300e-9;
                        % Set the carrier frequency to uplink
                        channel.CarrierFrequency = param.ULCarrierFreq;
                        %}
                        channel = nrCDLChannel; % CDL channel object
                        channel.DelayProfile = 'Custom';
                    
                        LinkDir = 1; % UL
                        CHParam = ConfigChannel(param, LinkDir, obj.siteIdx, ueIdx);
                        
                        channel.CarrierFrequency = param.ULCarrierFreq;
                        channel.HasLOSCluster = logical(CHParam.LOS);
                        channel.AngleSpreads = [CHParam.LSPs.ASD, CHParam.LSPs.ASA, CHParam.LSPs.ZSD, CHParam.LSPs.ZSA];
                        if CHParam.LOS
                            channel.KFactorFirstCluster = CHParam.LSPs.K;
                        end
                        channel.PathDelays = CHParam.tau_n';
                        channel.AveragePathGains = CHParam.P_n';
                        channel.AnglesAoA = CHParam.ADAngles.phi_n_AOA';
                        channel.AnglesAoD = CHParam.ADAngles.phi_n_AOD';
                        channel.AnglesZoA = CHParam.ADAngles.theta_n_ZOA';
                        channel.AnglesZoD = CHParam.ADAngles.theta_n_ZOD';
                        channel.XPR = CHParam.XPR;
                        %MXC_1 quick fix
                        %causes problems with 100M element array, change to zero for now
                        channel.NumStrongestClusters = 2;
%                         channel.NumStrongestClusters = 0;
                        %Mxc_1
                        %not relevant anymore due to quick fix
                        channel.ClusterDelaySpread = CHParam.c_DS;
                        
                        %MXC_1
                        
                        % Size of antenna array [M N P Mg Ng], where:
                        %    - M and N are the number of rows and columns in
                        %      the antenna array, respectively.
                        %    - P is the number of polarizations (1 or 2).
                        %    - Mg and Ng are the number of row and column
                        %      array panels, respectively.
                        
                        %MXC_2
                        obj.TxAntPanel = param.GNBTxAntPanelSize;
                        obj.RxAntPanel = param.GNBRxAntPanelSize;
                        
                        %UE Antenna Config
                        %channel.TransmitAntennaArray.Size = ueTxAntPanel(ueIdx, :);
                        channel.TransmitAntennaArray.Size = param.UETxAntPanelSize;
                        channel.TransmitAntennaArray.ElementSpacing = param.UETxAntElementSpacing;
                        channel.TransmitAntennaArray.PolarizationAngles = param.UETxAntPolarizationAngles;
                        %channel.TransmitAntennaArray.Orientation = [0; 0; 0];
                        channel.TransmitAntennaArray.Element = param.UEAntElement;
                        channel.TransmitAntennaArray.PolarizationModel = param.UEAntPolarizationModel;
                        channel.TransmitArrayOrientation = [0; 0; 0];
                        
                        %gNB Antenna Config
                        channel.ReceiveAntennaArray.Size = obj.RxAntPanel;
                        channel.ReceiveAntennaArray.ElementSpacing = param.GNBRxAntElementSpacing;
                        channel.ReceiveAntennaArray.PolarizationAngles = param.GNBRxAntPolarizationAngles;
                        %channel.ReceiveAntennaArray.Orientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
                        channel.ReceiveAntennaArray.Element = param.GNBAntElement;
                        channel.ReceiveAntennaArray.PolarizationModel = param.GNBAntPolarizationModel;
                        channel.ReceiveArrayOrientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
                        %MXC_2
                        
                        % UE Mobility
                        c = physconst('lightspeed'); % speed of light in m/s
                        fd = param.UEStates{obj.siteIdx,ueIdx}.Speed/c*param.ULCarrierFreq; % UE max Doppler frequency in Hz
                        obj.ChannelModel.MaximumDopplerShift = fd;
                        obj.ChannelModel.UTDirectionOfTravel = [param.UEStates{obj.siteIdx,ueIdx}.DirectionOfTravel; 90];
                        
                        channel.SampleRate = waveformInfo.SampleRate;
                        chInfo = info(channel);
                        obj.ChannelModel{ueIdx} = channel;
                        % Update the maximum delay caused due to CDL channel model
                        obj.MaxChannelDelay(ueIdx) = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
                    end
                end
            end
            
            % Set receiver noise figure
            if isfield(param, 'NoiseFigure')
                obj.NoiseFigure = param.NoiseFigure;
            end
            
            % Create reception buffer object
            if isfield(param, 'GNBRxBufferSize')
                obj.RxBuffer = hNRPhyRxBuffer('BufferSize', param.GNBRxBufferSize, 'NRxAnts', param.GNBRxAnts);
            else
                obj.RxBuffer = hNRPhyRxBuffer('NRxAnts', param.GNBRxAnts);
            end
        end
        
        function run(obj)
            %run Run the gNB Phy layer operations
            
            % Phy processing and transmission of PDSCH (along with its
            % DM-RS) and CSI-RS. It is assumed that MAC has already loaded
            % the Phy Tx context for anything scheduled to be transmitted
            % at the current symbol
            phyTx(obj);
            
            % Phy reception of PUSCH and sending decoded information to
            % MAC. Receive the PUSCHs which ended in the last symbol.
            % Reception as well as processing is done in the symbol after
            % the last symbol in PUSCH duration (till then the packets are
            % queued in Rx buffer). Phy calculates the last symbol of PUSCH
            % duration based on 'rxDataRequest' call from MAC (which comes
            % at the first symbol of PUSCH Rx time) and the PUSCH duration
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
        end
        
        function registerMACInterfaceFcn(obj, sendMACPDUFcn)
            %registerMACInterfaceFcn Register MAC interface functions at Phy, for sending information to MAC
            %   registerMACInterfaceFcn(OBJ, SENDMACPDUFCN) registers the
            %   callback function to send decoded MAC PDUs to MAC.
            %
            %   SENDMACPDUFCN Function handle provided by MAC to Phy for
            %   sending PDUs to MAC.
            
            obj.RxIndicationFcn = sendMACPDUFcn;
        end
        
        function registerInBandTxFcn(obj, sendWaveformFcn)
            %registerInBandTxFcn Register callback for transmission on PDSCH
            %
            %   SENDWAVEFORMFCN Function handle provided by packet
            %   distribution object for packet transmission
            
            obj.SendWaveformFcn = sendWaveformFcn;
        end
        
        function txDataRequest(obj, PDSCHInfo, macPDU)
            %txDataRequest Tx request from MAC to Phy for starting PDSCH transmission
            %  txDataRequest(OBJ, PDSCHINFO, MACPDU) sets the Tx context to
            %  indicate PDSCH transmission in the current symbol
            %
            %  PDSCHInfo is an object of type hNRPDSCHInfo, sent by MAC. It
            %  contains the information required by the Phy for the
            %  transmission.
            %
            %  MACPDU is the downlink MAC PDU sent by MAC for transmission.
            
            % Update the Tx context. There can be multiple simultaneous
            % PDSCH transmissions for different UEs
            obj.MacPDU{end+1} = macPDU;
            obj.PDSCHPDU{end+1} = PDSCHInfo;
        end
        
        function dlControlRequest(obj, pduType, dlControlPDU)
            %dlControlRequest Downlink control (non-data) transmission request from MAC to Phy
            %   dlControlRequest(OBJ, PDUTYPES, DLCONTROLPDUS) is a request from
            %   MAC for downlink transmission. MAC sends it at the start of
            %   a DL slot for all the scheduled non-data DL transmission in
            %   the slot (Data i.e. PDSCH is sent by MAC using
            %   txDataRequest interface of this class).
            %
            %   PDUTYPE is an array of packet types. Currently, only packet
            %   type 0 (CSI-RS) is supported.
            %
            %   DLCONTROLPDU is an array of DL control information PDUs. Each PDU
            %   is stored at the index corresponding to its type in
            %   PDUTYPE. Currently supported CSI-RS information PDU is an object of
            %   type nrCSIRSConfig.
            
            % Update the Tx context
            for i = 1:length(pduType)
                switch(pduType(i))
                    case obj.CSIRSPDUType
                        obj.CSIRSPDU = dlControlPDU{i};
                end
            end
        end
        
        function rxDataRequest(obj, puschInfo)
            %rxDataRequest Rx request from MAC to Phy for starting PUSCH reception
            %   rxDataRequest(OBJ, PUSCHINFO) is a request to start PUSCH
            %   reception. It starts a timer for PUSCH end time (which on
            %   triggering receives the complete PUSCH). The Phy expects
            %   the MAC to send this request at the start of reception
            %   time.
            %
            %   PUSCHInfo is an object of type hNRPUSCHInfo. It contains
            %   the information required by the Phy for the reception.
            
            symbolNumFrame = obj.CurrSlot*14 + obj.CurrSymbol; % Current symbol number w.r.t start of 10 ms frame
            
            % PUSCH to be read in the symbol after the last symbol in
            % PUSCH reception
            numPUSCHSym =  puschInfo.PUSCHConfig.SymbolAllocation(2);
            puschRxSymbolFrame = mod(symbolNumFrame + numPUSCHSym, obj.CarrierInformation.SymbolsPerFrame);
            
            % Add the PUSCH Rx information at the index corresponding to
            % the symbol just after PUSCH end time
            obj.DataRxContext{puschRxSymbolFrame+1}{end+1} = puschInfo;
        end
        
        function phyTx(obj)
            %phyTx Physical layer processing and transmission
            
            if isempty(obj.PDSCHPDU) && isempty(obj.CSIRSPDU)
                return; % No transmission (PDSCH or CSI-RS) is scheduled to start at the current symbol
            end
            
            % Calculate Tx waveform length in symbols
            numTxSymbols = 0;
            if ~isempty(obj.CSIRSPDU)
                % Transmit any CSI-RS scheduled in the slot at the slot start
                numTxSymbols = obj.CSIRSPDU.SymbolLocations(end)+1;
            end
            if ~isempty(obj.PDSCHPDU) % If PDSCH(s) start at the current symbol
                % Among all the PDSCHs scheduled to be transmitted now, get
                % the duration in symbols of the PDSCH which spans maximum
                % number of symbols
                for i = 1:length(obj.PDSCHPDU)
                    pdschLen = obj.PDSCHPDU{i}.PDSCHConfig.SymbolAllocation(2);
                    if(pdschLen > numTxSymbols)
                        numTxSymbols = pdschLen;
                    end
                end
            end
           
            % Initialize Tx grid
            numTxAnts = prod(obj.TxAntPanel);
            txGrid = zeros(obj.CarrierInformation.NRBsDL*12, obj.WaveformInfoDL.SymbolsPerSlot, numTxAnts);
           
            % Set carrier configuration object
            carrier = nrCarrierConfig;
            carrier.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
            carrier.NSizeGrid = obj.CarrierInformation.NRBsDL;
            carrier.NSlot = obj.CurrSlot;
            carrier.NFrame = obj.AFN;
            
            % Fill CSI-RS in the grid
            if ~isempty(obj.CSIRSPDU)
                csirsInd = nrCSIRSIndices(carrier, obj.CSIRSPDU);
                csirsSym = nrCSIRS(carrier, obj.CSIRSPDU);
                % Placing the CSI-RS in the Tx grid
                txGrid(csirsInd) = csirsSym;
                obj.CSIRSPDU = {};
            end
            
            % Fill PDSCH symbols in the grid
            if ~isempty(obj.PDSCHPDU)
                txGrid = populatePDSCH(obj, obj.PDSCHPDU, obj.MacPDU, txGrid);
            end
            
            % OFDM modulation
            txWaveform = nrOFDMModulate(carrier, txGrid);
            
            % Trim txWaveform to span only the transmission symbols
            slotNumSubFrame = mod(obj.CurrSlot, obj.WaveformInfoDL.SlotsPerSubframe);
            startSymSubframe = slotNumSubFrame*obj.WaveformInfoDL.SymbolsPerSlot + 1; % Start symbol of current slot in the subframe
            lastSymSubframe = startSymSubframe + obj.WaveformInfoDL.SymbolsPerSlot -1; % Last symbol of current slot in the subframe
            symbolLengths = obj.WaveformInfoDL.SymbolLengths(startSymSubframe : lastSymSubframe); % Length of symbols of current slot
            startSample = sum(symbolLengths(1:obj.CurrSymbol)) + 1;
            endSample = sum(symbolLengths(1:obj.CurrSymbol+numTxSymbols));
            txWaveform = txWaveform(startSample:endSample, :);
            
            % Apply Tx power and gain
            gain = 0;
            txWaveform = applyTxPowerLevelAndGain(obj, txWaveform, gain);
            
            % Construct packet information
            packetInfo.Waveform = txWaveform;
            packetInfo.Position = obj.Node.NodePosition;
            packetInfo.CarrierFreq = obj.CarrierInformation.DLFreq;
            packetInfo.TxPower = obj.TxPower;
            packetInfo.NTxAnts = numTxAnts;
            packetInfo.SampleRate = obj.WaveformInfoDL.SampleRate;
            
            % Waveform transmission by sending it to packet
            % distribution entity
            obj.SendWaveformFcn(packetInfo);
            
            % Clear the Tx contexts
            obj.PDSCHPDU = {};
            obj.MacPDU = {};
        end

        function storeReception(obj, waveformInfo)
            %storeReception Receive the incoming waveform and add it to the reception
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
            puschInfoList = obj.DataRxContext{symbolNumFrame + 1};
            currentTime = getCurrentTime(obj);
            
            for i = 1:length(puschInfoList) % For all PUSCH receptions which ended in the last symbol
                puschInfo = puschInfoList{i};
                startSymPUSCH = puschInfo.PUSCHConfig.SymbolAllocation(1);
                numSymPUSCH = puschInfo.PUSCHConfig.SymbolAllocation(2);
                % Calculate the symbol start index w.r.t start of 1 ms sub frame
                slotNumSubFrame = mod(puschInfo.NSlot, obj.WaveformInfoUL.SlotsPerSubframe);
                % Calculate PUSCH duration
                puschSymbolSet = startSymPUSCH : startSymPUSCH+numSymPUSCH-1;
                symbolSetSubFrame = (slotNumSubFrame * 14) + puschSymbolSet + 1;
                duration = 1e6 * (1/obj.WaveformInfoUL.SampleRate) * ...
                    sum(obj.WaveformInfoUL.SymbolLengths(symbolSetSubFrame)); % In microseconds
                
                % Convert channel delay into microseconds
                maxChannelDelay = 1e6 * (1/obj.WaveformInfoUL.SampleRate) * obj.MaxChannelDelay(puschInfo.PUSCHConfig.RNTI);
                
                % Get the received waveform
                duration = duration + maxChannelDelay;
                rxWaveform = getReceivedWaveform(obj.RxBuffer, currentTime + maxChannelDelay - duration, duration, obj.WaveformInfoUL.SampleRate);
                
                % Process the waveform and send the decoded information to MAC
                phyRxProcessing(obj, rxWaveform, puschInfo);
            end
           
            % Clear the Rx context
            obj.DataRxContext{symbolNumFrame + 1} = {};
        end
        
        function ulBLER = getULBLER(obj)
            %getULBLER Get block error statistics of the slot for each UE
            
            ulBLER = obj.ULBlkErr - obj.PrevCumulativeULBlkErr;
            obj.PrevCumulativeULBlkErr = obj.ULBlkErr;
        end

        function timestamp = getCurrentTime(obj)
            %getCurrentTime Return the current timestamp of node in microseconds

            % Calculate number of samples from the start of the current
            % frame to the current symbol
            numSubFrames = floor(obj.CurrSlot / obj.WaveformInfoUL.SlotsPerSubframe);
            numSlotSubFrame = mod(obj.CurrSlot, obj.WaveformInfoUL.SlotsPerSubframe);
            symbolNumSubFrame = numSlotSubFrame*obj.WaveformInfoUL.SymbolsPerSlot + obj.CurrSymbol;
            numSamples = (numSubFrames * sum(obj.WaveformInfoUL.SymbolLengths))...
                + sum(obj.WaveformInfoUL.SymbolLengths(1:symbolNumSubFrame));

            % Timestamp in microseconds
            timestamp = (obj.AFN * 0.01) + (numSamples *  1 / obj.WaveformInfoUL.SampleRate);
            timestamp = (1e6 * timestamp);
        end
    end

    methods (Access = private)
        function setWaveformProperties(obj, carrierInformation)
            %setWaveformProperties Set the UL and DL waveform properties
            %   setWaveformProperties(OBJ, CARRIERINFORMATION) sets the UL
            %   and DL waveform properties as per the information in
            %   CARRIERINFORMATION. CARRIERINFORMATION is a structure
            %   including the following fields:
            %       SubcarrierSpacing  - Subcarrier spacing used
            %       NRBsDL             - Downlink bandwidth in terms of
            %                            number of resource blocks
            %       NRBsUL             - Uplink bandwidth in terms of
            %                            number of resource blocks
            
            % Set the UL waveform properties
            obj.WaveformInfoUL = nrOFDMInfo(carrierInformation.NRBsUL, carrierInformation.SubcarrierSpacing);
            
            % Set the DL waveform properties
            obj.WaveformInfoDL = nrOFDMInfo(carrierInformation.NRBsDL, carrierInformation.SubcarrierSpacing);
        end
        
        function updatedSlotGrid = populatePDSCH(obj, pdschPDU, macPDU, txSlotGrid)
            %populatePDSCH Populate PDSCH symbols in the Tx grid and return the updated grid
            
            for i=1:length(pdschPDU) % For each PDSCH scheduled for this slot
                pdschInfo = pdschPDU{i};
                % Set transport block in the encoder. In case of empty MAC
                % PDU sent from MAC (indicating retransmission), no need to set transport
                % block as it is already buffered in DL-SCH encoder object
                if ~isempty(macPDU{i})
                    % A non-empty MAC PDU is sent by MAC which indicates new
                    % transmission
                    macPDUBitmap = de2bi(macPDU{i}, 8);
                    macPDUBitmap = reshape(macPDUBitmap', [], 1); % Convert to column vector
                    setTransportBlock(obj.DLSCHEncoders{pdschInfo.PDSCHConfig.RNTI}, macPDUBitmap, 0, pdschInfo.HARQID);
                end
                
                if ~isempty(obj.PacketLogger) % Packet capture enabled
                    % Log downlink packets
                    if isempty(macPDU{i})
                        tbID = 0; % Transport block id
                        macPDUBitmap = getTransportBlock(obj.DLSCHEncoders{pdschInfo.PDSCHConfig.RNTI}, tbID, pdschInfo.HARQID);
                        macPDUBitmap = reshape(macPDUBitmap, 8, [])';
                        macPacket = bi2de(macPDUBitmap);
                        logPackets(obj, pdschInfo, macPacket, 0);
                    else
                        logPackets(obj, pdschInfo, macPDU{i}, 0);
                    end
                end
                W = pdschInfo.PrecodingMatrix;
                
                % Calculate PDSCH and DM-RS information
                carrierConfig = nrCarrierConfig;
                carrierConfig.NSizeGrid = obj.CarrierInformation.NRBsDL;
                carrierConfig.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
                carrierConfig.NSlot = pdschInfo.NSlot;
                carrierConfig.NCellID = obj.CellConfig.NCellID;
                [pdschIndices, pdschIndicesInfo] = nrPDSCHIndices(carrierConfig, pdschInfo.PDSCHConfig);
                dmrsSymbols = nrPDSCHDMRS(carrierConfig, pdschInfo.PDSCHConfig);
                dmrsIndices = nrPDSCHDMRSIndices(carrierConfig, pdschInfo.PDSCHConfig);
                
                % Encode the DL-SCH transport blocks
                obj.DLSCHEncoders{pdschInfo.PDSCHConfig.RNTI}.TargetCodeRate = pdschInfo.TargetCodeRate;
                codedTrBlock = step(obj.DLSCHEncoders{pdschInfo.PDSCHConfig.RNTI}, pdschInfo.PDSCHConfig.Modulation, ...
                    pdschInfo.PDSCHConfig.NumLayers, pdschIndicesInfo.G, pdschInfo.RV,pdschInfo.HARQID);

                % PDSCH modulation and precoding
                pdschSymbols = nrPDSCH(carrierConfig, pdschInfo.PDSCHConfig, codedTrBlock);
                [pdschAntSymbols, pdschAntIndices] = hPRGPrecode(size(txSlotGrid), carrierConfig.NStartGrid, pdschSymbols, pdschIndices, W);
                txSlotGrid(pdschAntIndices) = pdschAntSymbols;
                
                % PDSCH DM-RS precoding and mapping
                [dmrsAntSymbols, dmrsAntIndices] = hPRGPrecode(size(txSlotGrid), carrierConfig.NStartGrid, dmrsSymbols, dmrsIndices, W);
                txSlotGrid(dmrsAntIndices) = dmrsAntSymbols;
            end
            updatedSlotGrid = txSlotGrid;
        end
        
        function rxWaveform = applyChannelModel(obj, pktInfo)
            %applyChannelModel Return the waveform after applying channel model
            
            rxWaveform = pktInfo.Waveform;
            % Check if channel model is specified between gNB and a
            % particular UE
            if ~isempty(obj.ChannelModel{pktInfo.RNTI})
                rxWaveform = [rxWaveform; zeros(obj.MaxChannelDelay(pktInfo.RNTI), size(rxWaveform,2))];
                rxWaveform = obj.ChannelModel{pktInfo.RNTI}(rxWaveform);
            end
            
            % Apply path loss on the waveform
            [rxWaveform, pathloss] = applyPathLoss(obj, rxWaveform, pktInfo);
            pktInfo.TxPower = pktInfo.TxPower - pathloss;
            
            % Apply receiver antenna gain
            rxWaveform = applyRxGain(obj, rxWaveform);
            pktInfo.TxPower = pktInfo.TxPower + obj.RxGain;
            
            % Add thermal noise to the waveform
            selfInfo.Temperature = obj.Temperature;
            selfInfo.Bandwidth = obj.CarrierInformation.ULBandwidth;
            rxWaveform = applyThermalNoise(obj, rxWaveform, pktInfo, selfInfo);
        end
        
        function phyRxProcessing(obj, rxWaveform, puschInfo)
            %phyRxProcessing Read the PUSCH as per the passed PUSCH information and send the decoded information to MAC
            
            % Get the Tx slot
            if obj.CurrSymbol == 0 % Current symbol is first in the slot hence transmission was done in the last slot
                if obj.CurrSlot > 0
                    txSlot = obj.CurrSlot-1;
                    txSlotAFN = obj.AFN; % Tx slot was in the current frame
                else
                    txSlot = obj.WaveformInfoUL.SlotsPerSubframe*10-1;
                    txSlotAFN = obj.AFN - 1; % Tx slot was in the previous frame
                end
                lastSym = obj.WaveformInfoUL.SymbolsPerSlot-1; % Last symbol number of the waveform
            else % Transmission was done in the current slot
                txSlot = obj.CurrSlot;
                txSlotAFN = obj.AFN; % Tx slot was in the current frame
                lastSym = obj.CurrSymbol - 1; % Last symbol number of the waveform
            end
            startSym = lastSym - puschInfo.PUSCHConfig.SymbolAllocation(2) + 1;
             
            % Carrier information
            carrier = nrCarrierConfig;
            carrier.SubcarrierSpacing = obj.CarrierInformation.SubcarrierSpacing;
            carrier.NSizeGrid = obj.CarrierInformation.NRBsUL;
            carrier.NSlot = txSlot;
            carrier.NFrame = txSlotAFN;
            carrier.NCellID = obj.CellConfig.NCellID;
            
            % Populate the received waveform at appropriate indices in the slot-length waveform
            slotNumSubFrame = mod(txSlot, obj.WaveformInfoUL.SlotsPerSubframe);
            startSymSubframe = slotNumSubFrame*obj.WaveformInfoUL.SymbolsPerSlot + 1; % Start symbol of tx slot in the subframe
            lastSymSubframe = startSymSubframe + obj.WaveformInfoUL.SymbolsPerSlot -1; % Last symbol of tx slot in the subframe
            symbolLengths = obj.WaveformInfoUL.SymbolLengths(startSymSubframe : lastSymSubframe); % Length of symbols of tx slot
            slotWaveform = zeros(sum(symbolLengths) + obj.MaxChannelDelay(puschInfo.PUSCHConfig.RNTI), 1);
            startSample = sum(symbolLengths(1:startSym)) + 1;
            slotWaveform(startSample : startSample+length(rxWaveform)-1) = rxWaveform;
            
            % Get PUSCH and DM-RS information
            [puschIndices, ~] = nrPUSCHIndices(carrier, puschInfo.PUSCHConfig);
            dmrsSymbols = nrPUSCHDMRS(carrier, puschInfo.PUSCHConfig);
            dmrsIndices = nrPUSCHDMRSIndices(carrier, puschInfo.PUSCHConfig);
            
            % Set TBS
            obj.ULSCHDecoders{puschInfo.PUSCHConfig.RNTI}.TransportBlockLength = puschInfo.TBS*8;
            
            % Practical synchronization. Correlate the received waveform
            % with the PUSCH DM-RS to give timing offset estimate 't' and
            % correlation magnitude 'mag'
            [t,mag] = nrTimingEstimate(carrier, slotWaveform, dmrsIndices, dmrsSymbols);
            obj.TimingOffset(puschInfo.PUSCHConfig.RNTI) = hSkipWeakTimingOffset(obj.TimingOffset(puschInfo.PUSCHConfig.RNTI), t, mag);
            offset = obj.TimingOffset(puschInfo.PUSCHConfig.RNTI);
            if(offset > obj.MaxChannelDelay(puschInfo.PUSCHConfig.RNTI))
                % Ignore the timing offset estimate resulting from weak correlation
                offset = 0;
            end
            
            slotWaveform = slotWaveform(1+offset:end, :);
            % Perform OFDM demodulation on the received data to recreate the
            % resource grid, including padding in the event that practical
            % synchronization results in an incomplete slot being demodulated
            rxGrid = nrOFDMDemodulate(carrier, slotWaveform);
            [K, L, R] = size(rxGrid);
            if (L < obj.WaveformInfoUL.SymbolsPerSlot)
                rxGrid = cat(2, rxGrid, zeros(K, obj.WaveformInfoUL.SymbolsPerSlot-L, R));
            end
            
            % Practical channel estimation between the received grid
            % and each transmission layer, using the PUSCH DM-RS for
            % each layer
            [estChannelGrid, noiseEst] = nrChannelEstimate(rxGrid, dmrsIndices, dmrsSymbols);
            
            % Get PUSCH resource elements from the received grid
            [puschRx,puschHest] = nrExtractResources(puschIndices,rxGrid,estChannelGrid);
            
            % Equalization
            [puschEq,csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEst);
            
            % Decode PUSCH physical channel
            [ulschLLRs,rxSymbols] = nrPUSCHDecode(carrier, puschInfo.PUSCHConfig, puschEq, noiseEst);
            
            csi = nrLayerDemap(csi);
            Qm = length(ulschLLRs) / length(rxSymbols);
            csi = reshape(repmat(csi{1}.',Qm,1),[],1);
            ulschLLRs = ulschLLRs .* csi;
            
            % Decode the UL-SCH transport channel
            obj.ULSCHDecoders{puschInfo.PUSCHConfig.RNTI}.TargetCodeRate = puschInfo.TargetCodeRate;
            [decbits, crcFlag] = step(obj.ULSCHDecoders{puschInfo.PUSCHConfig.RNTI}, ulschLLRs, ...
                puschInfo.PUSCHConfig.Modulation, puschInfo.PUSCHConfig.NumLayers, puschInfo.RV, puschInfo.HARQID);
            
            if puschInfo.RV == 1
                % The last redundancy version as per the order [0 3 2 1]
                % failed. Reset the soft buffer
                resetSoftBuffer(obj.ULSCHDecoders{puschInfo.PUSCHConfig.RNTI}, puschInfo.HARQID);
            end
            % Convert bit stream to byte stream
            decbits = (reshape(decbits, 8, []))';
            macPDU = bi2de(decbits);
            
            % Rx callback to MAC
            macPDUInfo = hNRRxIndicationInfo;
            macPDUInfo.RNTI = puschInfo.PUSCHConfig.RNTI;
            macPDUInfo.TBS = puschInfo.TBS;
            macPDUInfo.HARQID = puschInfo.HARQID;
            obj.RxIndicationFcn(macPDU, crcFlag, macPDUInfo); % Send PDU to MAC
            
            % Increment the number of erroneous packets received for UE
            obj.ULBlkErr(puschInfo.PUSCHConfig.RNTI, 1) = obj.ULBlkErr(puschInfo.PUSCHConfig.RNTI, 1) + crcFlag;
            % Increment the number of received packets for UE
            obj.ULBlkErr(puschInfo.PUSCHConfig.RNTI, 2) = obj.ULBlkErr(puschInfo.PUSCHConfig.RNTI, 2) + 1;
            
            if ~isempty(obj.PacketLogger) % Packet capture enabled
                logPackets(obj, puschInfo, macPDU, 1); % Log UL packets
            end
        end
        
        function waveformOut = applyTxPowerLevelAndGain(obj, waverformIn, gain)
            %applyTxPowerLevel Applies Tx power level to IQ samples
            
            % Apply Tx power to IQ samples.
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
            % Calculate the pathloss
            pathloss = fspl(distance, lambda);
            
            %}
            
            %apply pathloss according to 38.901
            %get distance and heights
            gNBPos = obj.Node.NodePosition;
            UEPos = txInfo.Position;
            
            %Calculate new UE position in wrap-around mode
            [~,UEPos] = obj.Node.DistanceCalculatorFcn(UEPos,gNBPos);
            
            % Calculate pathloss
            YUF = YusUtilityFunctions;
            % Construct the structure of the parameters for input for
            % calculating pathloss
            paramForPL.Scenario = obj.YusUtilityParameter.Scenario;
            paramForPL.ULCarrierFreq = obj.ChannelModel{txInfo.RNTI}.CarrierFrequency;
            LinkDir = 1; % 0 for DL, 1 for UL
            UEStat = obj.YusUtilityParameter.UEStat{obj.siteIdx,txInfo.RNTI};
            UEStat.UEPosition = UEPos;
            [pathloss, sigma_SF] = calculatePathloss(YUF, paramForPL, gNBPos, UEStat, LinkDir);
            pathloss = normrnd(pathloss,sigma_SF); % Add shadow fading to pathloss

            %{
            % The functions for calculating pathloss are moved to
            % YusUilityFuctions.m
            d_2D=norm(gNBPos(1:2)-UEPos(1:2));
            d_3D=norm(gNBPos(1:3)-UEPos(1:3));
            h_BS = gNBPos(3);
            h_UT=UEPos(3);
            
            LOS = obj.ChannelModel{txInfo.RNTI}.HasLOSCluster;
    
            %constants
            f_c = obj.ChannelModel{txInfo.RNTI}.CarrierFrequency / 1e9; % frequency is stored in Hz, need in GHz for pathloss calculations
            
            
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
            % Apply pathloss on IQ samples
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
            % Calculate SNR.
            SNR = pktInfo.TxPower - ((10*log10(totalnoise)) + 30);
            % Add noise
            waveformOut = awgn(waveformIn,SNR,pktInfo.TxPower-30, 'db');
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
            
            timestamp = round(obj.getCurrentTime());
            obj.PacketMetaData.HARQID = info.HARQID;
            obj.PacketMetaData.SlotNumber = info.NSlot;
            
            if linkDir % Uplink
                % Get frame number of previous slot i.e the Tx slot. Reception ended at the
                % end of previous slot
                if obj.CurrSlot > 0
                    prevSlotAFN = obj.AFN; % Previous slot was in the current frame
                else
                    % Previous slot was in the previous frame
                    prevSlotAFN = obj.AFN - 1;
                end
                obj.PacketMetaData.SystemFrameNumber = mod(prevSlotAFN, 1024);
                obj.PacketMetaData.LinkDir = obj.PacketLogger.Uplink;
                obj.PacketMetaData.RNTI = info.PUSCHConfig.RNTI;
            else % Downlink
                obj.PacketMetaData.SystemFrameNumber = mod(obj.AFN, 1024);
                obj.PacketMetaData.LinkDir = obj.PacketLogger.Downlink;
                obj.PacketMetaData.RNTI = info.PDSCHConfig.RNTI;
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