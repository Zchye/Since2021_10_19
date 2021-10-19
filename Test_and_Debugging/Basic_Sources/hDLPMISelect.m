function [PMISet,info] = hDLPMISelect(carrier,csirs,reportConfig,nLayers,H,varargin)
%hDLPMISelect PDSCH precoding matrix indicator calculation
%   [PMISET,INFO] = hDLPMISelect(CARRIER,CSIRS,REPORTCONFIG,NLAYERS,H)
%   returns the precoding matrix indicator (PMI) values, as defined in TS
%   38.214 Section 5.2.2.2, for the specified carrier configuration
%   CARRIER, CSI-RS configuration CSIRS, channel state information (CSI)
%   reporting configuration REPORTCONFIG, number of transmission layers
%   NLAYERS, and estimated channel information H.
%
%   CARRIER is a carrier specific configuration object, as described in
%   <a href="matlab:help('nrCarrierConfig')">nrCarrierConfig</a>. Only these object properties are relevant for this
%   function:
%
%   SubcarrierSpacing - Subcarrier spacing in kHz
%   CyclicPrefix      - Cyclic prefix type
%   NSizeGrid         - Number of resource blocks (RBs) in
%                       carrier resource grid
%   NStartGrid        - Start of carrier resource grid relative to common
%                       resource block 0 (CRB 0)
%   NSlot             - Slot number
%   NFrame            - System frame number
%
%   CSIRS is a CSI-RS specific configuration object to specify one or more
%   CSI-RS resources, as described in <a href="matlab:help('nrCSIRSConfig')">nrCSIRSConfig</a>. Only these object
%   properties are relevant for this function:
%
%   CSIRSType           - Type of a CSI-RS resource, whether it is ZP or NZP
%   CSIRSPeriod         - CSI-RS slot periodicity and offset
%   RowNumber           - Row number corresponding to a CSI-RS resource, as
%                         defined in TS 38.211 Table 7.4.1.5.3-1
%   Density             - CSI-RS resource frequency density
%   SymbolLocations     - Time-domain locations of a CSI-RS resource
%   SubcarrierLocations - Frequency-domain locations of a CSI-RS resource
%   NumRB               - Number of RBs allocated for a CSI-RS resource
%   RBOffset            - Starting RB index of CSI-RS allocation relative
%                         to carrier resource grid
%   It is recommended to use the same CSI-RS resource(s) that is used for
%   channel estimate for better results, because the resource elements (REs)
%   not containing the CSI-RS may have the interpolated channel estimates.
%   Note that the CDM lengths and the number of ports configured for all
%   the CSI-RS resources must be same.
%
%   REPORTCONFIG is a CSI reporting configuration structure with these
%   fields:
%   NSizeBWP        - Size of the bandwidth part (BWP) in terms of number
%                     of physical resource blocks (PRBs). It must be a
%                     scalar and the value must be in the range 1...275.
%                     Empty ([]) is also supported and it implies that the
%                     value of NSizeBWP is equal to the size of carrier
%                     resource grid
%   NStartBWP       - Starting PRB index of BWP relative to common resource
%                     block 0 (CRB 0). It must be a scalar and the value
%                     must be in the range 0...2473. Empty ([]) is also
%                     supported and it implies that the value of NStartBWP
%                     is equal to the start of carrier resource grid
%   PanelDimensions - Antenna panel configuration as a two-element vector
%                     in the form of [N1 N2]. N1 represents the number of
%                     antenna elements in horizontal direction and N2
%                     represents the number of antenna elements in vertical
%                     direction. Valid combinations of [N1 N2] are defined
%                     in TS 38.214 Table 5.2.2.2.12. It is not applicable
%                     when the number of CSI-RS ports is less than or equal to 2
%   PMIMode         - Optional. It represents the mode of PMI reporting. It
%                     must be a character array or string scalar. It must
%                     be one of {'Subband', 'Wideband'}. The default value
%                     is 'Wideband'
%   SubbandSize     - Subband size for PMI reporting, provided by the
%                     higher-layer parameter NSBPRB. It must be a positive
%                     scalar and must be one of two possible subband sizes,
%                     as defined in TS 38.214 Table 5.2.1.4-2. It is
%                     applicable only when the PMIMode is provided as
%                     'Subband' and the size of BWP is greater than or
%                     equal to 24 PRBs
%   PRGSize         - Optional. Precoding resource block group (PRG) size
%                     for CQI calculation, provided by the higher-layer
%                     parameter pdsch-BundleSizeForCSI. This field is
%                     applicable when the PMI reporting is needed for the
%                     CSI report quantity cri-RI-i1-CQI. It must be a
%                     scalar and it must be one of {2, 4}. Empty ([]) is
%                     also supported to represent that this field is not
%                     configured by higher layers. If it is present and not
%                     configured as empty, irrespective of the PMIMode,
%                     PRGSize is considered for the number of subbands
%                     calculation instead of SubbandSize and the function
%                     reports PMI for each PRG. The default value is []
%   CodebookMode    - Optional. It represents the codebook mode and it must
%                     be a scalar. The value must be one of {1, 2}. It is
%                     applicable only when the number of transmission
%                     layers is 1 or 2 and number of CSI-RS ports is
%                     greater than 2. The default value is 1
%   CodebookSubsetRestriction
%                   - Optional. It is a binary vector which represents the
%                     codebook subset restriction. When the number of
%                     CSI-RS ports is greater than 2, the length of the
%                     input vector must be N1*N2*O1*O2, where N1 and N2 are
%                     panel configurations obtained from PanelDimensions
%                     field and O1 and O2 are the respective DFT
%                     oversampling factors obtained from TS.38.214 Table
%                     5.2.2.2-2. When the number of CSI-RS ports is 2, the
%                     length of the input vector must be 6, as defined in
%                     TS 38.214 Section 5.2.2.2-1. The default value is
%                     empty ([]), which means there is no codebook subset
%                     restriction
%   i2Restriction   - Optional. It is a binary vector which represents the
%                     restricted i2 values in a codebook. Length of the
%                     input vector must be 16. First element of the input
%                     binary vector corresponds to i2 as 0, second element
%                     corresponds to i2 as 1, and so on. Binary value 1
%                     indicates that the precoding matrix associated with
%                     the respective i2 is unrestricted and 0 indicates
%                     that the precoding matrix associated with the
%                     respective i2 is restricted. This field is applicable
%                     only when the number of CSI-RS ports is greater than
%                     2. For a precoding matrices codebook, if the number
%                     of possible i2 values are less than 16, then only the
%                     required binary elements are considered and the
%                     trailing extra elements in the input vector are
%                     ignored. The default value is empty ([]), which means
%                     there is no i2 restriction
%
%   For the CodebookSubsetRestriction, when the transmission layers is one
%   of {3, 4} and the number of CSI-RS ports is greater than or equal to
%   16, the elements {mod(N2*O2*(2*l-1)+m,N1*O1*N2*O2)+1, N2*O2*(2*l)+m+1,
%   N2*O2*(2*l+1)+m+1} (1-based) are each associated with all the precoding
%   matrices based on vbarlm (l = 0....(N1*O1/2)-1, m = 0....N2*O2), as
%   defined in TS 38.214 Section 5.2.2.2.1. If one or more of the
%   associated binary values is zero, then all the precoding matrices based
%   on vbarlm are restricted. When the number of transmission layers is
%   other than 3, 4 or the number of CSI-RS ports is less than 16, the
%   element N2*O2*l+m+1 (1-based) is associated with all the precoding
%   matrices based on vlm (l = 0...N1*O1-1, m = N2*O2-1). If the associated
%   binary value is zero, then all the precoding matrices based on vlm are
%   restricted.
%
%   NLAYERS is a scalar representing the number of transmission layers. The
%   value must be in the range of 1...8.
%
%   H is the channel estimation matrix. It is of size
%   K-by-L-by-nRxAnts-by-Pcsirs, where K is the number of subcarriers in
%   the carrier resource grid, L is the number of OFDM symbols spanning one
%   slot, nRxAnts is the number of receive antennas, and Pcsirs is the
%   number of CSI-RS antenna ports.
%   Note that the number of transmission layers provided must be less than
%   or equal to min(nRxAnts,Pcsirs).
%
%   PMISET is an output structure with these fields:
%   i1 - Indicates wideband PMI (1-based). It is a three-element vector in
%        the form of [i11 i12 i13]. Note that i13 is not applicable when
%        the number of transmission layers is one of {1, 5, 6, 7, 8}. In
%        that case, function returns the value of i13 as 1
%   i2 - Indicates subband PMI (1-based). When PMIMode is specified as
%        'wideband', it is a scalar representing one i2 indication for the
%        entire band. When PMIMode is specified as 'subband' or when
%        PRGSize is configured as other than empty ([]), one subband
%        indication i2 is reported for each subband or PRG, respectively.
%        Length of the i2 vector in the latter case equals to the number of
%        subbands or PRGs
%   When the number of CSI-RS ports is greater than 2, the precoding matrix
%   is obtained by considering both i1 and i2 values based on TS 38.214
%   Tables 5.2.2.2.1-2 to 5.2.2.2.1-12. When the number of CSI-RS ports is
%   2, the precoding matrix is obtained by a single index (i2 field here)
%   based on TS 38.214 Table 5.2.2.2.1-1. In the latter case, the function
%   returns the i1 as [1 1 1] to support same indexing for all INFO fields
%   for both CSI-RS ports equal to 2 and greater than 2 cases. All the
%   values of i1 and i2 fields are returned as ones when the number of
%   CSI-RS ports is 1.
%
%   INFO is an output structure with these fields:
%   SINRPerRE      - It represents the linear signal to noise plus
%                    interference ratio (SINR) values in each RE within the
%                    BWP for all the layers and all the precoding matrices.
%                    It is a multidimensional array of size
%                       - N-by-L-by-nLayers-by-i2Length-by-i11Length-by-i12Length-by-i13Length
%                         when the number of CSI-RS ports is greater than 2
%                       - N-by-L-by-nLayers-by-i2Length when the number of
%                         CSI-RS ports is 2
%                       - N-by-L when the number of CSI-RS ports is 1
%                    N is the number of subcarriers in the BWP resource
%                    grid, L is the number of OFDM symbols spanning one
%                    slot, nLayers is the number of transmission layers,
%                    i2Length is the maximum number of possible i2 values
%                    and i11Length, i12Length, i13Length are the maximum
%                    number of possible i11, i12 and i13 values for the
%                    given report configuration respectively
%   SINRPerSubband - It represents the linear SINR values in each subband
%                    for all the layers. SINR value in each subband is
%                    formed by averaging SINRPerRE estimates across each
%                    subband (i.e. in the appropriate region of the N
%                    dimension and across the L dimension).
%                    It is a multidimensional array of size
%                       - numSubbands-by-nLayers-by-i2Length-by-i11Length-by-i12Length-by-i13Length
%                         when the number of CSI-RS ports is greater than 2
%                       - numSubbands-by-nLayers-by-i2Length when the
%                         number of CSI-RS ports is 2
%                       - numSubbands-by-1 when the number of CSI-RS ports is 1
%                    The numSubbands represents the number of subbands
%   W              - Multidimensional array containing type-1 single panel
%                    precoding matrices based on the CSI reporting
%                    configuration. It is of size
%                       - Pcsirs-by-nLayers-by-i2Length-by-i11Length-by-i12Length-by-i13Length
%                         when the number of CSI-RS ports is greater than 2
%                       - 2-by-nLayers-by-i2Length when the number of
%                         CSI-RS ports is 2
%                       - 1-by-1 with the value 1 when the number of CSI-RS
%                         ports is 1
%                    The Pcsirs represents the number of CSI-RS ports.
%                    Note that the restricted precoding matrices as per the
%                    report configuration are returned as all zeros
%
%   [PMISET,INFO] = hDLPMISelect(...,NVAR) specifies the estimated noise
%   variance at the receiver NVAR as a nonnegative scalar. By default, the
%   value of nVar is considered as 1e-10, if it is not given as input.
%
%   Note that i1 and i2 fields of PMISET and SINRPerRE and SINRPerSubband
%   fields of INFO are returned as array of NaNs for these cases:
%   - When CSI-RS is not present in the operating slot or in the BWP
%   - When all the precoding matrices in a codebook are restricted
%   Also note that the PMI i2 index is reported as NaNs in the subbands
%   where CSI-RS is not present.
%
%   Note that the function only supports PMI reporting with Type-1 single
%   panel.
%
%   % Example:
%   % This example demonstrates how to calculate PMI.
%
%   % Carrier configuration
%   carrier = nrCarrierConfig;
%
%   % CSI-RS configuration
%   csirs = nrCSIRSConfig;
%   csirs.CSIRSType = {'nzp','nzp'};
%   csirs.RowNumber = [4 4];
%   csirs.Density = {'one','one'};
%   csirs.SubcarrierLocations = {0 0};
%   csirs.SymbolLocations = {0,5};
%   csirs.NumRB = 52;
%   csirs.RBOffset = 0;
%   csirs.CSIRSPeriod = [4 0];
%
%   % Configure the number of transmit and receive antennas
%   nTxAnts = max(csirs.NumCSIRSPorts);
%   nRxAnts = nTxAnts;
%
%   % Configure the number of transmission layers
%   numLayers = 1;
%
%   % Generate CSI-RS indices and symbols
%   csirsInd = nrCSIRSIndices(carrier,csirs);
%   csirsSym = nrCSIRS(carrier,csirs);
%
%   % Resource element mapping
%   txGrid = nrResourceGrid(carrier,nTxAnts);
%   txGrid(csirsInd) = csirsSym;
%
%   % Get OFDM modulation related information
%   OFDMInfo = nrOFDMInfo(carrier);
%
%   % Perform OFDM modulation
%   txWaveform = nrOFDMModulate(carrier,txGrid);
%
%   % Configure the channel parameters.
%   channel = nrTDLChannel;
%   channel.NumTransmitAntennas = nTxAnts;
%   channel.NumReceiveAntennas = nRxAnts;
%   channel.SampleRate = OFDMInfo.SampleRate;
%   channel.DelayProfile = 'TDL-C';
%   channel.DelaySpread = 300e-9;
%   channel.MaximumDopplerShift = 5;
%   chInfo = info(channel);
%
%   % Calculate the maximum channel delay
%   maxChDelay = ceil(max(chInfo.PathDelays*OFDMInfo.SampleRate)) + chInfo.ChannelFilterDelay;
%
%   % Pass the time-domain waveform through the channel
%   rxWaveform = channel([txWaveform; zeros(maxChDelay,nTxAnts)]);
%
%   % Calculate the timing offset
%   offset = nrTimingEstimate(carrier,rxWaveform,csirsInd,csirsSym);
%
%   % Perform timing synchronization
%   rxWaveform = rxWaveform(1+offset:end,:);
%
%   % Add AWGN
%   SNRdB = 20;          % in dB
%   SNR = 10^(SNRdB/20); % Linear value
%   N0 = 1/(sqrt(2.0*channel.NumReceiveAntennas*double(OFDMInfo.Nfft))*SNR); % Noise standard deviation
%   noise = N0*complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
%   rxWaveform = rxWaveform + noise;
%   rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
%
%   % Perform the channel estimate
%   [H,nVar] = nrChannelEstimate(rxGrid,csirsInd,csirsSym,'CDMLengths',[2 1]);
%
%   % Configure the parameters related to CSI reporting
%   reportConfig.NStartBWP = 0;
%   reportConfig.NSizeBWP = 52;
%   reportConfig.PanelDimensions = [2 1];
%   reportConfig.PMIMode = 'Subband';
%   reportConfig.SubbandSize = 4;
%   reportConfig.PRGSize = [];
%   reportConfig.CodebookMode = 2;
%   reportConfig.CodebookSubsetRestriction = [];
%   reportConfig.i2Restriction = [];
%
%   % Calculate the PMI values
%   [PMISet,PMIInfo] = hDLPMISelect(carrier,csirs,reportConfig,numLayers,H,nVar)

%   Copyright 2020 The MathWorks, Inc.

    narginchk(5,6);
    if (nargin == 6)
        nVar = varargin{1};
    else
        % Consider a small noise variance value by default, if the noise
        % variance is not given
        nVar = 1e-10;
    end
    [reportConfig,csirsIndSubs,nVar] = validateInputs(carrier,csirs,reportConfig,nLayers,H,nVar);

    % If PRGSize is present, consider the subband size as PRG size
    if ~isempty(reportConfig.PRGSize)
        reportMode = 'Subband';
        subbandSize = reportConfig.PRGSize;
        ignoreBWPSize = true; % To ignore the BWP size for the validation of PRG size
    else
        reportMode = reportConfig.PMIMode;
        subbandSize = reportConfig.SubbandSize;
        ignoreBWPSize = false; % To consider the BWP size for the validation of subband size
    end
    subbandInfo = getDLCSISubbandInfo(reportConfig.NSizeBWP,reportConfig.NStartBWP,reportMode,subbandSize,ignoreBWPSize);

    numCSIRSPorts = csirs.NumCSIRSPorts(1);
    if numCSIRSPorts == 1
        % W is a scalar with the value 1, when the number of CSI-RS ports
        % is 1
        W = 1;
    else
        % W is a multidimensional matrix of size
        % Pcsirs-by-nLayers-by-i2Length-by-i11Length-by-i12Length-by-i13Length
        % or Pcsirs-by-nLayers-by-i2Length based on the number of CSI-RS
        % ports
        W = getPMIType1SinglePanelCodebook(reportConfig,nLayers);
    end
    % Get the size of W
    [~,~,i2Length,i11Length,i12Length,i13Length] = size(W);

    % Calculate the start of BWP relative to the carrier
    bwpStart = reportConfig.NStartBWP - carrier.NStartGrid;

    % Consider only the RE indices corresponding to the first CSI-RS port
    csirsIndSubs = csirsIndSubs(csirsIndSubs(:,3) == 1,:);
    csirsIndSubs_k = csirsIndSubs(:,1);
    csirsIndSubs_l = csirsIndSubs(:,2);

    % Extract the CSI-RS indices which are present in the BWP
    csirsIndSubs_k = csirsIndSubs_k((csirsIndSubs_k >= bwpStart*12 + 1) & csirsIndSubs_k <= (bwpStart + reportConfig.NSizeBWP)*12);
    csirsIndSubs_l = csirsIndSubs_l((csirsIndSubs_k >= bwpStart*12 + 1) & csirsIndSubs_k <= (bwpStart + reportConfig.NSizeBWP)*12);

    % Make the CSI-RS subscripts relative to BWP
    csirsIndSubs_k = csirsIndSubs_k - bwpStart*12;

    if isempty(csirsIndSubs_k) || ~any(W(:))
        % Report the outputs as all NaNs, if there are no CSI-RS resources
        % present in the BWP or the noise variance value is zero or all the
        % precoding matrices of the codebook are restricted
        PMISet.i1 = [NaN NaN NaN];
        PMISet.i2 = NaN(1,subbandInfo.NumSubbands);
        info.SINRPerRE = NaN(reportConfig.NSizeBWP*12,14,nLayers,i2Length,i11Length,i12Length,i13Length);
        info.SINRPerSubband = NaN(subbandInfo.NumSubbands,nLayers,i2Length,i11Length,i12Length,i13Length);
        info.W = W;
    else
        % Rearrange the channel matrix dimensions from
        % K-by-L-by-nRxAnts-by-Pcsirs to nRxAnts-by-Pcsirs-by-K-by-L
        H = permute(H,[3,4,1,2]);

        SINRPerRE = NaN(reportConfig.NSizeBWP*12,carrier.SymbolsPerSlot,nLayers,i2Length,i11Length,i12Length,i13Length);
        for reIdx = 1:numel(csirsIndSubs_k)
            % Calculate the linear SINR values in all CSI-RS REs for all
            % the layers by considering all unrestricted precoding matrices
            k = csirsIndSubs_k(reIdx);
            l = csirsIndSubs_l(reIdx);
            Htemp = H(:,:,k,l);
            for i11 = 1:i11Length
                for i12 = 1:i12Length
                    for i13 = 1:i13Length
                        for i2 = 1:i2Length
                            currentW = W(:,:,i2,i11,i12,i13);
                            if any(currentW)
                                % Calculate the linear SINR value for the
                                % current RE and for the current precoding
                                % matrix
                                SINRPerRE(k,l,:,i2,i11,i12,i13) = getPrecodedSINR(Htemp,nVar,currentW);
                            end
                        end
                    end
                end
            end
        end

        % Calculate the total SINR value for the entire grid
        % corresponding to each index to check which index gives the
        % maximum wideband SINR value
        if all(isnan(SINRPerRE(:)))
            PMISet.i1 = NaN(1,3);
            PMISet.i2 = NaN;
        else
            totalSINR = squeeze(sum(SINRPerRE,[1 2 3],'omitnan')); % Sum of SINRs across the BWP for all layers for each PMI index
            % Round the total SINR value to four decimals, to avoid the
            % fluctuations in the PMI output because of the minute
            % variations among the SINR values corresponding to different
            % PMI indices
            totalSINR = round(reshape(totalSINR,i2Length,i11Length,i12Length,i13Length),4,'decimal');
            [i2,i11,i12,i13] = ind2sub(size(totalSINR),find(totalSINR == max(totalSINR,[],'all'),1));
            PMISet.i1 = [i11 i12 i13];
            PMISet.i2 = i2;
        end

        % Consider the starting position of the first subband as 0, which
        % is the start of BWP
        subbandStart = 0;
        SubbandSINRs = NaN(subbandInfo.NumSubbands,nLayers,i2Length,i11Length,i12Length,i13Length);
        % Loop over all the subbands
        for SubbandIdx = 1:subbandInfo.NumSubbands
            % Extract the SINR values in the subband
            sinrValuesPerSubband = SINRPerRE((subbandStart*12 + 1):(subbandStart+ subbandInfo.SubbandSizes(SubbandIdx))*12,:,:,:,:,:,:);
            if all(isnan(sinrValuesPerSubband(:))) % CSI-RS is absent in the subband
                % Report i2 as NaN for the current subband as CSI-RS is not
                % present
                PMISet.i2(SubbandIdx) = NaN;
            else                                   % CSI-RS is present in the subband
                % Average the SINR per RE values across the subband for all
                % the PMI indices
                SubbandSINRs(SubbandIdx,:,:,:,:,:) = squeeze(mean(mean(sinrValuesPerSubband,'omitnan'),'omitnan'));

                % Add the subband SINR values across all the layers for
                % each PMI i2 index. Then round the total subband SINR
                % value to four decimals, to avoid the fluctuations in the
                % PMI output
                tempSubbandSINR = round(sum(SubbandSINRs(SubbandIdx,:,:,PMISet.i1(1),PMISet.i1(2),PMISet.i1(3)),2,'omitnan'),4,'decimal');
                % Report i2 index corresponds to the maximum SINR for
                % current subband
                [~,PMISet.i2(SubbandIdx)] = max(tempSubbandSINR);
            end
            % Compute the starting position of next subband
            subbandStart = subbandStart + subbandInfo.SubbandSizes(SubbandIdx);
        end
        SubbandSINRs = reshape(SubbandSINRs,[subbandInfo.NumSubbands,nLayers,i2Length,i11Length,i12Length,i13Length]);

        % Form the output structure
        info.SINRPerRE = SINRPerRE;         % SINR value per RE for all the layers for all PMI indices
        info.SINRPerSubband = SubbandSINRs; % SINR value per subband for all the layers for all PMI indices
        info.W = W;                         % PMI codebook containing the precoding matrices corresponding to all PMI indices
    end
end

function [reportConfigOut,csirsInd,nVar] = validateInputs(carrier,csirs,reportConfig,nLayers,H,nVar)
%   [REPORTCONFIGOUT,CSIRSIND] = validateInputs(CARRIER,CSIRS,REPORTCONFIG,NLAYERS,H,NVAR)
%   validates the inputs arguments and returns the validated CSI report
%   configuration structure REPORTCONFIGOUT along with the NZP-CSI-RS
%   indices CSIRSIND.

    fcnName = 'hDLPMISelect';
    validateattributes(carrier,{'nrCarrierConfig'},{'scalar'},fcnName,'CARRIER');
    % Validate 'csirs'
    validateattributes(csirs,{'nrCSIRSConfig'},{'scalar'},fcnName,'CSIRS');
    if ~isscalar(unique(csirs.NumCSIRSPorts))
        error('nr5g:hDLPMISelect:InvalidCSIRSPorts',...
            'All the CSI-RS resources must be configured to have the same number of CSI-RS ports.');
    end
    if ~iscell(csirs.CDMType)
        cdmType = {csirs.CDMType};
    else
        cdmType = csirs.CDMType;
    end
    if ~all(strcmpi(cdmType,cdmType{1}))
        error('nr5g:hDLPMISelect:InvalidCSIRSCDMTypes',...
            'All the CSI-RS resources must be configured to have the same CDM lengths.');
    end

    % Validate 'reportConfig' 
    % Validate 'NSizeBWP'
    if ~isfield(reportConfig,'NSizeBWP')
        error('nr5g:hDLPMISelect:NSizeBWPMissing','NSizeBWP field is mandatory.');
    end
    nSizeBWP = reportConfig.NSizeBWP;
    if ~(isnumeric(nSizeBWP) && isempty(nSizeBWP))
        validateattributes(nSizeBWP,{'double','single'},{'scalar','integer','positive','<=',275},fcnName,'the size of BWP');
    else
        nSizeBWP = carrier.NSizeGrid;
    end
    % Validate 'NStartBWP'
    if ~isfield(reportConfig,'NStartBWP')
        error('nr5g:hDLPMISelect:NStartBWPMissing','NStartBWP field is mandatory.');
    end
    nStartBWP = reportConfig.NStartBWP;
    if ~(isnumeric(nStartBWP) && isempty(nStartBWP))
        validateattributes(nStartBWP,{'double','single'},{'scalar','integer','nonnegative','<=',2473},fcnName,'the start of BWP');
    else
        nStartBWP = carrier.NStartGrid;
    end
    if nStartBWP < carrier.NStartGrid
        error('nr5g:hDLPMISelect:InvalidNStartBWP',...
            ['The starting resource block of BWP ('...
            num2str(nStartBWP) ') must be greater than '...
            'or equal to the starting resource block of carrier ('...
            num2str(carrier.NStartGrid) ').']);
    end
    % Check whether BWP is located within the limits of carrier or not
    if (nSizeBWP + nStartBWP)>(carrier.NStartGrid + carrier.NSizeGrid)
        error('nr5g:hDLPMISelect:InvalidBWPLimits',['The sum of starting resource '...
            'block of BWP (' num2str(nStartBWP) ') and the size of BWP ('...
            num2str(nSizeBWP) ') must be less than or equal to '...
            'the sum of starting resource block of carrier ('...
            num2str(carrier.NStartGrid) ') and size of the carrier ('...
            num2str(carrier.NSizeGrid) ').']);
    end
    reportConfigOut.NStartBWP = nStartBWP;
    reportConfigOut.NSizeBWP = nSizeBWP;

    % Validate 'PanelDimensions'
    N1 = 1;
    N2 = 1;
    O1 = 1;
    O2 = 1;
    NumCSIRSPorts = csirs.NumCSIRSPorts(1);
    if NumCSIRSPorts > 2
        if ~isfield(reportConfig,'PanelDimensions')
            error('nr5g:hDLPMISelect:PanelDimensionsMissing',...
                'PanelDimensions field is mandatory.');
        end
        validateattributes(reportConfig.PanelDimensions,...
            {'double','single'},{'vector','numel',2},fcnName,'PanelDimensions field');
        N1 = reportConfig.PanelDimensions(1);
        N2 = reportConfig.PanelDimensions(2);
        Pcsirs = 2*prod(reportConfig.PanelDimensions);
        if Pcsirs ~= NumCSIRSPorts
            error('nr5g:hDLPMISelect:InvalidPanelDimensions',...
                ['For the configured number of CSI-RS ports (' num2str(NumCSIRSPorts)...
                '), the given panel configuration [' num2str(N1) ' ' num2str(N2)...
                '] is not valid. Note that, two times the product of panel dimensions ('...
                num2str(2*N1*N2) ') must be equal to the number of CSI-RS ports (' num2str(NumCSIRSPorts) ').']);
        end
        % Supported panel configurations and oversampling factors, as
        % defined in TS 38.214 Table 5.2.2.2.1-2
        panelConfigs = [2     2     4     3     6     4     8     4     6    12     4     8    16   % N1
                        1     2     1     2     1     2     1     3     2     1     4     2     1   % N2
                        4     4     4     4     4     4     4     4     4     4     4     4     4   % O1
                        1     4     1     4     1     4     1     4     4     1     4     4     1]; % O2
        configIdx = find(panelConfigs(1,:) == N1 & panelConfigs(2,:) == N2,1);
        if isempty(configIdx)
            error('nr5g:hDLPMISelect:InvalidPanelConfiguration',['The given panel configuration ['...
                num2str(reportConfig.PanelDimensions(1)) ' ' num2str(reportConfig.PanelDimensions(2)) '] is not valid for the given CSI-RS configuration. '...
                'For a number of CSI-RS ports, the panel configuration should be one of the possibilities from TS 38.214 Table 5.2.2.2.1-2.']);
        end

        % Extract the oversampling factors
        O1 = panelConfigs(3,configIdx);
        O2 = panelConfigs(4,configIdx);
    end
    reportConfigOut.PanelDimensions = [N1 N2];
    reportConfigOut.OverSamplingFactors = [O1 O2];

    % Validate 'PMIMode'
    if isfield(reportConfig,'PMIMode')
        validatestring(reportConfig.PMIMode,{'Wideband','Subband'},fcnName,'PMIMode field');
        reportConfigOut.PMIMode = reportConfig.PMIMode;
    else
        reportConfigOut.PMIMode = 'Wideband';
    end

    % Validate 'SubbandSize'
    NSBPRB = [];
    if strcmpi(reportConfigOut.PMIMode,'Subband')
        if nSizeBWP >= 24
            if ~isfield(reportConfig,'SubbandSize')
                error('nr5g:hDLPMISelect:SubbandSizeMissing',...
                    ['For the subband mode, SubbandSize field is '...
                    'mandatory when the size of BWP is more than 24 PRBs.']);
            end
            validateattributes(reportConfig.SubbandSize,{'double','single'},...
                {'real','scalar'},fcnName,'SubbandSize field');
            NSBPRB = reportConfig.SubbandSize;
        end
    end
    reportConfigOut.SubbandSize = NSBPRB;

    % Validate 'PRGSize'
    if isfield(reportConfig,'PRGSize')
        if ~(isnumeric(reportConfig.PRGSize) && isempty(reportConfig.PRGSize))
            validateattributes(reportConfig.PRGSize,{'double','single'},...
                {'real','scalar'},fcnName,'PRGSize field');
        end
        if ~(isempty(reportConfig.PRGSize) || any(reportConfig.PRGSize == [2 4]))
            error('nr5g:hDLPMISelect:InvalidPRGSize',...
                ['PRGSize value (' num2str(reportConfig.PRGSize) ') must be [], 2, or 4.']);
        end
        reportConfigOut.PRGSize = reportConfig.PRGSize;
    else
        reportConfigOut.PRGSize = [];
    end

    if strcmpi(reportConfigOut.PMIMode,'Subband') && isempty(reportConfigOut.PRGSize)
        if nSizeBWP >= 24
            % Validate the subband size, based on the size of BWP
            % BWP size ranges
            nSizeBWPRange = [24  72;
                             73  144;
                             145 275];
            % Possible values of subband size
            nSBPRBValues = [4  8;
                            8  16;
                            16 32];
            bwpRangeCheck = (nSizeBWP >= nSizeBWPRange(:,1)) & (nSizeBWP <= nSizeBWPRange(:,2));
            validNSBPRBValues = nSBPRBValues(bwpRangeCheck,:);
            if ~any(NSBPRB == validNSBPRBValues)
                error('nr5g:hDLPMISelect:InvalidSubbandSize',['For the configured BWP size (' num2str(nSizeBWP) ...
                    '), subband size (' num2str(NSBPRB) ') must be ' num2str(validNSBPRBValues(1)) ...
                    ' or ' num2str(validNSBPRBValues(2)) '.']);
            end
        end
    end

    % Validate 'CodebookMode'
    if isfield(reportConfig,'CodebookMode')
        validateattributes(reportConfig.CodebookMode,{'numeric'},...
            {'scalar','integer','positive','<=',2},fcnName,'CodebookMode field');
        reportConfigOut.CodebookMode = reportConfig.CodebookMode;
    else
        reportConfigOut.CodebookMode = 1;
    end

    % Validate 'CodebookSubsetRestriction'
    if NumCSIRSPorts > 2
        codebookLength = N1*O1*N2*O2;
        codebookSubsetRestriction = ones(1,codebookLength);
        if isfield(reportConfig,'CodebookSubsetRestriction') &&...
                ~isempty(reportConfig.CodebookSubsetRestriction)
            codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
            validateattributes(codebookSubsetRestriction,...
                {'numeric'},{'vector','binary','numel',codebookLength},fcnName,'CodebookSubsetRestriction field');
        end
    elseif NumCSIRSPorts == 2
        codebookSubsetRestriction = ones(1,6);
        if isfield(reportConfig,'CodebookSubsetRestriction') &&...
                ~isempty(reportConfig.CodebookSubsetRestriction)
            codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
            validateattributes(codebookSubsetRestriction,{'numeric'},{'vector','binary','numel',6},fcnName,'CodebookSubsetRestriction field');
        end
    else
        codebookSubsetRestriction = 1;
    end
    reportConfigOut.CodebookSubsetRestriction = codebookSubsetRestriction;

    % Validate 'i2Restriction'
    i2Restriction = ones(1,16);
    if NumCSIRSPorts > 2
        if isfield(reportConfig,'i2Restriction') &&  ~isempty(reportConfig.i2Restriction)
            validateattributes(reportConfig.i2Restriction,...
                {'numeric'},{'vector','binary','numel',16},fcnName,'i2Restriction field');
            i2Restriction = reportConfig.i2Restriction;
        end
    end
    reportConfigOut.i2Restriction = i2Restriction;

    % Validate 'nLayers'
    validateattributes(nLayers,{'numeric'},{'scalar','integer','positive','<=',8},fcnName,'NLAYERS');

    % Validate 'H'
    validateattributes(H,{'double','single'},{},fcnName,'H');
    validateattributes(numel(size(H)),{'double'},{'>=',2,'<=',4},fcnName,'number of dimensions of H');

    % Ignore zero-power (ZP) CSI-RS resources, as they are not used for CSI
    % estimation
    if ~iscell(csirs.CSIRSType)
        csirs.CSIRSType = {csirs.CSIRSType};
    end
    numZPCSIRSRes = sum(strcmpi(csirs.CSIRSType,'zp'));
    tempInd = nrCSIRSIndices(carrier,csirs,"IndexStyle","subscript","OutputResourceFormat","cell");
    tempInd = tempInd(numZPCSIRSRes+1:end)';
    csirsInd = zeros(0,3);
    if ~isempty(tempInd)
        csirsInd = cell2mat(tempInd);
    end
    if ~isempty(csirsInd)
        K = carrier.NSizeGrid*12;
        L = carrier.SymbolsPerSlot;
        validateattributes(H,{class(H)},{'size',[K L NaN NumCSIRSPorts]},fcnName,'H');

        % Validate 'nLayers'
        nRxAnts = size(H,3);
        maxNLayers = min(nRxAnts,NumCSIRSPorts);
        if nLayers > maxNLayers
            error('nr5g:hDLPMISelect:InvalidNumLayers',...
                ['The given antenna configuration (' ...
                num2str(NumCSIRSPorts) 'x' num2str(nRxAnts)...
                ') supports only up to (' num2str(maxNLayers) ') layers.']);
        end
    end

    % Validate 'nVar'
    validateattributes(nVar,{'double','single'},{'scalar','real','nonnegative','finite'},fcnName,'NVAR');
    % Clip nVar to a small noise variance to avoid +/-Inf outputs
    if nVar < 1e-10
        nVar = 1e-10;
    end
end

function  W = getPMIType1SinglePanelCodebook(reportConfig,nLayers)
%   W = getPMIType1SinglePanelCodebook(REPORTCONFIG,NLAYERS) returns type-1
%   single panel precoding matrices W, as defined in TS 38.214 Tables
%   5.2.2.2.1-1 to 5.2.2.2.1-12 by considering these inputs:
%
%   REPORTCONFIG is a CSI reporting configuration structure with these
%   fields:
%   PanelDimensions            - Antenna panel configuration as a
%                                two-element vector ([N1 N2]). It is
%                                not applicable for CSI-RS ports less
%                                than or equal to 2
%   OverSamplingFactors        - DFT oversampling factors corresponds to
%                                the panel configuration
%   CodebookMode               - Codebook mode. Applicable only when the
%                                number of MIMO layers is 1 or 2 and
%                                number of CSI-RS ports is greater than 2
%   CodebookSubsetRestriction  - Binary vector for vlm or vbarlm restriction
%   i2Restriction              - Binary vector for i2 restriction
%
%   NLAYERS      - Number of transmission layers
%
%   W            - Multidimensional array containing unrestricted type-1
%                  single panel precoding matrices. It is of size
%                  Pcsirs-by-nLayers-by-i2Length-by-i11Length-by-i12Length-by-i13Length
%
%   Note that the restricted precoding matrices are returned as all zeros.

    panelDimensions           = reportConfig.PanelDimensions;
    codebookMode              = reportConfig.CodebookMode;
    codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
    i2Restriction             = reportConfig.i2Restriction;

    % Create a function handle to compute the co-phasing factor value
    % according to TS 38.214 Section 5.2.2.2, considering the co-phasing
    % factor index
    phi = @(x)exp(1i*pi*x/2);

    % Get the number of CSI-RS ports using the panel dimensions
    Pcsirs = 2*panelDimensions(1)*panelDimensions(2);
    if Pcsirs == 2
        % Codebooks for 1-layer and 2-layer CSI reporting using antenna
        % ports 3000 to 3001, as defined in TS 38.214 Table 5.2.2.2.1-1
        if nLayers == 1
            W(:,:,1) = 1/sqrt(2).*[1; 1];
            W(:,:,2) = 1/sqrt(2).*[1; 1i];
            W(:,:,3) = 1/sqrt(2).*[1; -1];
            W(:,:,4) = 1/sqrt(2).*[1; -1i];
            restrictedIndices = find(~codebookSubsetRestriction);
            restrictedIndices = restrictedIndices(restrictedIndices <= 4);
            if ~isempty(restrictedIndices)
                restrictedSet = logical(sum(restrictedIndices == [1;2;3;4],2));
                W(:,:,restrictedSet) = 0;
            end
        elseif nLayers == 2
            W(:,:,1) = 1/2*[1 1;1 -1];
            W(:,:,2) = 1/2*[1 1; 1i -1i];
            restrictedIndices = find(~codebookSubsetRestriction);
            restrictedIndices = restrictedIndices(restrictedIndices > 4);
            if ~isempty(restrictedIndices)
                restrictedSet = logical(sum(restrictedIndices == [5;6],2));
                W(:,:,restrictedSet) = 0;
            end
        end
    elseif Pcsirs > 2
        N1 = panelDimensions(1);
        N2 = panelDimensions(2);
        O1 = reportConfig.OverSamplingFactors(1);
        O2 = reportConfig.OverSamplingFactors(2);

        % Select the codebook based on the number of layers, panel
        % configuration, and the codebook mode
        switch nLayers
            case 1 % Number of layers is 1
                % Codebooks for 1-layer CSI reporting using antenna ports
                % 3000 to 2999+P_CSIRS, as defined in TS 38.214 Table
                % 5.2.2.2.1-5
                if codebookMode == 1
                    i11_length = N1*O1;
                    i12_length = N2*O2;
                    i2_length = 4;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length);
                    % Loop over all the values of i11, i12, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i2 = 0:i2_length-1
                                l = i11;
                                m = i12;
                                n = i2;
                                bitIndex = N2*O2*l+m;
                                [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                                if ~(lmRestricted || i2Restricted)
                                    vlm = getVlm(N1,N2,O1,O2,l,m);
                                    phi_n = phi(n);
                                    W(:,:,i2+1,i11+1,i12+1) = (1/sqrt(Pcsirs))*[vlm ;...
                                                                                phi_n*vlm];
                                end
                            end
                        end
                    end
                else % codebookMode == 2
                    i11_length = N1*O1/2;
                    i12_length = N2*O2/2;
                    if N2 == 1
                        i12_length = 1;
                    end
                    i2_length = 16;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length);
                    % Loop over all the values of i11, i12, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i2 = 0:i2_length-1
                                floor_i2by4 = floor(i2/4);
                                if N2 == 1
                                    l = 2*i11 + floor_i2by4;
                                    m = 0;
                                else % N2 > 1
                                    lmAddVals = [0 0; 1 0; 0 1;1 1];
                                    l = 2*i11 + lmAddVals(floor_i2by4+1,1);
                                    m = 2*i12 + lmAddVals(floor_i2by4+1,2);
                                end
                                n = mod(i2,4);
                                bitIndex = N2*O2*l+m;
                                [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                                if ~(lmRestricted || i2Restricted)
                                    vlm = getVlm(N1,N2,O1,O2,l,m);
                                    phi_n = phi(n);
                                    W(:,:,i2+1,i11+1,i12+1) = (1/sqrt(Pcsirs))*[vlm;...
                                                                                phi_n*vlm];
                                end
                            end
                        end
                    end
                end

            case 2 % Number of layers is 2
                % Codebooks for 2-layer CSI reporting using antenna ports
                % 3000 to 2999+P_CSIRS, as defined in TS 38.214 Table
                % 5.2.2.2.1-6

                % Compute i13 parameter range and corresponding k1 and k2,
                % as defined in TS 38.214 table 5.2.2.2.1-3
                if (N1 > N2) && (N2 > 1)
                    i13_length = 4;
                    k1 = [0 O1 0 2*O1];
                    k2 = [0 0 O2 0];
                elseif N1 == N2
                    i13_length = 4;
                    k1 = [0 O1 0 O1];
                    k2 = [0 0 O2 O2];
                elseif (N1 == 2) && (N2 == 1)
                    i13_length = 2;
                    k1 = O1*(0:1);
                    k2 = [0 0];
                else
                    i13_length = 4;
                    k1 = O1*(0:3);
                    k2 = [0 0 0 0] ;
                end

                if codebookMode == 1
                    i11_length = N1*O1;
                    i12_length = N2*O2;
                    i2_length = 2;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length,i13_length);
                    % Loop over all the values of i11, i12, i13, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    l = i11;
                                    m = i12;
                                    n = i2;
                                    lPrime = i11+k1(i13+1);
                                    mPrime = i12+k2(i13+1);
                                    bitIndex = N2*O2*l+m;
                                    [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                                    if ~(lmRestricted || i2Restricted)
                                        vlm = getVlm(N1,N2,O1,O2,l,m);
                                        vlPrime_mPrime = getVlm(N1,N2,O1,O2,lPrime,mPrime);
                                        phi_n = phi(n);
                                        W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(2*Pcsirs))*[vlm        vlPrime_mPrime;...
                                                                phi_n*vlm  -phi_n*vlPrime_mPrime];
                                    end
                                end
                            end
                        end
                    end
                else % codebookMode == 2
                    i11_length = N1*O1/2;
                    if N2 == 1
                        i12_length = 1;
                    else
                        i12_length = N2*O2/2;
                    end
                    i2_length = 8;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length,i13_length);
                    % Loop over all the values of i11, i12, i13, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    floor_i2by2 = floor(i2/2);
                                    if N2 == 1
                                        l = 2*i11 + floor_i2by2;
                                        lPrime = 2*i11 + floor_i2by2 + k1(i13+1);
                                        m = 0;
                                        mPrime = 0;
                                    else % N2 > 1
                                        lmAddVals = [0 0; 1 0; 0 1;1 1];
                                        l = 2*i11 + lmAddVals(floor_i2by2+1,1);
                                        lPrime =  2*i11 + k1(i13+1) + lmAddVals(floor_i2by2+1,1);
                                        m = 2*i12 + lmAddVals(floor_i2by2+1,2);
                                        mPrime =  2*i12 + k2(i13+1) + lmAddVals(floor_i2by2+1,2);
                                    end
                                    n = mod(i2,2);
                                    bitIndex = N2*O2*l+m;
                                    [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                                    if ~(lmRestricted || i2Restricted)
                                        vlm = getVlm(N1,N2,O1,O2,l,m);
                                        vlPrime_mPrime = getVlm(N1,N2,O1,O2,lPrime,mPrime);
                                        phi_n = phi(n);
                                        W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                            (1/sqrt(2*Pcsirs))*[vlm        vlPrime_mPrime;...
                                                                phi_n*vlm  -phi_n*vlPrime_mPrime];
                                    end
                                end
                            end
                        end
                    end
                end

            case {3,4} % Number of layers is 3 or 4
                if (Pcsirs < 16)
                    % For the number of CSI-RS ports less than 16, compute
                    % i13 parameter range, corresponding k1 and k2,
                    % according to TS 38.214 table 5.2.2.2.1-4
                    if (N1 == 2) && (N2 == 1)
                        i13_length = 1;
                        k1 = O1;
                        k2 = 0;
                    elseif (N1 == 4) && (N2 == 1)
                        i13_length = 3;
                        k1 = O1*(1:3);
                        k2 = [0 0 0];
                    elseif (N1 == 6) && (N2 == 1)
                        i13_length = 4;
                        k1 = O1*(1:4);
                        k2 = [0 0 0 0];
                    elseif (N1 == 2) && (N2 == 2)
                        i13_length = 3;
                        k1 = [O1 0 O1];
                        k2 = [0 O2 O2];
                    elseif (N1 == 3) && (N2 == 2)
                        i13_length = 4;
                        k1 = [O1 0 O1 2*O1];
                        k2 = [0 O2 O2 0];
                    end

                    % For 3 and 4 layers the procedure for computation of W
                    % is same, other than the dimensions of W. Compute W
                    % for either case accordingly
                    i11_length = N1*O1;
                    i12_length = N2*O2;
                    i2_length = 2;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length,i13_length);
                    % Loop over all the values of i11, i12, i13, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    l = i11;
                                    lPrime = i11+k1(i13+1);
                                    m = i12;
                                    mPrime = i12+k2(i13+1);
                                    n = i2;
                                    bitIndex = N2*O2*l+m;
                                    [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                                    if ~(lmRestricted || i2Restricted)
                                        vlm = getVlm(N1,N2,O1,O2,l,m);
                                        vlPrime_mPrime = getVlm(N1,N2,O1,O2,lPrime,mPrime);
                                        phi_n = phi(n);
                                        phi_vlm = phi_n*vlm;
                                        phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                                        if nLayers == 3
                                            % Codebook for 3-layer CSI
                                            % reporting using antenna ports
                                            % 3000 to 2999+P_CSIRS, as
                                            % defined in TS 38.214 Table
                                            % 5.2.2.2.1-7
                                            W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                                (1/sqrt(3*Pcsirs))*[vlm      vlPrime_mPrime      vlm;...
                                                                    phi_vlm  phi_vlPrime_mPrime  -phi_vlm];
                                        else
                                            % Codebook for 4-layer CSI
                                            % reporting using antenna ports
                                            % 3000 to 2999+P_CSIRS, as
                                            % defined in TS 38.214 Table
                                            % 5.2.2.2.1-8
                                            W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                                (1/sqrt(4*Pcsirs))*[vlm      vlPrime_mPrime      vlm       vlPrime_mPrime;...
                                                                    phi_vlm  phi_vlPrime_mPrime  -phi_vlm  -phi_vlPrime_mPrime];
                                        end
                                    end
                                end
                            end
                        end
                    end
                else % Number of CSI-RS ports is greater than or equal to 16
                    i11_length = N1*O1/2;
                    i12_length = N2*O2;
                    i13_length = 4;
                    i2_length = 2;
                    W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length,i13_length);
                    % Loop over all the values of i11, i12, i13, and i2
                    for i11 = 0:i11_length-1
                        for i12 = 0:i12_length-1
                            for i13 = 0:i13_length-1
                                for i2 = 0:i2_length-1
                                    theta = exp(1i*pi*i13/4);
                                    l = i11;
                                    m = i12;
                                    n = i2;
                                    phi_n = phi(n);
                                    bitValues = [mod(N2*O2*(2*l-1)+m,N1*O1*N2*O2), N2*O2*(2*l)+m, N2*O2*(2*l+1)+m];
                                    [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitValues,i2,i2Restriction);
                                    if ~(lmRestricted || i2Restricted)
                                        vbarlm = getVbarlm(N1,N2,O1,O2,l,m);
                                        theta_vbarlm = theta*vbarlm;
                                        phi_vbarlm = phi_n*vbarlm;
                                        phi_theta_vbarlm = phi_n*theta*vbarlm;
                                        if nLayers == 3
                                            % Codebook for 3-layer CSI
                                            % reporting using antenna ports
                                            % 3000 to 2999+P_CSIRS, as
                                            % defined in TS 38.214 Table
                                            % 5.2.2.2.1-7
                                            W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                                (1/sqrt(3*Pcsirs))*[vbarlm            vbarlm             vbarlm;...
                                                                    theta_vbarlm      -theta_vbarlm      theta_vbarlm;...
                                                                    phi_vbarlm        phi_vbarlm         -phi_vbarlm;...
                                                                    phi_theta_vbarlm  -phi_theta_vbarlm  -phi_theta_vbarlm];
                                        else
                                            % Codebook for 4-layer CSI
                                            % reporting using antenna ports
                                            % 3000 to 2999+P_CSIRS, as
                                            % defined in TS 38.214 Table
                                            % 5.2.2.2.1-8
                                            W(:,:,i2+1,i11+1,i12+1,i13+1) = ...
                                                (1/sqrt(4*Pcsirs))*[vbarlm            vbarlm             vbarlm             vbarlm;...
                                                                    theta_vbarlm      -theta_vbarlm      theta_vbarlm       -theta_vbarlm;...
                                                                    phi_vbarlm        phi_vbarlm         -phi_vbarlm        -phi_vbarlm;...
                                                                    phi_theta_vbarlm  -phi_theta_vbarlm  -phi_theta_vbarlm  phi_theta_vbarlm];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end

            case {5,6} % Number of layers is 5 or 6
                i11_length = N1*O1;
                if N2 == 1
                    i12_length = 1;
                else % N2 > 1
                    i12_length = N2*O2;
                end
                i2_length = 2;
                W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length);
                % Loop over all the values of i11, i12, and i2
                for i11 = 0:i11_length-1
                    for i12 = 0:i12_length-1
                        for i2 = 0:i2_length-1
                            if N2 == 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+2*O1;
                                m = 0;
                                mPrime = 0;
                                m_dPrime = 0;
                            else % N2 > 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+O1;
                                m = i12;
                                mPrime = i12;
                                m_dPrime = i12+O2;
                            end
                            n = i2;
                            bitIndex = N2*O2*l+m;
                            [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                            if ~(lmRestricted || i2Restricted)
                                vlm = getVlm(N1,N2,O1,O2,l,m);
                                vlPrime_mPrime = getVlm(N1,N2,O1,O2,lPrime,mPrime);
                                vlDPrime_mDPrime = getVlm(N1,N2,O1,O2,l_dPrime,m_dPrime);
                                phi_n = phi(n);
                                phi_vlm = phi_n*vlm;
                                phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                                if nLayers == 5
                                    % Codebook for 5-layer CSI reporting
                                    % using antenna ports 3000 to
                                    % 2999+P_CSIRS, as defined in TS 38.214
                                    % Table 5.2.2.2.1-9
                                    W(:,:,i2+1,i11+1,i12+1) = ...
                                        1/(sqrt(5*Pcsirs))*[vlm       vlm        vlPrime_mPrime   vlPrime_mPrime    vlDPrime_mDPrime;...
                                                            phi_vlm   -phi_vlm   vlPrime_mPrime   -vlPrime_mPrime   vlDPrime_mDPrime];
                                else
                                    % Codebook for 6-layer CSI reporting
                                    % using antenna ports 3000 to
                                    % 2999+P_CSIRS, as defined in TS 38.214
                                    % Table 5.2.2.2.1-10
                                    W(:,:,i2+1,i11+1,i12+1) = ...
                                        1/(sqrt(6*Pcsirs))*[vlm       vlm        vlPrime_mPrime       vlPrime_mPrime        vlDPrime_mDPrime   vlDPrime_mDPrime;...
                                                            phi_vlm   -phi_vlm   phi_vlPrime_mPrime   -phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime];
                                end
                            end
                        end
                    end
                end

            case{7,8} % Number of layers is 7 or 8
                if N2 == 1
                    i12_length = 1;
                    if N1 == 4
                        i11_length = N1*O1/2;
                    else % N1 > 4
                        i11_length = N1*O1;
                    end
                else % N2 > 1
                    i11_length = N1*O1;
                    if (N1 == 2 && N2 == 2) || (N1 > 2 && N2 > 2)
                        i12_length = N2*O2;
                    else % (N1 > 2 && N2 == 2)
                        i12_length = N2*O2/2;
                    end
                end
                i2_length = 2;
                W = zeros(Pcsirs,nLayers,i2_length,i11_length,i12_length);
                % Loop over all the values of i11, i12, and i2
                for i11 = 0:i11_length-1
                    for i12 = 0:i12_length-1
                        for i2 = 0:i2_length-1
                            if N2 == 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11+2*O1;
                                l_tPrime = i11+3*O1;
                                m = 0;
                                mPrime = 0;
                                m_dPrime = 0;
                                m_tPrime = 0;
                            else % N2 > 1
                                l = i11;
                                lPrime = i11+O1;
                                l_dPrime = i11;
                                l_tPrime = i11+O1;
                                m = i12;
                                mPrime = i12;
                                m_dPrime = i12+O2;
                                m_tPrime = i12+O2;
                            end
                            n = i2;
                            bitIndex = N2*O2*l+m;
                            [lmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,i2,i2Restriction);
                            if ~(lmRestricted || i2Restricted)
                                vlm = getVlm(N1,N2,O1,O2,l,m);
                                vlPrime_mPrime = getVlm(N1,N2,O1,O2,lPrime,mPrime);
                                vlDPrime_mDPrime = getVlm(N1,N2,O1,O2,l_dPrime,m_dPrime);
                                vlTPrime_mTPrime = getVlm(N1,N2,O1,O2,l_tPrime,m_tPrime);
                                phi_n = phi(n);
                                phi_vlm = phi_n*vlm;
                                phi_vlPrime_mPrime = phi_n*vlPrime_mPrime;
                                if nLayers == 7
                                    % Codebook for 7-layer CSI reporting
                                    % using antenna ports 3000 to
                                    % 2999+P_CSIRS, as defined in TS 38.214
                                    % Table 5.2.2.2.1-11
                                    W(:,:,i2+1,i11+1,i12+1) = ...
                                        1/(sqrt(7*Pcsirs))*[vlm       vlm        vlPrime_mPrime       vlDPrime_mDPrime   vlDPrime_mDPrime    vlTPrime_mTPrime   vlTPrime_mTPrime;...
                                                            phi_vlm   -phi_vlm   phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime   vlTPrime_mTPrime   -vlTPrime_mTPrime];
                                else
                                    % Codebook for 8-layer CSI reporting
                                    % using antenna ports 3000 to
                                    % 2999+P_CSIRS, as defined in TS 38.214
                                    % Table 5.2.2.2.1-12
                                    W(:,:,i2+1,i11+1,i12+1) = ...
                                        1/(sqrt(8*Pcsirs))*[vlm       vlm        vlPrime_mPrime       vlPrime_mPrime        vlDPrime_mDPrime   vlDPrime_mDPrime    vlTPrime_mTPrime   vlTPrime_mTPrime;...
                                                            phi_vlm   -phi_vlm   phi_vlPrime_mPrime   -phi_vlPrime_mPrime   vlDPrime_mDPrime   -vlDPrime_mDPrime   vlTPrime_mTPrime   -vlTPrime_mTPrime];
                                end
                            end
                        end
                    end
                end
        end
    end
end

function vlm = getVlm(N1,N2,O1,O2,l,m)
%   VLM = getVlm(N1,N2,O1,O2,L,M) computes vlm vector according to
%   TS 38.214 Section 5.2.2.2 considering the panel configuration [N1, N2],
%   DFT oversampling factors [O1, O2], and vlm indices L and M.

    um = exp(2*pi*1i*m*(0:N2-1)/(O2*N2));
    ul = exp(2*pi*1i*l*(0:N1-1)/(O1*N1)).';
    vlm =  reshape((ul.*um).',[],1);
end

function vbarlm = getVbarlm(N1,N2,O1,O2,l,m)
%   VBARLM = getVbarlm(N1,N2,O1,O2,L,M) computes vbarlm vector according to
%   TS 38.214 Section 5.2.2.2 considering the panel configuration [N1, N2],
%   DFT oversampling factors [O1, O2], and vbarlm indices L and M.

    % Calculate vbarlm (DFT vector required to compute the precoding matrix)
    um = exp(2*pi*1i*m*(0:N2-1)/(O2*N2));
    ul = exp(2*pi*1i*l*(0:(N1/2)-1)/(O1*N1/2)).';
    vbarlm = reshape((ul.*um).',[],1);
end

function [vlmRestricted,i2Restricted] = isRestricted(codebookSubsetRestriction,bitIndex,n,i2Restriction)
%   [VLMRESTRICTED,I2RESTRICTED] = isRestricted(CODEBOOKSUBSETRESTRICTION,BITINDEX,N,I2RESTRICTION)
%   returns the status of vlm or vbarlm restriction and i2 restriction for
%   a codebook index set, as defined in TS 38.214 Section 5.2.2.2.1 by
%   considering these inputs:
%
%   CODEBOOKSUBSETRESTRICTION - Binary vector for vlm or vbarlm restriction
%   BITINDEX                  - Bit index or indices (0-based) associated
%                               with all the precoding matrices based on
%                               vlm or vbarlm
%   N                         - Co-phasing factor index
%   I2RESTRICTION             - Binary vector for i2 restriction

    % Get the restricted index positions from the codebookSubsetRestriction
    % binary vector
    restrictedIdx = reshape(find(~codebookSubsetRestriction)-1,1,[]);
    vlmRestricted = false;
    if any(sum(restrictedIdx == bitIndex(:),2))
        vlmRestricted = true;
    end

    restrictedi2List = find(~i2Restriction)-1;
    i2Restricted = false;
    % Update the i2Restricted flag, if the precoding matrices based on vlm
    % or vbarlm are restricted
    if any(restrictedi2List == n)
        i2Restricted = true;
    end
end

function info = getDLCSISubbandInfo(nSizeBWP,nStartBWP,reportingMode,NSBPRB,ignoreBWPSize)
%   INFO = getDLCSISubbandInfo(NSIZEBWP,NSTARTBWP,REPORTINGMODE,NSBPRB,FLAG)
%   returns the subband or PRG information INFO, as defined in TS 38.214
%   Table 5.2.1.4-2 by considering these inputs:
%
%   NSIZEBWP      - Size of BWP in terms of number of PRBs
%   NSTARTBWP     - Starting PRB index of BWP relative to CRB 0
%   REPORTINGMODE - PMI reporting mode
%   NSBPRB        - Subband size
%   IGNOREBWPSIZE - Logical flag to consider (false) or ignore (true) the
%                   BWP size for the calculation of subband or PRG sizes

    % Get the subband information
    if strcmpi(reportingMode,'Wideband') || (~ignoreBWPSize && nSizeBWP < 24)
        % According to TS 38.214 Table 5.2.1.4-2, if the size of BWP is
        % less than 24 PRBs, the division of BWP into subbands is not
        % applicable. In this case, the number of subbands is considered as
        % 1 and the subband size is considered as the size of BWP
        numSubbands = 1;
        NSBPRB = nSizeBWP;
        subbandSizes = NSBPRB;
    else
        % Calculate the size of first subband
        firstSubbandSize = NSBPRB - mod(nStartBWP,NSBPRB);

        % Calculate the size of last subband
        if mod(nStartBWP + nSizeBWP,NSBPRB) ~= 0
            lastSubbandSize = mod(nStartBWP + nSizeBWP,NSBPRB);
        else
            lastSubbandSize = NSBPRB;
        end

        % Calculate the number of subbands
        numSubbands = (nSizeBWP - (firstSubbandSize + lastSubbandSize))/NSBPRB + 2;

        % Form a vector with each element representing the size of a subband
        subbandSizes = NSBPRB*ones(1,numSubbands);
        subbandSizes(1) = firstSubbandSize;
        subbandSizes(end) = lastSubbandSize;
    end
    % Place the number of subbands and subband sizes in the output
    % structure
    info.NumSubbands = numSubbands;
    info.SubbandSizes = subbandSizes;
end

function sinr = getPrecodedSINR(H,nVar,W)
%   SINR = getPrecodedSINR(H,NVAR,W) returns the linear SINR values for an
%   RE across all the layers by considering the channel matrix H, estimated
%   noise variance NVAR, and precoding matrix W.

    % Calculate the SINR values as per LMMSE method
    noise = nVar*eye(size(W,2)); % Noise variance multiplied by identity matrix
    den = noise/((W'*H')*H*W+noise);
    sinr = real(((1./diag(den))-1));
end