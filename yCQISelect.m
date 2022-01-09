function [CQI,PMISet,CQIInfo,PMIInfo] = yCQISelect(carrier,varargin)
    
    narginchk(5,7)
    % Extract the input arguments 
    [reportConfig,nLayers,H,nVar,SINRTable,isCSIRSObjSyntax,nTxAnts,csirsInd,out] = parseInputs(carrier,varargin);
    if ~isempty(out)
        % If the out is not empty, follow the syntax with CSI-RS
        % configuration object
        csirs = out;
    end

    % Validate the input arguments
    [reportConfig,SINRTable,nVar] = validateInputs(carrier,reportConfig,nLayers,H,nVar,SINRTable,nTxAnts,csirsInd,isCSIRSObjSyntax);
 
    % Calculate the number of subbands and size of each subband for the
    % given CQI configuration
    CQISubbandInfo = getDLCSISubbandInfo(reportConfig.NSizeBWP,reportConfig.NStartBWP,reportConfig.CQIMode,reportConfig.SubbandSize,false);
    % Calculate the number of subbands and size of each subband for the
    % given PMI configuration. The parameters are useful in deducing the
    % SINR values in RB level granularity when PMI Mode is 'Subband'
    PMISubbandInfo = getDLCSISubbandInfo(reportConfig.NSizeBWP,reportConfig.NStartBWP,reportConfig.PMIMode,reportConfig.SubbandSize,false);

    if isfield(reportConfig,'PRGSize') && ~isempty(reportConfig.PRGSize)
        % Compute the PRG related information, if PRGSize parameter is
        % present and configured with a value other than empty
        PRGInfo = getDLCSISubbandInfo(reportConfig.NSizeBWP,reportConfig.NStartBWP,'Subband',reportConfig.PRGSize,true);
    end

    % Calculate the number of codewords for the given number of layers. For
    % number of layers greater than 4, there are two codewords, else one
    % codeword
    numCodewords = ceil(nLayers/4);

    if ~isCSIRSObjSyntax
        % Calculate the SINR and CQI values according to the syntax with
        % CSI-RS indices. It supports the computation of the CQI values for
        % SISO case

        % Consider W as 1 in SISO case
        W = 1;

        % Calculate the start of BWP relative to the carrier
        bwpStart = reportConfig.NStartBWP - carrier.NStartGrid;
        % Consider only the unique positions for all CSI-RS ports to avoid
        % repetitive calculation
        csirsInd = unique(csirsInd);

        % Convert CSI-RS indices to subscripts in 1-based notation
        [csirsIndSubs_k,csirsIndSubs_l,~] = ind2sub([carrier.NSizeGrid*12 carrier.SymbolsPerSlot nTxAnts],csirsInd);
        % Consider the CSI-RS indices present only in the BWP
        csirsIndSubs_k = csirsIndSubs_k((csirsIndSubs_k >= bwpStart*12 + 1) & csirsIndSubs_k <= (bwpStart + reportConfig.NSizeBWP)*12);
        csirsIndSubs_l = csirsIndSubs_l((csirsIndSubs_k >= bwpStart*12 + 1) & csirsIndSubs_k <= (bwpStart + reportConfig.NSizeBWP)*12);
        % Make the CSI-RS subscripts relative to BWP start
        csirsIndSubs_k = csirsIndSubs_k - bwpStart*12;
        if isempty(csirsIndSubs_k) || (nVar == 0)
            % Report PMI related outputs as all NaNs, if there are no
            % CSI-RS resources present in the BWP or the noise variance
            % value is zero
            PMISet.i1 = [NaN NaN NaN];
            PMISet.i2 = NaN(1,PMISubbandInfo.NumSubbands);

            PMIInfo.SINRPerSubband = NaN(PMISubbandInfo.NumSubbands,nLayers);
            PMIInfo.SINRPerRE = NaN(reportConfig.NSizeBWP*12,carrier.SymbolsPerSlot,nLayers);
            PMIInfo.W = W;
        else
            sigma = sqrt(nVar);
            K = reportConfig.NSizeBWP*12;
            L = carrier.SymbolsPerSlot;
            % Create an SINR grid of NaNs to calculate SINR for each RE
            SINRsperRE = NaN(K,L);
            % Loop over all the REs in which CSI-RS is present
            for reIdx = 1: numel(csirsIndSubs_k)
                prgIdx = csirsIndSubs_k(reIdx);
                l = csirsIndSubs_l(reIdx);
                Htemp = H(prgIdx,l);
                % Compute the SINR values at each subcarrier location where
                % CSI-RS is present
                SINRsperRE(prgIdx,l) = yPrecodedSINR(Htemp,sigma,W);
            end

            % Consider the PMI indices as all ones for SISO case
            PMI.i1 = [1 1 1];
            PMI.i2 = 1*ones(1,CQISubbandInfo.NumSubbands);

            % Compute the SINR values in subband level granularity
            % according to CQI mode
            SINRperSubbandperCW = getSubbandSINR(SINRsperRE,PMI,CQISubbandInfo); % Corresponds to single codeword

            % Compute wideband SINR as a mean of subband SINR values and
            % place it in position 1
            SINRperSubbandperCW = [entropicMean(SINRperSubbandperCW); SINRperSubbandperCW];

            % Get the SINR value per RB spanning one slot
            SINRsperRBperCW = getSINRperRB(SINRsperRE,PMI,CQISubbandInfo.SubbandSizes);

            % This syntax does not consider the PMI mode. The PMISet and
            % PMIInfo output are returned by considering the PMI mode as
            % 'Wideband'
            PMISet.i1 = PMI.i1;
            PMISet.i2 = 1;

            PMIInfo.SINRPerRE = SINRsperRE;
            PMIInfo.SINRPerSubband = SINRperSubbandperCW(1,:);
            PMIInfo.W = W;
        end
    else
        % Calculate the SINR and CQI values according to the syntax with
        % the CSI-RS configuration object

        % Get the PMI and SINR values from the PMI selection function
        [PMISet,PMIInfo] = hDLPMISelect(carrier,csirs,reportConfig,nLayers,H,nVar);

        SINRperSubband = NaN(CQISubbandInfo.NumSubbands,nLayers);
        if isfield(reportConfig,'PRGSize') && ~isempty(reportConfig.PRGSize)
            % When PRGSize field is configured as other than empty, the CQI
            % computation is done by choosing one random i2 value from all
            % the i2 values corresponding to the PRGs spanning the subband
            % or the wideband based on the CQI mode, as defined in TS
            % 38.214 Section 5.2.1.4.2
            rng(0); % Set RNG state for repeatability
            randomi2 = zeros(1,CQISubbandInfo.NumSubbands);
            if strcmpi(reportConfig.CQIMode,'Subband')
                % Map the PRGs to subbands
                index = 1;
                thisSubbandSize = CQISubbandInfo.SubbandSizes(1);
                % Get the starting position of each PRG with respect to the
                % current subband. It helps to compute the number of PRGs
                % in the respective subband
                startPRG = ones(1,CQISubbandInfo.NumSubbands+1);
                for prgIdx = 1:numel(PRGInfo.SubbandSizes)
                    if (thisSubbandSize - PRGInfo.SubbandSizes(prgIdx) == 0) && (index < CQISubbandInfo.NumSubbands)
                        % Go to the next subband index and replace the
                        % current subband size
                        index = index + 1;
                        thisSubbandSize = CQISubbandInfo.SubbandSizes(index);
                        % Mark the corresponding PRG index as the start of
                        % subband
                        startPRG(index) = prgIdx + 1;
                    else
                        thisSubbandSize = thisSubbandSize - PRGInfo.SubbandSizes(prgIdx);
                    end
                end
                % Append the total number of PRGs + 1 value to the
                % startPRG vector. The value points to the last PRG at the
                % end of the BWP, to know the number of PRGs in the last
                % subband
                startPRG(index+1) = PRGInfo.NumSubbands+1;
                % Loop over all the subbands and choose an i2 value
                % randomly from the i2 values corresponding to all the PRGs
                % spanning each subband
                for idx = 2:numel(startPRG)
                    i2Set = PMISet.i2(startPRG(idx-1):startPRG(idx)-1);
                    randomi2(idx-1) = i2Set(randi(numel(i2Set)));
                    if ~isnan(randomi2(idx-1))
                        SINRperSubband(idx-1,:) = entropicMean(PMIInfo.SINRPerSubband(startPRG(idx-1):startPRG(idx)-1,:,randomi2(idx-1),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3)));
                    end
                end
            else
                % Choose an i2 value randomly from the i2 values
                % corresponding to all the PRGs in the BWP
                randomi2 = PMISet.i2(randi(size(PMISet.i2,2)));
                SINRperSubband(:,:) = entropicMean(PMIInfo.SINRPerSubband(:,:,randomi2,PMISet.i1(1),PMISet.i1(2),PMISet.i1(3)));
            end
            randomPMISet.i1 = PMISet.i1;
            randomPMISet.i2 = randomi2;
            % Get the SINR values in RB level granularity, based on the
            % random i2 values selected. These values are not directly used
            % for CQI computation. These are just for information purpose
            SINRsperRBperCW = getSINRperRB(PMIInfo.SINRPerRE,randomPMISet,CQISubbandInfo.SubbandSizes);
        else
            % If PRGSize is not configured, the output from PMI selection
            % function is either in wideband or subband level granularity
            % based on the PMIMode

            % Get the SINR values corresponding to the PMISet in RB level
            % granularity. These values are not directly used for CQI
            % computation. These are just for information purpose
            SINRsperRBperCW = getSINRperRB(PMIInfo.SINRPerRE,PMISet,PMISubbandInfo.SubbandSizes);

            % Deduce the SINR values for the CQI computation based on the
            % CQI mode, as the SINRPerSubband field in the PMI information
            % output has the SINR values according to the PMIMode
            if strcmpi(reportConfig.PMIMode,'Wideband')
                % If PMI mode is 'Wideband', only one i2 value is reported
                % and the SINR values are obtained for the entire BWP in
                % the SINRPerSubband field of PMIInfo output. In this case
                % compute the SINR values corresponding to subband or
                % wideband based on the CQI mode. Choose the same i2 value
                % for all subbands
                PMI = PMISet;
                PMI.i2 = PMISet.i2.*ones(CQISubbandInfo.NumSubbands,1);
                SINRperSubband = getSubbandSINR(PMIInfo.SINRPerRE,PMI,CQISubbandInfo);
            else
                % If PMI mode is 'Subband', one i2 value is reported per
                % subband and the SINR values are obtained in subband level
                % granularity from PMI selection function. Extract the SINR
                % values accordingly
                for subbandIdx = 1:numel(PMISet.i2)
                    if ~isnan(PMISet.i2(subbandIdx))
                        SINRperSubband(subbandIdx,:) = PMIInfo.SINRPerSubband(subbandIdx,:,PMISet.i2(subbandIdx),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3));
                    end
                end
            end
        end

        SINRperSubbandperCW = zeros(size(SINRperSubband,1),numCodewords);
        for subbandIdx = 1:size(SINRperSubband,1)
            % Get the SINR values per layer and calculate the SINR values
            % corresponding to each codeword
            layerSINRs = squeeze(SINRperSubband(subbandIdx,:));

            if ~any(isnan(layerSINRs))
                %YXC begin
                % Effective SINR
                %codewordSINRs = cellfun(@sum,nrLayerDemap(layerSINRs));
                f = @(x) log2(1+x);
                finv = @(x) 2.^x - 1;
                codewordSINRs = cellfun(@(x)finv(sum(f(x))),nrLayerDemap(layerSINRs));
                %YXC end
            else
                % If the linear SINR values of the codeword are NaNs, which
                % implies, there are no CSI-RS resources in the current
                % subband. So, the SINR values for the codewords are
                % considered as NaNs for the particular subband
                codewordSINRs = NaN(1,numCodewords);
            end
            SINRperSubbandperCW(subbandIdx,:) = codewordSINRs;
        end

        % Compute the wideband SINR value as a mean of the subband SINRs,
        % if either CQI or PMI are configured in subband mode
        if size(SINRperSubbandperCW,1) > 1
            SINRperSubbandperCW = [entropicMean(SINRperSubbandperCW); SINRperSubbandperCW];
        end
    end

    if all(isnan(PMISet.i1)) && all(isnan(PMISet.i2))
        % If PMISet contains only NaN values, it means that there are no
        % CSI-RS indices present in the slot or the value of nVar is zero
        if CQISubbandInfo.NumSubbands == 1
            % Convert the numSubbands to 0 to report only the wideband CQI
            % index in case of wideband mode
            numSubbands = 0;
        else
            numSubbands = CQISubbandInfo.NumSubbands;
        end
        % Report CQI and the CQI information structure parameters as NaN
        CQI = NaN(numSubbands+1,numCodewords);
        CQIInfo.SINRPerSubbandPerCW = NaN(numSubbands+1,numCodewords);
        CQIInfo.SINRPerRBPerCW = NaN(reportConfig.NSizeBWP,carrier.SymbolsPerSlot,numCodewords);
        CQIInfo.SubbandCQI = NaN(numSubbands+1,numCodewords);
    else
        % Get the CQI value
        CQIForAllSubbands = arrayfun(@(x)getCQI(x,SINRTable),SINRperSubbandperCW);

        % Compute the subband differential CQI value in case of subband
        % mode
        if strcmpi(reportConfig.CQIMode,'Subband')
            % Map the subband CQI values to their subband differential
            % value as defined in TS 38.214 Table 5.2.2.1-1. According to
            % this table, a subband differential CQI value is reported for
            % each subband based on the offset level, where the offset
            % level = subband CQI index - wideband CQI index
            CQIdiff = CQIForAllSubbands(2:end,:) - CQIForAllSubbands(1,:);

            % If the CQI value in any subband is NaN, consider the
            % corresponding subband differential CQI as NaN. It indicates
            % that there are no CSI-RS resources present in that particular
            % subband
            CQIOffset(isnan(CQIdiff)) = NaN;
            CQIOffset(CQIdiff == 0) = 0;
            CQIOffset(CQIdiff == 1) = 1;
            CQIOffset(CQIdiff >= 2) = 2;
            CQIOffset(CQIdiff <= -1) = 3;

            CQIOffset = reshape(CQIOffset,[],numCodewords);
            % Form an output CQI array to include wideband CQI value
            % followed by subband differential values
            CQI = [CQIForAllSubbands(1,:); CQIOffset];
        else
            % In 'Wideband' CQI mode, report only the wideband CQI index
            CQI = CQIForAllSubbands(1,:);
        end

        % Form the output CQI information structure
        CQIInfo.SINRPerRBPerCW = SINRsperRBperCW;
        CQIInfo.SINRPerSubbandPerCW = SINRperSubbandperCW;
        if strcmpi(reportConfig.CQIMode,'Wideband')
            % Output wideband CQI value, if CQIMode is 'Wideband'
            CQIInfo.SubbandCQI = CQIForAllSubbands(1,:);
        else
            % Output wideband CQI value followed by subband CQI values, if
            % CQIMode is 'Subband'
            CQIInfo.SubbandCQI = CQIForAllSubbands;
        end
    end
end
function CQI = getCQI(linearSINR,SINRTable)
%   CQI = getCQI(LINEARSINR,SINRTABLE) returns the maximum CQI value that
%   corresponds to 90 percent throughput by considering these inputs:
%
%   LINEARSINR - The SINR value in linear scale for which the CQI value has
%                to be computed
%   SINRTABLE  - The SINR lookup table using which the CQI value is reverse
%                mapped

    % Convert the SINR values to decibels
    SINRatRxdB  = 10*log10(linearSINR);

    % The measured SINR value is compared with the SINRs in the lookup
    % table. The CQI index corresponding to the maximum SINR value from the
    % table, which is less than the measured value is reported by the UE
    cqiTemp = find(SINRTable(SINRTable <= SINRatRxdB),1,'last');
    if all(isnan(SINRatRxdB))
        CQI = NaN;
    elseif isempty(cqiTemp)
        % If there is no CQI value that corresponds to 90 percent
        % throughput, CQI value is chosen as 0
        CQI = 0;
    else
        CQI = cqiTemp;
    end
end

function SINRsperRBperCW = getSINRperRB(SINRsperRE,PMISet,subbandSizes)
%   SINRSPERRBPERCW = getSINRperRB(SINRSPERRE,PMISET,SUBBANDSIZES) returns
%   the SINR values corresponding to the PMISet in RB level granularity
%   spanning one slot, by considering these inputs:
%
%   SINRSPERRE   - The SINR values per RE for all PMI indices
%   PMISET       - The PMI value reported
%   SUBBANDSIZES - The array representing size of each subband

    numSubbands = numel(PMISet.i2);
    % Get SINR values per RE based on the PMI values
    start = 0;
    SINRsperRECQI = NaN(size(SINRsperRE,1),size(SINRsperRE,2),size(SINRsperRE,3));
    for idx = 1:numSubbands
        if ~isnan(PMISet.i2(idx))
            SINRsperRECQI(start*12 + 1:(start + subbandSizes(idx))*12,:,:) = SINRsperRE(start*12 + 1:(start + subbandSizes(idx))*12,:,:,PMISet.i2(idx),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3));
        end
        start = start + subbandSizes(idx);
    end

    % Calculate SINR value per RE per each codeword
    nLayers = size(SINRsperRECQI,3);
    numCodewords = ceil(nLayers/4);
    SINRsperREperCW = NaN(size(SINRsperRECQI,1),size(SINRsperRECQI,2),numCodewords);
    for k = 1:size(SINRsperRECQI,1)
        for l = 1:size(SINRsperRECQI,2)
            if(~isnan(SINRsperRECQI(k,l,:)))
                SINRsperREperCW(k,l,:) = cellfun(@sum,nrLayerDemap(squeeze(SINRsperRECQI(k,l,:))'));
            end
        end
    end

    % Calculate the SINR value per RB by averaging the SINR values per
    % RE within RB spanning one slot
    SINRsperRBperCW = zeros(size(SINRsperREperCW,1)/12,size(SINRsperREperCW,2),size(SINRsperREperCW,3));
    for RBidx = 0:(size(SINRsperREperCW,1)/12)-1
        % Consider the mean of SINR values over each RB
        RBSINRs = SINRsperREperCW((1:12)+RBidx*12,:,:);
        SINRsperRBperCW(RBidx+1,:,:) = mean(RBSINRs,1,'omitnan');
    end
end

function SubbandSINRs = getSubbandSINR(SINRsperRE,PMISet,SubbandInfo)
%   SUBBANDSINRS = getSubbandSINR(SINRSPERRE,PMISET,SUBBANDINFO) returns
%   the SINR values per subband by averaging the SINR values across all the
%   REs within the subband spanning one slot, corresponding to the reported
%   PMI indices, by considering these inputs:
%
%   SINRSPERRE     - SINR values per RE for all PMI indices
%   PMISET         - The PMI indices according to which SINRs must be
%                    extracted
%   SUBBANDINFO    - Subband information related structure with these 
%   fields:
%      NumSubbands  - Number of subbands
%      SubbandSizes - Size of each subband

    SubbandSINRs = NaN(SubbandInfo.NumSubbands,size(SINRsperRE,3));
    % Consider the starting position of first subband as start of BWP
    subbandStart = 0;
    for SubbandIdx = 1:SubbandInfo.NumSubbands
        if ~isnan(PMISet.i2(SubbandIdx))
            SubbandSINRs(SubbandIdx,:) = squeeze(entropicMean(entropicMean(SINRsperRE((subbandStart*12 + 1):(subbandStart+ SubbandInfo.SubbandSizes(SubbandIdx))*12,:,:,PMISet.i2(SubbandIdx),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3)))));
        end
        subbandStart = subbandStart+ SubbandInfo.SubbandSizes(SubbandIdx);
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

function [reportConfigOut,SINRTable,nVar] = validateInputs(carrier,reportConfig,nLayers,H,nVar,SINRTable,numCSIRSPorts,csirsInd,isCSIRSObjSyntax)
%   [REPORTCONFIGOUT,SINRTABLE,NVAR] = validateInputs(CARRIER,REPORTCONFIG,NLAYERS,H,NVAR,SINRTABLE,NUMCSIRSPORTS,CSIRSIND,ISCSIRSOBJSYNTAX)
%   validates the inputs arguments and returns the validated CSI report
%   configuration structure REPORTCONFIGOUT along with the SINR lookup
%   table SINRTABLE for SINR to CQI mapping, and the noise variance NVAR.

    fcnName = 'hCQISelect';
    % Validate 'reportConfig'
    % Validate 'NSizeBWP'
    if ~isfield(reportConfig,'NSizeBWP')
        error('nr5g:hCQISelect:NSizeBWPMissing','NSizeBWP field is mandatory.');
    end
    nSizeBWP = reportConfig.NSizeBWP;
    if ~(isnumeric(nSizeBWP) && isempty(nSizeBWP))
        validateattributes(nSizeBWP,{'double','single'},{'scalar','integer','positive','<=',275},fcnName,'the size of BWP');
    else
        nSizeBWP = carrier.NSizeGrid;
    end
    % Validate 'NStartBWP'
    if ~isfield(reportConfig,'NStartBWP')
        error('nr5g:hCQISelect:NStartBWPMissing','NStartBWP field is mandatory.');
    end
    nStartBWP = reportConfig.NStartBWP;
    if ~(isnumeric(nStartBWP) && isempty(nStartBWP))
        validateattributes(nStartBWP,{'double','single'},{'scalar','integer','nonnegative','<=',2473},fcnName,'the start of BWP');
    else
        nStartBWP = carrier.NStartGrid;
    end
    if nStartBWP < carrier.NStartGrid
        error('nr5g:hCQISelect:InvalidNStartBWP',...
            ['The starting resource block of BWP ('...
            num2str(nStartBWP) ') must be greater than '...
            'or equal to the starting resource block of carrier ('...
            num2str(carrier.NStartGrid) ').']);
    end
    % BWP must lie within the limits of carrier
    if (nSizeBWP + nStartBWP)>(carrier.NStartGrid + carrier.NSizeGrid)
        error('nr5g:hCQISelect:InvalidBWPLimits',['The sum of starting resource '...
            'block of BWP (' num2str(nStartBWP) ') and the size of BWP ('...
            num2str(nSizeBWP) ') must be less than or equal to '...
            'the sum of starting resource block of carrier ('...
            num2str(carrier.NStartGrid) ') and size of the carrier ('...
            num2str(carrier.NSizeGrid) ').']);
    end
    reportConfigOut.NStartBWP = nStartBWP;
    reportConfigOut.NSizeBWP = nSizeBWP;

    % Check for the presence of panel dimensions parameter and add it to
    % the reportConfig structure, if present. This parameter is used to
    % obtain the precoding matrices and is validated in hDLPMISelect
    % function
    if isfield(reportConfig,'PanelDimensions')
        reportConfigOut.PanelDimensions = reportConfig.PanelDimensions;
    end
    % Check if CQI Mode is specified. Otherwise, by default, consider
    % 'Wideband' mode
    if isfield(reportConfig,'CQIMode')
        validatestring(reportConfig.CQIMode,{'Wideband','Subband'},...
            fcnName,'CQIMode field');
        reportConfigOut.CQIMode = reportConfig.CQIMode;
    else
        reportConfigOut.CQIMode = 'Wideband';
    end
    % Check if PMI Mode is specified. Otherwise, by default, consider
    % 'Wideband' mode
    if isfield(reportConfig,'PMIMode')
        validatestring(reportConfig.PMIMode,{'Wideband','Subband'},...
            fcnName,'PMIMode field');
        reportConfigOut.PMIMode = reportConfig.PMIMode;
    else
        reportConfigOut.PMIMode = 'Wideband';
    end

    % Validate 'SubbandSize'
    NSBPRB = [];
    if strcmpi(reportConfigOut.CQIMode,'Subband') || strcmpi(reportConfigOut.PMIMode,'Subband')
        if nSizeBWP >= 24
            if isCSIRSObjSyntax
                fieldName = 'SubbandSize';
            else
                fieldName = 'NSBPRB';
            end
            if ~isfield(reportConfig,'SubbandSize')
                error('nr5g:hCQISelect:SubbandSizeMissing',...
                    ['For the subband mode, ' fieldName ' field is '...
                    'mandatory when the size of BWP is more than 24 PRBs.']);
            else
                validateattributes(reportConfig.SubbandSize,{'double','single'},...
                    {'real','scalar'},fcnName,fieldName);
                NSBPRB = reportConfig.SubbandSize;
            end
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
            error('nr5g:hCQISelect:InvalidPRGSize',...
                ['PRGSize value (' num2str(reportConfig.PRGSize) ') must be [], 2, or 4.']);
        end
        reportConfigOut.PRGSize = reportConfig.PRGSize;
    else
        reportConfigOut.PRGSize = [];
    end

    if (strcmpi(reportConfigOut.CQIMode,'Subband') ||...
            strcmpi(reportConfigOut.PMIMode,'Subband')) && isempty(reportConfigOut.PRGSize)
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
                error('nr5g:hCQISelect:InvalidSubbandSize',['For the configured BWP size (' num2str(nSizeBWP) ...
                    '), subband size (' num2str(NSBPRB) ') must be ' num2str(validNSBPRBValues(1)) ...
                    ' or ' num2str(validNSBPRBValues(2)) '.']);
            end
        end
    end

    % Consider CodebookMode field, if present, for PMI selection
    if isfield(reportConfig,'CodebookMode')
        reportConfigOut.CodebookMode = reportConfig.CodebookMode;
    end

    % Consider CodebookSubsetRestriction field, if present, for
    % PMI selection
    if isfield(reportConfig,'CodebookSubsetRestriction')
        reportConfigOut.CodebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
    end

    % Consider i2Restriction field, if present, for PMI selection
    if isfield(reportConfig,'i2Restriction')
        reportConfigOut.i2Restriction = reportConfig.i2Restriction;
    end

    % Validate 'nLayers'
    validateattributes(nLayers,{'numeric'},{'scalar','integer','positive','<=',8},fcnName,'NLAYERS');

    % Validate the channel estimate and its dimensions
    validateattributes(numel(size(H)),{'double'},{'>=',2,'<=',4},fcnName,'number of dimensions of H');
    if ~isempty(csirsInd)
        nRxAnts = size(H,3);
        K = carrier.NSizeGrid*12;
        L = carrier.SymbolsPerSlot;
        if ~isCSIRSObjSyntax
            thirdDim = 1;
        else
            thirdDim = NaN;
        end
        validateattributes(H,{class(H)},{'size',[K L thirdDim numCSIRSPorts]},fcnName,'H');
        
        % Validate 'nLayers'
        maxNLayers = min(nRxAnts,numCSIRSPorts);
        if nLayers > maxNLayers
            error('nr5g:hCQISelect:InvalidNumLayers',...
                ['The given antenna configuration (' ...
                num2str(numCSIRSPorts) 'x' num2str(nRxAnts)...
                ') supports only up to (' num2str(maxNLayers) ') layers.']);
        end
    end

    % Validate noise variance
    validateattributes(nVar,{'double','single'},{'scalar','real','nonnegative','finite'},fcnName,'NVAR');
    % Clip nVar to a small noise variance to avoid +/-Inf outputs
    if nVar < 1e-10
        nVar = 1e-10;
    end
    % Validate SINRTable
    if isempty(SINRTable)
        % Default SINRTable is generated by running simulations for 100
        % frames for the CSI reference resource as defined in TS 38.214
        % Section 5.2.2.5, for an AWGN channel for SISO scenario,
        % considering 0.1 BLER condition and TS 38.214 Table 5.2.2.1-2
        SINRTable =  [-5.84  -4.20  -2.08  -0.23  1.66  3.08  5.03  7.02...
                       9.01  10.99  12.99  15.01  16.51  18.49  19.99];
    else
        if isCSIRSObjSyntax
            syntaxString = 'SINRTable';
        else
            syntaxString = 'SINR90pc field';
        end
        validateattributes(SINRTable,{'double','single'},{'vector','real','numel',15},fcnName,syntaxString);
    end
end

function [reportConfig,nLayers,H,nVar,SINRTable,isCSIRSObjSyntax,nTxAnts,csirsInd,out] = parseInputs(carrier,varargin)
%   [REPORTCONFIG,NLAYERS,H,NVAR,SINRTABLE,ISCSIRSOBJSYNTAX,NTXANTS,CSIRSIND,OUT] = parseInputs(CARRIER,CSIRS,REPORTCONFIG,NLAYERS,H,NVAR,SINRTABLE) 
%   returns the parsed arguments and other required parameters for the
%   syntax with CSI-RS configuration object by considering these inputs:
%   CARRIER      - Carrier configuration object
%   CSIRS        - CSI-RS configuration object
%   REPORTCONFIG - Structure of CSI reporting configuration
%   NLAYERS      - Number of transmission layers
%   H            - Estimated channel information 
%   NVAR         - Estimated noise variance
%   SINRTABLE    - SINR lookup table for SINR to CQI mapping
%
%   [REPORTCONFIG,NLAYERS,H,NVAR,SINRTABLE,ISCSIRSOBJSYNTAX,NTXANTS,CSIRSIND,CSIRS] = parseInputs(CARRIER,BWP,CQICONFIG,CSIRSIND,H,NVAR) 
%   returns the parsed arguments and other required parameters for the
%   syntax with CSI-RS indices by considering these inputs:
%   CARRIER      - Carrier configuration object
%   BWP          - Structure of BWP dimensions
%   CQICONFIG    - Structure of CQI reporting configuration
%   CSIRSIND     - CSI-RS indices
%   H            - Estimated channel information
%   NVAR         - Estimated noise variance

    % Validate the carrier configuration object
    validateattributes(carrier,{'nrCarrierConfig'},{'scalar'},'hCQISelect','CARRIER');
    variableInputArgs = varargin{1};
    if isstruct(variableInputArgs{1})
        % If the first variable argument is a structure, the syntax with
        % CSI-RS indices input is considered. This syntax supports CQI
        % computation only for SISO case. Move the required set of
        % parameters that can adapt into the syntax with CSI-RS
        % configuration object and bind them accordingly, in order to
        % enable easy validation

        % Extract bwp from the first variable input argument
        bwp = variableInputArgs{1};

        % Check if the size and start of BWP fields are present in the bwp
        % structure.
        if ~isfield(bwp,'NSizeBWP')
            error('nr5g:hCQISelect:NSizeBWPMissing','NSizeBWP field is mandatory.');
        end
        if ~isfield(bwp,'NStartBWP')
            error('nr5g:hCQISelect:NStartBWPMissing','NStartBWP field is mandatory.');
        end

        % Extract the CQI configuration related parameter, from the second
        % variable input argument
        reportConfig = variableInputArgs{2};
        % Bind the BWP dimensions into the reportConfig structure
        reportConfig.NStartBWP = bwp.NStartBWP;
        reportConfig.NSizeBWP = bwp.NSizeBWP;

        % Extract the CSI-RS indices from the third variable input argument
        csirsInd = variableInputArgs{3};
        % Validate CSI-RS indices
        validateattributes(csirsInd,{'numeric'},{'positive','integer'},'hCQISelect','CSIRSIND');

        % Extract the channel estimation matrix from fourth variable
        % argument
        H = variableInputArgs{4};

        % Extract the noise variance from fifth variable input argument.
        % For this syntax, the noise variance is a mandatory input. So
        % default value is not considered here
        nVar = variableInputArgs{5};

        % Consider the number of transmission layers as 1 for SISO case
        nLayers = 1;

        % Extract the subband size value NSBPRB, if present, and store it
        % as SubbandSize field (as in the syntax with CSI-RS configuration
        % object) in reportConfig structure
        if isfield(reportConfig,'NSBPRB')
            reportConfig.SubbandSize = reportConfig.NSBPRB;
        end

        % Extract the SINR lookup table SINR90pc, if present, and store it
        % as SINRTable
        if isfield(reportConfig,'SINR90pc')
            SINRTable = reportConfig.SINR90pc;
        else
            % If SINR lookup table is not configured, consider SINRTable as
            % empty
            SINRTable = [];
        end

        % Consider the number of transmit antennas as 1 for SISO case
        nTxAnts = 1;
        % Consider a variable to denote if this syntax considers CSI-RS
        % configuration object
        isCSIRSObjSyntax = false;

        % Assign the output argument out as [], in case of the syntax with
        % CSI-RS indices as an input
        out = [];
    elseif isa(variableInputArgs{1},'nrCSIRSConfig')
        % If the first variable input argument is a CSI-RS configuration
        % object, consider the input arguments according to the syntax with
        % CSI-RS configuration object as an input

        % Extract the CSI-RS configuration object as csirs from the first
        % variable input argument
        csirs = variableInputArgs{1};

        % Validate the CSI-RS resources used for CQI computation. All the
        % CSI-RS resources used for the CQI computation must have same CDM
        % lengths and same number of ports according to TS 38.214 Section
        % 5.2.2.3.1
        validateattributes(csirs,{'nrCSIRSConfig'},{'scalar'},'hCQISelect','CSIRS');
        if ~isscalar(unique(csirs.NumCSIRSPorts))
            error('nr5g:hCQISelect:InvalidCSIRSPorts','All the CSI-RS resources must be configured to have same number of CSI-RS ports.');
        else
            % If all the CSI-RS resources have the same number of CSI-RS
            % ports, get the value as number of transmit antennas
            nTxAnts = unique(csirs.NumCSIRSPorts);
        end

        if ~iscell(csirs.CDMType)
            cdmType = {csirs.CDMType};
        else
            cdmType = csirs.CDMType;
        end
        % Check if all the CSI-RS resources are configured to have same
        % CDM lengths
        if ~all(strcmpi(cdmType,cdmType{1}))
            error('nr5g:hCQISelect:InvalidCSIRSCDMTypes','All the CSI-RS resources must be configured to have same CDM lengths.');
        end

        % Ignore zero-power (ZP) CSI-RS resources, as they are not used for
        % CSI estimation
        if ~iscell(csirs.CSIRSType)
            csirs.CSIRSType = {csirs.CSIRSType};
        end
        numZPCSIRSRes = sum(strcmpi(csirs.CSIRSType,'zp'));
        % Calculate the CSI-RS indices
        tempInd = nrCSIRSIndices(carrier,csirs,"IndexStyle","subscript","OutputResourceFormat","cell");
        tempInd = tempInd(numZPCSIRSRes+1:end)';
        csirsInd = cell2mat(tempInd);

        % Extract the CSI reporting related configuration from second
        % variable input argument
        reportConfig = variableInputArgs{2};

        % Extract the number of transmission layers value as nLayers from
        % the third variable input argument
        nLayers = variableInputArgs{3};

        % Extract the channel estimation matrix from the fourth variable
        % input argument
        H = variableInputArgs{4};

        % Get the number of variable input arguments
        numVarInputArgs = length(variableInputArgs);
        % Extract the noise variance and SINR lookup table
        if numVarInputArgs == 4
            nVar = 1e-10;
            SINRTable = [];
        elseif numVarInputArgs == 5
            nVar = variableInputArgs{5};
            SINRTable = [];
        elseif numVarInputArgs == 6
            nVar = variableInputArgs{5};
            SINRTable = variableInputArgs{6};
        end

        % Consider a variable to denote if this syntax considers CSI-RS
        % configuration object
        isCSIRSObjSyntax = true;
        % Assign the output argument out as csirs object
        out = csirs;
    else
        error('nr5g:hCQISelect:InvalidInputsToTheHelper','The second input argument can be either a structure or a CSI-RS configuration object.');
    end
end

function y = entropicMean(x)
    % Entropic mean of SINR's
    f = @(x) log2(1+x);
    finv = @(x) 2.^x - 1;
    y = finv(mean(f(x),1,'omitnan'));
end