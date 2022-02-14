function SINR_plotting_trail(YUO)
%MXC_2
%load('YUO');
%within the object, CQIInfoSet is a cell array that stores all downlink CQI information
%its dimension is number of slots in the simulation X number of sites in the simulation X total number of UEs

%get dimensions
% [numSlots, numSites, numUEs] = size(YUO.CQIInfoSet, 1, 2, 3);

% CQI_info_set_raw = YUO.CQIInfoSet; %numSlots X numSites X numUEs cell array

%YXC begin
% DMRSSINR_raw = YUO.DMRSSINR;
%YXC end

%
%YusUtilityObj, last dimension too big, param.NumUEsCell is sufficient
%

%each structure is a CQIInfo structure from hCQISelect.m

%   CQIINFO is an output structure for the CQI information 
%   within its fields:
%   SINRPerSubbandPerCW - It represents the linear SINR values in each
%                         subband for all the codewords. It is a
%                         two-dimensional matrix of size
%                            - 1-by-numCodewords, when CQI reporting mode
%                              is 'Wideband'
%                            - (numSubbands + 1)-by-numCodewords, when
%                              CQI reporting mode is 'Subband'
%                         Each column of the matrix contains wideband SINR
%                         value (the average SINR value across all
%                         subbands) followed by the SINR values of each
%                         subband. The SINR value in each subband is taken
%                         as an average of SINR values of all the REs
%                         across the particular subband spanning one slot


% UE_Wideband_SINR=zeros(numSites,numUEs,numSlots);
%UE_Wideband_SINR is an array with size equal number of sites X number of UEs X number of slots
%use this variable to store the wideband SINR of each UE during each slot

%YXC begin
% DMRSSINR = zeros(numSites,numUEs,numSlots);
%YXC end

%findout how many samples was collected for each UE
%count how many times this UE was updated
% UE_Wideband_SINR_counter=zeros(numSites,numUEs);

%YXC begin
% DMRSSINR_counter = zeros(numSites,numUEs);
%YXC end


% for i = 1:numSites
%     for j = 1:numUEs
%         for k = 1:numSlots
%             if ~isempty(CQI_info_set_raw{k,i,j})
%                 Wideband_SINR_temp_1 = CQI_info_set_raw(k,i,j); %1 X 1 cell
%                 Wideband_SINR_temp_2 = Wideband_SINR_temp_1{1}; %1 X 1 struct
%                 Wideband_SINR_temp_3 = Wideband_SINR_temp_2.SINRPerSubbandPerCW; %(numSubbands + 1) X numCodewords double
%                 UE_Wideband_SINR(i,j,k) = Wideband_SINR_temp_3(1); % 1 X 1 double
%             
%                 UE_Wideband_SINR_counter(i,j) = UE_Wideband_SINR_counter(i,j) + 1;
%             end
%             %YXC begin
%             if ~isempty(DMRSSINR_raw{k,i,j})
%                 DMRSSINR(i,j,k) = DMRSSINR_raw{k,i,j};
%                 DMRSSINR_counter(i,j) = DMRSSINR_counter(i,j)+1;
%             end
%             %YXC end
%         end
%     end
% end

% UE_average_SINR = (sum(UE_Wideband_SINR, 3))./UE_Wideband_SINR_counter;

%YXC begin
% UEAverageDMRSSINR = (sum(DMRSSINR,3))./DMRSSINR_counter;
%YXC end

% UE_average_SINR_dB = 20*log10(UE_average_SINR);
% 10 or 20?
%YXC begin
% UEAverageDMRSSINR_dB = 20*log10(UEAverageDMRSSINR);
%YXC end

%UE_average_SINR_dB = 10*log10(UE_average_SINR);
% 10 not 20 because SINR is the ratio of powers and powers should be
% applied with coefficient 10 when converted from linear to dB.

% UE_average_SINR_dB = reshape(UE_average_SINR_dB,[],1);


% [f,x]=ecdf(UE_average_SINR_dB);

%YXC begin
% [f_dmrs,x_dmrs] = ecdf(UEAverageDMRSSINR_dB(:));
%YXC end

% figure('name','CDF UEs DL SINR [dB]');
% plot(x, f);
% title('SINR CDF','FontSize',12);
% xlabel('UE Average SINR [dB]','FontSize',12);
% ylabel('C.D.F','FontSize',12);
% hold on
% plot(x_dmrs,f_dmrs)
% legend('CSI-RS','DMRS')

%YXC begin
% figure('name','CDF UEs DL DMRS SINR [dB]')
% plot(x_dmrs,f_dmrs)
% title('DMRS SINR CDF','FontSize',12)
% xlabel('UE Average DMRS SINR [dB]','FontSize',12)
% ylabel('C.D.F','FontSize',12)
%YXC end

    if all(isnan(YUO.GodSINR(:))) % Check if there is valid data at all
        disp('No valid data in YUO.GodSINR')
        return
    end
    % Process God SINR
    GodSINR = YUO.GodSINR; % Linear SINR numeric array, NaN's are filled in missing values
    GodSINR = permute(GodSINR,[2,3,1]); % Rearrage to NumSites-by-NumUEs-by-NumSlots
    entro = @(x) log2(1+x); % Function handle for computing entropies
    entroinv = @(x) 2.^x-1; % Function handle for computing the inverses of entropies
    lin2db = @(x) 10*log10(x); % Convert sinr from linear to logarithm
    CalcCDF = @(sinr) ecdf(reshape(lin2db(sinr),[],1)); % Function handle for getting CDF of SINR
    % Construct 3 types of SINRs
    sinr = cell(1,3);
    sinr{1} = GodSINR; % No averaging
    sinr{2} = mean(GodSINR,3,'omitnan'); % Linearly averaging
    sinr{3} = entroinv(mean(entro(GodSINR),3,'omitnan')); % Entropically averaging
    
    % Plot God SINRs
    figure
    title('God SINRs vs 3GPP SINRs')
    hold on
    for ii = 1:length(sinr)
        [f,x] = CalcCDF(sinr{ii});
        plot(x,f,'LineWidth',2)
    end
    
    % Plot 3GPP SINRs
    load RMaA700MBDLGeom.mat RMaA700MBDLGeom
    NCat = size(RMaA700MBDLGeom,2); % Number of categories
    for ii = 1:NCat
        [f,x] = ecdf(RMaA700MBDLGeom{:,ii});
        plot(x,f)
    end
    
    % Configure ticks in x-axis and show grid
    xt = xticks;
    xtMx = ceil(max(xt));
    xtMn = floor(min(xt));
    xtStep = floor((xtMx-xtMn)/50); % Place 51 ticks on x-axis
    NewTicks = xtMn:xtStep:xtMx;
    xticks(NewTicks)
    grid on
    
    % Insert legend
    lgd = {'No Averaging', 'Linearly Averaging', 'Entropically Averaging'};
    EmptyChar = @(x) ''; % Suppress 3GPP's part of the legend
    lgd = [lgd, cellfun(EmptyChar, cell(1,NCat),'UniformOutput',false)];
    legend(lgd)
end