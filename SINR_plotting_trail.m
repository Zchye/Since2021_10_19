function SINR_plotting_trail(YUO)
%MXC_2
%load('YUO');
%within the object, CQIInfoSet is a cell array that stores all downlink CQI information
%its dimension is number of slots in the simulation X number of sites in the simulation X total number of UEs

%get dimensions
[numSlots, numSites, numUEs] = size(YUO.CQIInfoSet, 1, 2, 3);

CQI_info_set_raw = YUO.CQIInfoSet; %numSlots X numSites X numUEs cell array

%YXC begin
DMRSSINR_raw = YUO.DMRSSINR;
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


UE_Wideband_SINR=zeros(numSites,numUEs,numSlots);
%UE_Wideband_SINR is an array with size equal number of sites X number of UEs X number of slots
%use this variable to store the wideband SINR of each UE during each slot

%YXC begin
DMRSSINR = zeros(numSites,numUEs,numSlots);
%YXC end

%findout how many samples was collected for each UE
%count how many times this UE was updated
UE_Wideband_SINR_counter=zeros(numSites,numUEs);

%YXC begin
DMRSSINR_counter = zeros(numSites,numUEs);

CQIOld = YUO.CQIOld;
OldLinSINR = zeros(numSites,numUEs,numSlots);
OLSCounter = zeros(numSites,numUEs);
%YXC end


for i = 1:numSites
    for j = 1:numUEs
        for k = 1:numSlots
            if ~isempty(CQI_info_set_raw{k,i,j})
                Wideband_SINR_temp_1 = CQI_info_set_raw(k,i,j); %1 X 1 cell
                Wideband_SINR_temp_2 = Wideband_SINR_temp_1{1}; %1 X 1 struct
                Wideband_SINR_temp_3 = Wideband_SINR_temp_2.SINRPerSubbandPerCW; %(numSubbands + 1) X numCodewords double
                UE_Wideband_SINR(i,j,k) = Wideband_SINR_temp_3(1); % 1 X 1 double
            
                UE_Wideband_SINR_counter(i,j) = UE_Wideband_SINR_counter(i,j) + 1;
            end
            %YXC begin
            if ~isempty(DMRSSINR_raw{k,i,j})
                DMRSSINR(i,j,k) = DMRSSINR_raw{k,i,j};
                DMRSSINR_counter(i,j) = DMRSSINR_counter(i,j)+1;
            end
            
            if ~isempty(CQIOld{k,i,j})
                OldLinSINR(i,j,k) = CQIOld{k,i,j}.SINRPerSubbandPerCW(1);
                OLSCounter(i,j) = OLSCounter(i,j) + 1;
            end
            %YXC end
        end
    end
end


% Calculate entropic average
f = @(x) log2(1+x);
finv = @(x) 2.^x - 1;
UE_average_SINR = finv((sum(f(UE_Wideband_SINR), 3))./UE_Wideband_SINR_counter);

UE_average_SINR_dB = 10*log10(UE_average_SINR);
% 10 not 20 because SINR is the ratio of powers and powers should be
% applied with coefficient 10 when converted from linear to dB.

UE_average_SINR_dB = reshape(UE_average_SINR_dB,[],1);

[f,x]=ecdf(UE_average_SINR_dB);

figure('name','CDF UEs DL SINR [dB]');
plot(x, f);
title('SINR CDF','FontSize',12);
xlabel('UE Average SINR [dB]','FontSize',12);
ylabel('C.D.F','FontSize',12);
if any(DMRSSINR_counter)
    UEAverageDMRSSINR = (sum(DMRSSINR,3))./DMRSSINR_counter;
    UEAverageDMRSSINR_dB = 10*log10(UEAverageDMRSSINR);
    [f_dmrs,x_dmrs] = ecdf(UEAverageDMRSSINR_dB(:));
    hold on
    plot(x_dmrs,f_dmrs)
    legend('CSI-RS','DMRS')
end
if any(OLSCounter)
    LinAveSINR = sum(OldLinSINR,3)./OLSCounter;
    LASdB = 10*log10(LinAveSINR);
    hold on
    [f_old, x_old] = ecdf(LASdB(:));
    plot(x_old, f_old)
    legend('Eff CSI-RS', 'Old CSI-RS')
end

end