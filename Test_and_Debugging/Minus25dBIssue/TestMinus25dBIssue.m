
load simParameters.mat

param = simParameters;
clear simParameters

% siteIdx = 2;
% ueIdx = 2;

numUEsCluster = 2;

f = waitbar(0,'Testing');

for siteIdx = 1:57
    for ueIdx = 1:numUEsCluster

        channel = nrCDLChannel; % CDL channel object
        channel.DelayProfile = 'Custom';

        LinkDir = 1; % UL
        CHParam = ConfigChannel(param, LinkDir, siteIdx, ueIdx);

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

        channel.NumStrongestClusters = 1;
        channel.ClusterDelaySpread = CHParam.c_DS;

        obj.TxAntPanel = param.GNBAntPanelSize;
        obj.RxAntPanel = param.GNBAntPanelSize;

        %UE Antenna Config
        %channel.TransmitAntennaArray.Size = ueTxAntPanel(ueIdx, :);
        channel.TransmitAntennaArray.Size = param.UEAntPanelSize;
        channel.TransmitAntennaArray.ElementSpacing = param.UEAntElementSpacing;
        channel.TransmitAntennaArray.PolarizationAngles = param.UEAntPolarizationAngles;
        %channel.TransmitAntennaArray.Orientation = [0; 0; 0];
        channel.TransmitAntennaArray.Element = param.UEAntElement;
        channel.TransmitAntennaArray.PolarizationModel = param.UEAntPolarizationModel;
        channel.TransmitArrayOrientation = [0; 0; 0];

        %gNB Antenna Config
        channel.ReceiveAntennaArray.Size = obj.RxAntPanel;
        channel.ReceiveAntennaArray.ElementSpacing = param.GNBAntElementSpacing;
        channel.ReceiveAntennaArray.PolarizationAngles = param.GNBAntPolarizationAngles;
        %channel.ReceiveAntennaArray.Orientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
        channel.ReceiveAntennaArray.Element = param.GNBAntElement;
        channel.ReceiveAntennaArray.PolarizationModel = param.GNBAntPolarizationModel;
        channel.ReceiveArrayOrientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];

        waveformInfo = nrOFDMInfo(param.NumRBs, param.SCS);

        channel.SampleRate = waveformInfo.SampleRate;
        chInfo = info(channel);
        
        x = ((siteIdx-1)*numUEsCluster+ueIdx)/(57*numUEsCluster);
        waitbar(x,f);
    end
end

close(f);