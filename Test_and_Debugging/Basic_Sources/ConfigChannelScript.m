tStart = tic;

simParameters.Scenario='UMa';
simParameters.NumClusters=1;
simParameters.NumUEsCell=25; %Vary for testing
simParameters.DuplexMode=1;%0 - FDD; 1 - TDD
simParameters.InterSiteDistance=200;%meters
simParameters.gNBHeight=25;%meters
simParameters.UEHeight=1.5;%meters
simParameters.AntennaDowntilt=6;%degrees
simParameters.AntennaSlant=0;%degrees
simParameters.DLCarrierFreq=30e9;%system DL centre frequency in Hz
load FastFadingTabs.mat;%Load 3GPP TR 38.901 Table 7.5-6, 7.5-7, 7.5-8, 7.5-9, 7.5-10, 7.5-11 to workspace
simParameters.FastFadingTabs=FastFadingTabs;%Functions cannot use the variables from workspace directly

%Generate 3D locations of gNB and UE, as well as gNB bearing angles, in the
%global coordinate system, as per 3GPP TR 38.901 section 7.5 step 1 e) and
%c)
[gNBBearing, ...
    gNBCoordinates, ...
    ueCoordinates, ...
    cellCoordinates] = hMacrocellTopology(simParameters);

simParameters.gNBBearing = gNBBearing;
simParameters.gNBCoordinates = gNBCoordinates;
simParameters.ueCoordinates = ueCoordinates;

[DLTab,~] = ConfigChannel(simParameters);

tEnd = toc(tStart);
disp(['Time elasped: ',num2str(tEnd),'s']);