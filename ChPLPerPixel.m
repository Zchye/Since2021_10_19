classdef ChPLPerPixel < matlab.System
    %Channels and Pathlosses per pixel
    properties
        %Channels is a cell array stores channel objects of the links
        %between a pixel and all gNBs
        Channels
        
        %Pathlosses is a number array stores pathlosses of the links
        %bewteen a pixel and all gNBs
        Pathlosses
        
        Nsites
        
        param
        
        TxAntPanel
        
        RxAntPanel
    end
    methods
        function obj = ChPLPerPixel(param)
            %Constructor
            obj.param = param;
            obj.Nsites = param.NumSitesPerCluster;
            obj.Channels = cell(obj.Nsites,1);
            obj.Pathlosses = zeros(obj.Nsites,1);
            for siteIdx = 1:obj.Nsites
                %For test
                obj.Channels{siteIdx} = nrCDLChannel;
%                 obj.Channels{siteIdx} = @(x) x;
                %For test
                % To use the custom profile, uncomment the following line
                obj.Channels{siteIdx}.DelayProfile = 'Custom';
                %
%                 obj.Channels{siteIdx}.DelayProfile = 'CDL-C';
            end
        end
        
        function updateChPL(obj,PixelPos,varargin)
            % Update the channel objects and pathlosses
            % Optianl input - The site index to be updated; default for
            % updatingall channels.
            
            %For test
            if isempty(varargin)
                for siteIdx = 1:obj.Nsites
%                     configOneChannel(obj,siteIdx,PixelPos);
                    configOneChannel(obj,siteIdx,PixelPos,'PL');
                    %For test
%                     configFSPL(obj,siteIdx,PixelPos);
                end
            else
                mysite = varargin{1};
                %For test
%                 configOneChannel(obj,mysite,PixelPos);
                configOneChannel(obj,mysite,PixelPos,'PL');
%                 configFSPL(obj,mysite,PixelPos);
            end
        end
        
        function configOneChannel(obj,siteIdx,PixelPos,varargin)
            % Configure one channel object
            % optianl input - PL for 38901 pathloss; default for no
            % pathloss.
            CHParam = ConfigChannel(obj.param, 0, siteIdx, PixelPos);
            release(obj.Channels{siteIdx});
            obj.Channels{siteIdx}.TransmitArrayOrientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
            obj.Channels{siteIdx}.CarrierFrequency = obj.param.DLCarrierFreq;
            obj.Channels{siteIdx}.HasLOSCluster = logical(CHParam.LOS);
            obj.Channels{siteIdx}.AngleSpreads = [CHParam.LSPs.ASD, CHParam.LSPs.ASA, CHParam.LSPs.ZSD, CHParam.LSPs.ZSA];
            if CHParam.LOS
                obj.Channels{siteIdx}.KFactorFirstCluster = CHParam.LSPs.K;
            end 
            obj.Channels{siteIdx}.PathDelays = CHParam.tau_n';
            obj.Channels{siteIdx}.AveragePathGains = CHParam.P_n';
            obj.Channels{siteIdx}.AnglesAoA = CHParam.ADAngles.phi_n_AOA';
            obj.Channels{siteIdx}.AnglesAoD = CHParam.ADAngles.phi_n_AOD';
            obj.Channels{siteIdx}.AnglesZoA = CHParam.ADAngles.theta_n_ZOA';
            obj.Channels{siteIdx}.AnglesZoD = CHParam.ADAngles.theta_n_ZOD';
            obj.Channels{siteIdx}.XPR = CHParam.XPR;
            obj.Channels{siteIdx}.NumStrongestClusters = 2;
            obj.Channels{siteIdx}.ClusterDelaySpread = CHParam.c_DS;
            obj.TxAntPanel = obj.param.UETxAntPanelSize;
            obj.RxAntPanel = obj.param.UERxAntPanelSize;
            obj.Channels{siteIdx}.TransmitAntennaArray.Size = obj.param.GNBTxAntPanelSize;
            obj.Channels{siteIdx}.TransmitAntennaArray.ElementSpacing = obj.param.GNBTxAntElementSpacing;
            obj.Channels{siteIdx}.TransmitAntennaArray.PolarizationAngles = obj.param.GNBTxAntPolarizationAngles;
            obj.Channels{siteIdx}.TransmitAntennaArray.Element = obj.param.GNBAntElement;
            obj.Channels{siteIdx}.TransmitAntennaArray.PolarizationModel = obj.param.GNBAntPolarizationModel;
            obj.Channels{siteIdx}.TransmitArrayOrientation = [CHParam.bearing; CHParam.downtilt; CHParam.slant];
            obj.Channels{siteIdx}.ReceiveAntennaArray.Size = obj.RxAntPanel;
            obj.Channels{siteIdx}.ReceiveAntennaArray.ElementSpacing = obj.param.UERxAntElementSpacing;
            obj.Channels{siteIdx}.ReceiveAntennaArray.PolarizationAngles = obj.param.UERxAntPolarizationAngles;
            obj.Channels{siteIdx}.ReceiveAntennaArray.Element = obj.param.UEAntElement;
            obj.Channels{siteIdx}.ReceiveAntennaArray.PolarizationModel = obj.param.UEAntPolarizationModel;
            obj.Channels{siteIdx}.ReceiveArrayOrientation = [0; 0; 0];
            waveformInfo = nrOFDMInfo(obj.param.NumRBs, obj.param.SCS);
            obj.Channels{siteIdx}.SampleRate = waveformInfo.SampleRate;            
            
            if ~isempty(varargin)
                if isequal(varargin{1},'PL')
                    obj.Pathlosses(siteIdx) = CHParam.Pathloss;
                else
                    error('To use 38901 pathloss, please specify the optional input to be "PL".');
                end
            end
                
        end
        
        function configFSPL(obj,siteIdx,PixelPos)
            %Configure free space pathloss
            R = norm(obj.param.GNBPositions(siteIdx,:)-PixelPos);
            lambda = physconst('LightSpeed')/obj.param.DLCarrierFreq;
            obj.Pathlosses(siteIdx) = fspl(R,lambda);
        end
        
        function y = applyChannel(obj,idx,x)
            %Apply a channel to an input signal
            %For test
            %If the channe is of nrCDLChannel, append zeros to the signal
            %so that the delays can be simulated
            chInfo = info(obj.Channels{idx});
            MaxChannelDelay = ceil(max(chInfo.PathDelays*obj.Channels{idx}.SampleRate)) + chInfo.ChannelFilterDelay;
            x = [x;zeros(MaxChannelDelay,size(x,2))];
            y = obj.Channels{idx}(x);
        end
        
        function y = applyPL(obj,idx,x)
            %Apply pathloss to a signal
            s = 10^(-obj.Pathlosses(idx)/20);
            y = s*x;
        end
    end
    
    methods (Access = protected)
        function y = stepImpl(obj,idx,x,varargin)
            %The step function
            %Optional input paramter: 0 - Channel only; 1 - Pathloss only;
            %default - apply both channel and pathloss
            mark = -1; % -1 corresponds to no optional inputs
            if ~isempty(varargin)
                mark = varargin{1};
            end
            
            switch mark
                case 0
                    y = applyChannel(obj,idx,x);
                case 1
                    y = applyPL(obj,idx,x);
                otherwise 
                    y = applyChannel(obj,idx,x);
                    y = applyPL(obj,idx,y);
            end
        end
    end
end