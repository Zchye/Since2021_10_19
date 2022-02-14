classdef YusUtilityObj < handle
    %Yu's utility object
    properties
        %Stores the structure CQIInfo per slot per gNB per UE
        CQIInfoSet
        
        % Stores the SINR's computed from DMRS in linear scale
        DMRSSINR
        
        GodSINR
        
        %The triple (slotNum,siteIdx,ueIdx) identifies a CQIInfo
        Triple = cell(3,1);
        
        % Stores throughputs
        Throughput
        
        % Metrics step size
        MetricsStepSize
        
        % Indices of cells of interest
        cellOfInterestIdx
    end
    methods
        function obj = YusUtilityObj(param, numSlotsSim)
            %Constructor
            d1 = numSlotsSim;
            d2 = param.NumSitesPerCluster;
            %MXC_2
            d3 = param.NumUEsCell;
            obj.CQIInfoSet = cell(d1,d2,d3);
            obj.DMRSSINR = cell(d1,d2,d3);
            % obj. Throughput is a 4-dimensional array. The first index
            % stores link direction, 1 for DL, 2 for UL. The second index
            % stores UE indices and extra two item for cell throughput and
            % theoretical peak throughput. The third index stores site
            % indices. The fourth index stores slot numbers.
            obj.MetricsStepSize = param.MetricsStepSize;
            d4 = param.NumMetricsSteps;
            obj.Throughput = zeros(2, d3+2, d2, d4);
            obj.GodSINR = nan(d1,d2,d3);
        end
        
        function pushSlotNum(obj,slotNum)
            %Push the slot number into the Triple
            obj.Triple{1} = slotNum;
        end
        
        function pushSiteIdx(obj,siteIdx)
            %Push the site index into the Triple
            obj.Triple{2} = siteIdx;
        end
        
        function pushUEIdx(obj,ueIdx)
            %Push the UE index into the Triple            
            obj.Triple{3} = ueIdx;
        end
        
        function pushTriple(obj, slotNum, siteIdx, ueIdx)
            %Push the indices into the Triple
            pushSlotNum(obj,slotNum);
            pushSiteIdx(obj,siteIdx);
            pushUEIdx(obj,ueIdx);
        end
        
        function storeCQIInfo(obj,cqiInfo)
            %Push the cqiInfo returned by hCQISelect into the Triple           
%             if ~all(cellfun(@isempty,obj.Triple)==0)
%                 %Throw an error if there are still emtpy cells in the
%                 %Triple
%                 error('Triple is not ready for being stored in CQIInfoSet');
%             end
            
            idx = cellfun(@(x) x, obj.Triple); % Cast the cell array Triple to a number array
            obj.CQIInfoSet{idx(1),idx(2),idx(3)} = cqiInfo;
            
%             obj.Triple = cell(3,1);
        end
        
        function storeDMRSSINR(obj, SINR)
            % Stores the SINR's computed from DMRS in linear scale in
            % DMRSSINR
            
            idx = cellfun(@(x) x, obj.Triple); % Cast the cell array Triple to a number array
            obj.DMRSSINR{idx(1),idx(2),idx(3)} = SINR;
        end
        
        function storeThroughput(obj,throughputServed, siteIdx, slotNum)
            % Stores throughput in obj.Throughput
            d4Idx = slotNum/obj.MetricsStepSize;
            obj.Throughput(:,:,siteIdx,d4Idx) = throughputServed;
        end
        
        function storeGodSINR(obj, SINR)
            idx = cellfun(@(x) x, obj.Triple); % Cast the cell array Triple to a number array
            obj.GodSINR(idx(1),idx(2),idx(3)) = SINR;
        end
        
        function SaveFile(obj)
            % Save the simulation data
            YUO.CQIInfoSet = obj.CQIInfoSet;
            YUO.DMRSSINR = obj.DMRSSINR;
            YUO.Throughput = obj.Throughput;
            YUO.GodSINR = obj.GodSINR;
            save('outputYUO.mat','YUO')
        end
    end
end