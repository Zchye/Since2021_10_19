classdef YusUtilityObj < handle
    %Yu's utility object
    properties
        %Stores the structure CQIInfo per slot per gNB per UE
        CQIInfoSet
        
        %Stores pmiSet returned by hCQISelect inside hNRUEPhy.m
        PMISet
        
        %The triple (slotNum,siteIdx,ueIdx) identifies a CQIInfo
        Triple = cell(3,1);
    end
    methods
        function obj = YusUtilityObj(param, numSlotsSim)
            %Constructor
            d1 = numSlotsSim;
            d2 = param.NumSitesPerCluster;
            %MXC_2
            d3 = param.NumUEsCell;
            obj.CQIInfoSet = cell(d1,d2,d3);  
            obj.PMISet = cell(d1,d2,d3);
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
            %Push the cqiInfo returned by hCQISelect into the CQIInfoSet           
            if ~all(cellfun(@isempty,obj.Triple)==0)
                %Throw an error if there are still emtpy cells in the
                %Triple
                error('Triple is not ready for being stored in CQIInfoSet');
            end
            
            idx = cellfun(@(x) x, obj.Triple); % Cast the cell array Triple to a number array
            obj.CQIInfoSet{idx(1),idx(2),idx(3)} = cqiInfo;
            
%             obj.Triple = cell(3,1);
        end
        
        function storePMISet(obj,pmiSet)
            %Store pmiSet in obj.PMISet
            idx = cellfun(@(x) x, obj.Triple); % Cast the cell array Triple to a number array
            obj.PMISet{idx(1),idx(2),idx(3)} = pmiSet;
        end
    end
end