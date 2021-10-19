classdef hNRNodeDistanceCalculator < handle
%hNRNodeDistanceCalculator Computes the distance between nodes in the network
%   This class computes the distance between the transmitter and receiver
%   nodes for use in toroidal wrap-around modeling scenarios. With wrap-around
%   enabled, the node separation is taken as the minimum of the Euclidean
%   distance and six other wrap-around distances.

%   Copyright 2020 The MathWorks, Inc.
    properties
        % InterSiteDistance Distance between 2 adjacent gNBs in meters
        InterSiteDistance(1, 1){mustBeGreaterThan(InterSiteDistance, 0), mustBeFinite} = 1732;

        % EnableWrapAround Flag to enable/disable wrap-around distance calculations
        EnableWrapAround(1, 1){mustBeNumericOrLogical} = true;
    end

    properties(Access = private)
        % NumSitesPerCluster Number of gNBs in the cluster
        NumSitesPerCluster = 19;

        % NumCellsPerSite Number of cells connected to each gNB
        NumCellsPerSite = 1;        
    end

    methods        
        function obj = hNRNodeDistanceCalculator(param)
            %hNRNodeDistanceCalculator Construct an instance of the class
            % hNRNodeDistanceCalculator(PARAM) creates an object with specific
            % cluster configurations as per the fields of PARAM structure.
            % PARAM is a structure with the following fields
            %   InterSiteDistance  - Distance between adjacent gNBs
            %   EnableWrapAround   - Enable/disable cell wrap-around

            if isfield(param, 'InterSiteDistance')
                obj.InterSiteDistance = param.InterSiteDistance;
            end
            if isfield(param, 'EnableWrapAround')
                obj.EnableWrapAround = param.EnableWrapAround;
            end
        end

        function distance = getDistance(obj, txPosition, rxPosition)
            %getDistance Calculate the distance between the nodes
            %   DISTANCE = GETDISTANCE(OBJ, TXPOSITION, RXPOSITION) returns
            %   the distance between the transmitter and receiver nodes.
            %   TXPOSITION is a 1-by-3 vector containing the transmiiter
            %   position in cartesian coordinates. RXPOSITION is a 1-by-3
            %   vector containing the receiver position in cartesian
            %   coordinates.
            %
            %   DISTANCE contains the calculated distance.

            distance = norm(rxPosition - txPosition);
            if obj.EnableWrapAround
                % Calculate the six wrap-around distances
                wrapAroundDist = zeros(6,1);
                wrapAroundDist(1) = norm(rxPosition - (txPosition - [1.732*obj.InterSiteDistance, 4*obj.InterSiteDistance, 0]));
                wrapAroundDist(2) = norm(rxPosition - (txPosition - [-1.732*obj.InterSiteDistance, -4*obj.InterSiteDistance, 0]));
                wrapAroundDist(3) = norm(rxPosition - (txPosition - [2.598*obj.InterSiteDistance, -3.5*obj.InterSiteDistance, 0]));
                wrapAroundDist(4) = norm(rxPosition - (txPosition - [-2.598*obj.InterSiteDistance, 3.5*obj.InterSiteDistance, 0]));
                wrapAroundDist(5) = norm(rxPosition - (txPosition - [4.33*obj.InterSiteDistance, 0.5*obj.InterSiteDistance, 0]));
                wrapAroundDist(6) = norm(rxPosition - (txPosition - [-4.33*obj.InterSiteDistance, -0.5*obj.InterSiteDistance, 0]));
                distance = min([distance; wrapAroundDist]);
            end
        end
    end
end