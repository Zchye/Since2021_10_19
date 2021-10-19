classdef (Abstract) hWirelessNode < handle
%hWirelessNode Base class for nodes of wireless networks
%
%   hWirelessNode properties:
%
%   NodeID                  - Node identifier
%   NodePosition            - Node position
%   DistanceCalculatorFcn   - Reference to distance calculator function
%   Mobility                - Reference to mobility model object

%   Copyright 2020 The MathWorks, Inc.

    properties
        %NodeID Node identifier
        %   Specify this property as an integer. This is the unique
        %   identifier for this node in the network.
        NodeID = 1

        %NodePosition Node position
        %   Specify this property as a row vector with 3 elements. This
        %   property identifies the position of the node in the network.
        NodePosition = [0 0 0]

        %DistanceCalculatorFcn Function handle for distance calculations
        %   Specify this property as a function handle to the distance
        %   calculator function. If empty, Euclidean distance is used for
        %   node distance computation.
        DistanceCalculatorFcn

        %Mobility Node mobility model
        %   Specify this property as a mobility object. This object is invoked
        %   to determine latest position, whenever position is queried
        Mobility
    end
 
    methods
        % Constructor
        function obj = hWirelessNode(varargin)
            % Name-value pairs
            for idx = 1:2:nargin
                obj.(varargin{idx}) = varargin{idx+1};
            end
        end

        % Set Node ID
        function set.NodeID(obj, value)
            validateattributes(value, {'numeric'}, {'scalar', ...
                'positive', 'integer'}, mfilename, 'NodeID');
            obj.NodeID = value;
        end

        % Set node position
        function set.NodePosition(obj, value)
            validateattributes(value, {'numeric'}, {'row', 'numel', 3}, ...
                mfilename, 'NodePosition');
            obj.NodePosition = value;
        end

        % Get node position
        function value = get.NodePosition(obj)

            if ~isempty(obj.Mobility)
                value = Position(obj.Mobility, obj.getCurrentTime);
            else
                value = obj.NodePosition;
            end
        end

        % Get distance of node from transmitter node
        function value = getNodeDistance(obj, txPosition)

            validateattributes(txPosition, {'numeric'}, {'row', 'numel', 3});
            if ~isempty(obj.DistanceCalculatorFcn)
                value = obj.DistanceCalculatorFcn(txPosition, obj.NodePosition);
            else
                % Use Euclidean distance for distance calculation if
                % function handle is not specified
                value = norm(txPosition - obj.NodePosition);
            end
        end

        % Register function handle for calculating node distances
        function set.DistanceCalculatorFcn(obj, distanceCalculatorFcn)

            if isa(distanceCalculatorFcn, 'function_handle')
                obj.DistanceCalculatorFcn = distanceCalculatorFcn;
            else
                error('hWirelessNode:InvalidFunctionHandle', 'Input argument must be a proper function handle');
            end
        end

        % Set the mobility model
        function set.Mobility(obj, mobilityModelObj)
            % Set the mobility model
            %
            % MOBILITYMODELOBJ - Handle of the mobility object

            if isvalid(mobilityModelObj)
                obj.Mobility = mobilityModelObj;
            else
                error('hWirelessNode:InvalidObjectHandle', 'Input argument must be a valid object handle');
            end
        end
    end

    methods(Abstract)
        currentime = getCurrentTime(obj);
    end
end