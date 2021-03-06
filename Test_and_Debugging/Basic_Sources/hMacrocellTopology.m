function [gNBBearing, gNBCoordinates, ueCoordinates, cellCoordinates] = hMacrocellTopology(param)
%hMacrocellTopology Calculate gNB and UE positions for the network
%
%   [GNBCOORDINATES, UECOORDINATES] = hMacrocellTopology(PARAM) returns the gNB
%   and UE positions for the given 19-cell topology. The gNB is assumed to
%   be at the center of each cell. PARAM is a structure with the following
%   fields:
%   NumClusters       - Number of clusters. Specify NumClusters
%                       value as a positive integer in the range [1, 7]
%   InterSiteDistance - Distance between two adjacent gNBs. Specify
%                       InterSiteDistance as a numeric positive scalar.
%   NumUEsCell        - Number of UEs per cell. Specify NumUEsCell as a 
%                       numeric positive scalar.
%YXC begin
%   gNBHeight         - The antenna height shared by all gNBs.
%   UEHeight          - The antenna height shared by all UEs.
%YXC end
%
%   GNBCOORDINATES is an matrix of size NumClusters*19-by-3. The 'i'th row
%   contains the gNB coordinates of the cell with cluster number given by
%   floor(i/19) and cell ID given by modulus(i-1, 19).
%   UECOORDINATES is a cell array of length NumClusters*19. Each element at
%   the index 'i' is a matrix of size NumUEsCell-by-3 with each row
%   representing the position of UEs in the cell with cell ID 'i-1'.

%   Copyright 2021 The MathWorks, Inc.

% Validate the number of clusters in the network
validateattributes(param.NumClusters, {'numeric'}, {'integer', 'scalar', '>=', 1, '<=',7}, 'param.NumClusters', 'NumClusters')
% Validate the intersite distance values
validateattributes(param.InterSiteDistance, {'numeric'}, {'real', 'finite', 'scalar', '>', 0}, 'param.InterSiteDistance', 'InterSiteDistance')
% Validate the number of UEs in each cell
validateattributes(param.NumUEsCell, {'numeric'}, {'nonempty', 'integer', 'scalar','>', 0, '<=', 65519}, 'param.NumUEsCell', 'NumUEsCell')

%MXC
%new topology, cell radius is one third of ISD
%cellRadius = param.InterSiteDistance/sqrt(3);
cellRadius = param.InterSiteDistance/3;
%MXC

% x-coordinate and y-coordinate of cluster centres normalized with respect to the cell radius.
xClusterCentre = [0, 9/2, 15/2, 3, -9/2, -15/2, -3];
yClusterCentre = [0, 7*sqrt(3)/2, -sqrt(3)/2, -4*sqrt(3), -7*sqrt(3)/2, sqrt(3)/2, 4*sqrt(3)];

%MXC 
%new cell and gNB locations
%{
% xGNBPos and yGNBPos are x and y coordinates of base stations in a cluster respectively
xGNBPos = [0, sqrt(3)*cosd(30:60:360), 3*cosd(0), 2*sqrt(3)*cosd(30), 3*cosd(60), 2*sqrt(3)*cosd(90), 3*cosd(120), 2*sqrt(3)*cosd(150), 3*cosd(180), 2*sqrt(3)*cosd(210), 3*cosd(240), 2*sqrt(3)*cosd(270), 3*cosd(300), 2*sqrt(3)*cosd(330)];
yGNBPos = [0, sqrt(3)*sind(30:60:360), 3*sind(0), 2*sqrt(3)*sind(30), 3*sind(60), 2*sqrt(3)*sind(90), 3*sind(120), 2*sqrt(3)*sind(150), 3*sind(180), 2*sqrt(3)*sind(210), 3*sind(240), 2*sqrt(3)*sind(270), 3*sind(300), 2*sqrt(3)*sind(330)];

% Minimum offset required in x and y directions to ensure positive
% coordinates for all the vertices
if param.NumClusters == 1
    x = 4;
    y = 5*sqrt(3)/2;
else
    x = 11.5;
    y = 13*sqrt(3)/2;
end

numCellsPerCluster = 19;
%}

%position of 19 sites assuming hexagon with side length 1
xSitePos = [0, 3*cosd(0:60:300), sqrt(27)*cosd(30:60:330), 6*cosd(0:60:300)];
ySitePos = [0, 3*sind(0:60:300), sqrt(27)*sind(30:60:330), 6*sind(0:60:300)];

%gNB locations, all three gNBs within a site are at the same location
xGNBPos = repelem(xSitePos,3);
yGNBPos = repelem(ySitePos,3);

%since gNBs are on the edge of cells, need cell center for UE generation
xCellPos = xGNBPos;
yCellPos = yGNBPos;

xCellPos(1:3:end) = xCellPos(1:3:end)+cosd(0);
xCellPos(2:3:end) = xCellPos(2:3:end)+cosd(120);
xCellPos(3:3:end) = xCellPos(3:3:end)+cosd(240);

yCellPos(1:3:end) = yCellPos(1:3:end)+sind(0);
yCellPos(2:3:end) = yCellPos(2:3:end)+sind(120);
yCellPos(3:3:end) = yCellPos(3:3:end)+sind(240);

% Minimum offset required in x and y directions to ensure positive
% coordinates for all the vertices
if param.NumClusters == 1
    x = 8;
    y = 4*sqrt(3);
else
    error('did not implement wrap-around yet')
end

numCellsPerCluster = 57;
%MXC

%YXC begin
%Assign antenna bearing angles to each gNB.
temp=[0;120;240];
temp=repelem(temp,1,numCellsPerCluster/3);
gNBBearing=temp(:);
%YXC end

% gNB positions
gNBCoordinates = zeros(numCellsPerCluster, 3);
%MXC
cellCoordinates = zeros(numCellsPerCluster, 2);
%MXC
for clstIdx = 1:param.NumClusters
    for gNBIdx = 1:numCellsPerCluster
        index = (clstIdx - 1)*numCellsPerCluster + gNBIdx;
        % gNB position
        %YXC begin
        pos = [cellRadius*(xClusterCentre(clstIdx)+xGNBPos(gNBIdx)+x), cellRadius*(yClusterCentre(clstIdx)+yGNBPos(gNBIdx)+y), param.gNBHeight];
        %YXC end
        gNBCoordinates(index, :) = pos;
        
        %MXC
        %cell position
        cellPos = cellRadius*[xClusterCentre(clstIdx)+xCellPos(gNBIdx)+x, yClusterCentre(clstIdx)+yCellPos(gNBIdx)+y];
        cellCoordinates(index, :) = cellPos;
        %MXC
    end
end


% UE positions initialization
ueCoordinates = cell(numCellsPerCluster, 1);
for idx = 1:numCellsPerCluster
    ueCoordinates{idx, 1} = zeros(param.NumUEsCell, 3);
end
numGNBs = numCellsPerCluster*param.NumClusters;
% Distance of vertices from centre of a cell
% MXC: This is later used to check if the UE falls into the hexagon; Need 7 points to draw a closed hexagon
xDistance = cellRadius*cosd(0:60:360);
yDistance = cellRadius*sind(0:60:360);

%MXC: Use gNBIdx as index here is okay because each gNb corresponds to a particular cell
for gNBIdx = 1:numGNBs
    xCellPos = cellCoordinates(gNBIdx, 1);
    yCellPos = cellCoordinates(gNBIdx, 2);
    % Vertices of a cell in a cluster
    x = xCellPos + xDistance;
    y = yCellPos + yDistance;
    for ueIdx = 1:param.NumUEsCell
        flag = 1;
        while flag
            %MXC
            %fixed issue of UEs distributed closer to the center of the cell
            %r = cellRadius*rand;
            r = cellRadius*sqrt(rand);
            %MXC
            theta = 360*rand;
            % Get a random point
            px = r*cosd(theta) + xCellPos;
            py = r*sind(theta) + yCellPos;
            % Check if the random point is inside the cell
            [in, on] = inpolygon(px, py, x, y);
            
            %check if satisfy min gNB UE separation
            UE_gNB_distance = norm(gNBCoordinates(gNBIdx,1:2)-[px py]);
            
            if in && ~on && UE_gNB_distance>param.minUEgNBDistance
                flag = 0;
            end
            
        end
        % Store position of UE
        %YXC begin
        pos = [px, py, param.UEHeight];
        %YXC end
        ueCoordinates{gNBIdx, 1}(ueIdx, :) = pos;
    end
end



