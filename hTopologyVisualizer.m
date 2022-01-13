function hTopologyVisualizer(param)
%hTopologyVisualizer Visualize the network topology for the given configuration
%
%   hTopologyVisualizer(PARAM) creates a macrocell network according to the
%   configuration specified by PARAM. PARAM is a structure
%   with the following fields:
%   GNBPositions      - Position of gNBs in each cell
%   UEPositions       - Positions of UEs in each cell
%   InterSiteDistance - Distance between two adjacent gNBs
%   NumUEsCell        - Number of UEs in each cell

%   Copyright 2020 The MathWorks, Inc.

%MXC
%changed due to different topology
sideLength = param.InterSiteDistance/sqrt(3);
% sideLength = param.InterSiteDistance/3;
%MXC

% Using the screen width and height, calculate figure width and height
resolution = get(0,'ScreenSize');
screenWidth = resolution(3);
screenHeight = resolution(4);
figureWidth = screenWidth * 0.90;
figureHeight = screenHeight * 0.85;
fig = figure('Name', 'Topology Visualization', 'Position', [screenWidth*0.05 screenHeight*0.05 figureWidth figureHeight],...
                'Visible', 'on');
ax = axes(fig);
ax.XLabel.String = 'Meters';
ax.YLabel.String = 'Meters';
hold(ax,'on');
pbaspect(ax, [1 1 1]); 
numgNBs = size(param.GNBPositions, 1);
% Distance of vertices from centre of a cell
xDistance = sideLength*cosd(0:60:360);
yDistance = sideLength*sind(0:60:360);

%MXC
%gNB position does not correspond to cell center any more
 %{
for gNBIdx = 1:numgNBs
    xGNBPos = param.GNBPositions(gNBIdx, 1);
    yGNBPos = param.GNBPositions(gNBIdx, 2);
    % Vertices of a cell in a cluster
    x = xGNBPos + xDistance;
    y = yGNBPos + yDistance;
    if gNBIdx <= 19 % Highlighting centre cluster
        plot(ax, x, y, 'k', 'LineWidth', 2.5, 'Tag', 'line');
    else
        plot(ax, x, y, 'k', 'Tag', 'line');
    end
    
    % Plot gNB position
    text(xGNBPos, yGNBPos, num2str(mod(gNBIdx-1, 19)), 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    for ueIdx = 1:param.NumUEsCell
        % Plot UE positions
        uePlot = plot(ax, param.UEPositions{gNBIdx,1}(ueIdx,1), param.UEPositions{gNBIdx,1}(ueIdx,2),...
                     'r.', 'MarkerSize', 12, 'Tag', ['UE' num2str(ueIdx)]); 
    end
end
legend(uePlot, 'UE Positions');
title('Network Topology Visualization');
hold(ax,'off');
%}
for gNBIdx = 1:numgNBs
    xCellPos = param.CellPositions(gNBIdx, 1);
    yCellPos = param.CellPositions(gNBIdx, 2);
    % Vertices of a cell in a cluster
    x = xCellPos + xDistance;
    y = yCellPos + yDistance;
    if gNBIdx <= 21 % Highlighting centre 2 rings
        plot(ax, x, y, 'k', 'LineWidth', 2.5, 'Tag', 'line');
    else
        plot(ax, x, y, 'k', 'Tag', 'line');
    end
    gNBPlot = plot(ax, param.GNBPositions(gNBIdx,1), param.GNBPositions(gNBIdx,2), 'b.', 'MarkerSize', 24);
    
    % Plot cell position
    text(xCellPos, yCellPos, num2str(mod(gNBIdx-1, 57)), 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    for ueIdx = 1:param.NumUEsCell
        % Plot UE positions
        uePlot = plot(ax, param.UEPositions{gNBIdx,1}(ueIdx,1), param.UEPositions{gNBIdx,1}(ueIdx,2),...
                     'r.', 'MarkerSize', 6, 'Tag', ['UE' num2str(ueIdx)]); 
    end
end
legend([uePlot gNBPlot], {'UE Positions', 'gNB Positions'});
title('Network Topology Visualization');
hold(ax,'off');
%MXC
end