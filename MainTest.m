% MainTest for testing heatmap

% Files involved:
%   ConfigChannel.m
%   ChPLPerPixel.m

% Time the program
tStart = tic;

load simParameters.mat

MapCenter = simParameters.GNBPositions(1,:); % The center of the heatmap
CellRad = simParameters.InterSiteDistance/3; % Cell radius
Margin = 10; % The width of margins in the heatmap [m]
% Get the range of the heatmap. For now we focus on the 3 cells at the
%center
Xmin = MapCenter(1) - 1.5 * CellRad - Margin;
Xmax = MapCenter(1) + 2 * CellRad + Margin;
Ymin = MapCenter(2) - sqrt(3) * CellRad - Margin;
Ymax = MapCenter(2) + sqrt(3) * CellRad + Margin;

% Deploy pixels
prepwb = waitbar(0,'Deploying pixels...'); % The waitbar displaying the
%progress of pixel deployment.
Npix = 120; % The number of pixels on each dimension
X = linspace(Xmin, Xmax, Npix);
Y = linspace(Ymin, Ymax, Npix);
% Get the coordinates of vertices of a hexagon relative to its center
Xvert = CellRad*cosd(0:60:360);
Yvert = CellRad*sind(0:60:360);
% Preallocate a cell array to store the indices of pixels that are contianed
%in some hexagons.
% Each cell contains a triple with the first and second entry being the
%indices of X and Y coordinates and the third entry being the site index.
InHexPixels = cell(Npix*Npix,1);
IdxPtr = 1; % Index pointer pointing to the first empty cell in InHexPixels
% Populate InHexPixels
wblen = 0; % waitbar length
for xx = 1:Npix
    for yy = 1:Npix
        for siteIdx = 1:3
            xv = simParameters.CellPositions(siteIdx,1)+Xvert;
            yv = simParameters.CellPositions(siteIdx,2)+Yvert;
            isin = inpolygon(X(xx),Y(yy),xv,yv);
            if isin
                InHexPixels{IdxPtr} = [xx,yy,siteIdx];
                IdxPtr = IdxPtr+1;
                break
            end
            % Update waitbar
            wblen = ((xx-1)*Npix+yy)/(Npix^2);
            waitbar(wblen,prepwb);
        end
    end
end
% Remove empty cells
InHexPixels(IdxPtr:end) = [];
% Pixel deployment finished, close the waitbar
close(prepwb);

ChPLs = ChPLPerPixel(simParameters); % Instantiate ChPLPerPixel
NumGNBAnts = prod(simParameters.GNBTxAntPanelSize); % Number of gNB antennas
SINRs = zeros(Npix,Npix); % To store SINRs of each pixel

% Comput SINR for each pixel
MLStart = tic; % Time the computation, since it could run for hours or even days
wb = waitbar(0,'Starting main loop...');
for ii = 1:IdxPtr-1
    % Loop through all pixels that are contained in some hexagons
    
    % Get the pixel's information
    xx = InHexPixels{ii}(1);
    yy = InHexPixels{ii}(2);
    ThisSite = InHexPixels{ii}(3); % The index of the site where this pixel is in
    PixelPos = [X(xx),Y(yy),simParameters.UEHeight];
    
    % Comput SINR for this pixel
    updateChPL(ChPLs,PixelPos); % Update the channels and pathlosses for this pixel
    NoiseEnergy = 0; % Noise energy temporarily set to zero
    ips = ones(1,NumGNBAnts); % Impulse signal
    ipr = cell(simParameters.NumSitesPerCluster,1); % Cell array for storing impulse response
    % Filter signals coming from each gNB
    for siteIdx = 1:simParameters.NumSitesPerCluster
        if siteIdx == ThisSite % Apply channel and pathloss to the signal from ThisSite
            ipr{siteIdx} = ChPLs(siteIdx,ips);
            SignalEnergy = norm(ipr{ThisSite},'fro')^2; % Get the signal energy
        else % Filter signals from other gNBs (with or without channel)
            ipr{siteIdx} = ChPLs(siteIdx,ips);
        end
    end
    % Calculate the impulse response energy of the signal from each gNB and
    % store in a vector
    iprEnergy = cellfun(@(x) norm(x,'fro')^2, ipr);
    SumEnergy = sum(iprEnergy);
    InterfEnergy = SumEnergy - SignalEnergy; % Get interference energy
    % Get SINRs
    SINRs(xx,yy) = NaN2Nil(SignalEnergy/(InterfEnergy+NoiseEnergy));
    %Update waitbar
    wblen = ii/(IdxPtr-1); % Calculate waitbar length
    elapsedtime = toc(MLStart); % Get elapsed time
    totaltime = elapsedtime/wblen; % Estimate total time for the main loop
    ela_s = seconds(elapsedtime); % Cast the time to datetime data type
    tl_s = seconds(totaltime); % Same as above
    msg = [datestr(ela_s,'dd:HH:MM:SS'),...
        '/',datestr(tl_s,'dd:HH:MM:SS'),...
        ' Pixels:',num2str(ii),'/',num2str(IdxPtr-1)]; % Construct the message
    %to be shown in the waitbar dialogue
    waitbar(wblen,wb,msg);    
end
% Main loop finished, close the waitbar
close(wb);

% Display the total time for computing SINRs
tEnd = toc(tStart);
TE = seconds(tEnd);
TE.Format = 'dd:hh:mm:ss';
fprintf('Total time for %d-by-%d pixels:\n',Npix,Npix)
disp(TE)

% Plot the SINR heatmap and the SINR CDF
SINRs_dB = 10*log10(SINRs); % From linear to dB
SINRs_dB = flip(SINRs_dB',1); % Orient x-axis to East and y-axis to North
% Plot SINR heatmap
figure
image(X,Y,SINRs_dB,'CDataMapping','scaled') % Scale the SINRs(dB) to increase contrast
colorbar
title('SINR Heatmap')
% Plot SINR CDF
yuv = SINRs_dB(:); % Rearrange to a vector to fit the input signature of ecdf()
yuv(yuv==-inf) = [];
figure
ecdf(yuv)
xlabel('SINR')
ylabel('CDF')

function y = NaN2Nil(x)
% Return zero when input is NaN, otherwise simply pass it to the output
    if isnan(x)
        y = 0;
    else
        y = x;
    end
end
