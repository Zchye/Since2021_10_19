function [antsym,antind] = yPRGPrecode(siz,nstartgrid,portsym,portind,F)
    % Modified from 5G toolbox built-in function hPRGPrecode.m
    if nstartgrid
        % Ignore this input argument
    end
    K = siz(1); % K, number of subcarriers
    L = siz(2); % L, number of symbols
    nLayer = size(F,4); % nLayer, numerb of layers
    PortGridSize = [K,L,nLayer]; % Size of port grid
    PortGrid = zeros(PortGridSize); % Initialize port grid
    PortGrid(portind(:)) = portsym(:); % Reconstruct port grid
    BFGrid = zeros(K,L, size(F,3)); % Initialize beamformed grid
    % Apply precoder
    for k = 1:K
        for l = 1:L
            BFGrid(k,l,:,:) = squeeze(F(k,l,:,:))*reshape(PortGrid(k,l,:),[],1);
        end
    end
    
    [antsym, antind] = nrExtractResources(portind, BFGrid);
end