function [antsym,antind] = yPRGPrecode(siz,nstartgrid,portsym,portind,F)
    %
    if nstartgrid
        % Ignore this input argument
    end
    [K,L] = size(siz,1,2); % K, number of subcarriers; L, number of symbols
    nLayer = size(F,4); % nLayer, numerb of layers
    PortGridSize = [K,L,nLayer]; % Size of port grid
    PortGrid = zeros(PortGridSize); % Port grid
    PortGrid(portind(:)) = portsym(:); % Reconstruct port grid
    
end