function popedGrid = AntVir(AntGrid, NumPorts, Sym, SymLinInd)
    % Antenna virtualization
    % ANTGRID Antenna element resource grid to be populated
    % NUMPORTS Number of antenna ports for certain reference signal
    % SYM Symbols
    % SYMLININD The linear indices of symbols in port resource grid
    
    NumAnt = size(AntGrid,3); % Number of antenna elements
    NumAntSubArray = NumAnt/NumPorts; % Number of antenna elements in a subarray
    
    % Convert 1-based linear indices to subscripts
    PortGridSize = [size(AntGrid,1,2), NumPorts]; % Size of antenna port resource grid
    [PortSub1,PortSub2, PortSub3] = ind2sub(PortGridSize, SymLinInd); % The conversion
    
    % Normalize energy by NumAntSubArray
    NormSym = Sym/sqrt(NumAntSubArray);
    
    % Antenna port mapping
    popedGrid = AntGrid;
    for ii = 1:length(SymLinInd)
        for PortIdx = 1:NumPorts
            if PortSub3(ii) == PortIdx
                % Populate a copy of one port to adjacent
                % antennas
                for subp = 0:(NumAntSubArray-1)
                    popedGrid(PortSub1(ii), PortSub2(ii), NumAntSubArray*PortSub3(ii) - subp) = NormSym(ii);
                end
            end
        end
    end
end