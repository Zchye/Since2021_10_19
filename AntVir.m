function popedGrid = AntVir(InputStruct)
    % Antenna virtualization
    % InputStruct should contain following fields: 
    %
    % PortDims Dimensions of antenna ports, an array of dv-by-dh, where dv
    % is the size of vertical dimension, dh is that of horizontal one,
    % pertaining to ports
    % AntDims Dimensions of antenna elements (single panel only), an array
    % of dv-by-dh, pertaining to antenna elements
    % AntGrid Antenna element resource grid to be populated
    % NumPorts Number of antenna ports for certain reference signal
    % Sym Symbols
    % SymLinInd The linear indices of symbols in port resource grid
    
    % Parse InputStruct
    PortDims = InputStruct.PortDims(:)'; % Reshape to row vector
    AntDims = InputStruct.AntDims(:)'; % Reshape to row vector
    if any(mod(AntDims, PortDims)) % If PortDims does not divide AntDims
        error('The size of port dimensions ([%s]) must divides the size of antenna element dimensions ([%s])',...
            num2str(PortDims),num2str(AntDims));
    end
    AntGrid = InputStruct.AntGrid;
    NumPorts = InputStruct.NumPorts;
    Sym = InputStruct.Sym;
    SymLinInd = InputStruct.SymLinInd;
    
    NumAnt = size(AntGrid,3); % Number of antenna elements
    NumAntSubArray = NumAnt/NumPorts; % Number of antenna elements in a subarray
    SymLinInd = SymLinInd(:)'; % Reshape the array of linear indices to a row vector
    Sym = reshape(Sym, 1, []); % Reshape the array of symbols to a row vector
    
    % Convert 1-based linear indices to subscripts
    PortGridSize = [size(AntGrid,1,2), NumPorts]; % Size of antenna port resource grid
    [PortSub1,PortSub2, PortSub3] = ind2sub(PortGridSize, SymLinInd); % The conversion
    
    % Create map container that maps indices of port to indices of antenna
    % elements
    PEMap = P2EMapping(PortDims, AntDims);
    
    % Normalize energy by NumAntSubArray
    NormSym = Sym/sqrt(NumAntSubArray);
    
    % Antenna ports to antenna elements mapping
    popedGrid = AntGrid;
    for ii = 1:length(SymLinInd)
        AntPages = PEMap(PortSub3(ii));
        popedGrid(PortSub1(ii), PortSub2(ii), AntPages) = NormSym(ii)*ones(1,1,length(AntPages));
    end
end

function PEMap = P2EMapping(PortDims, AntDims)
    % Returns a map container that maps ports to antenna elements
    SubArrayDims = AntDims./PortDims; % Size of dimensions of an subarray
    SubArraySize = prod(SubArrayDims); % Number of elements in a subarray
    [SubArraySub1, SubArraySub2] = ind2sub(SubArrayDims, 1:SubArraySize); % Relate subarray indices to subscripts
    XPolJump = prod(AntDims); % Index jump between cross polarizations
    NumPorts = prod(PortDims); % Number of ports
    [PortSub1, PortSub2] = ind2sub(PortDims, 1:NumPorts); % Relate port indices to subscripts
    PEMap = containers.Map(1:NumPorts,cell(1,NumPorts)); % Initialize map containers
    for portIdx = 1:NumPorts
        ThisSub1 = (PortSub1(portIdx) - 1)*SubArrayDims(1) + SubArraySub1;
        ThisSub2 = (PortSub2(portIdx) - 1)*SubArrayDims(2) + SubArraySub2;
        tempInd = sub2ind(AntDims, ThisSub1(:), ThisSub2(:));
        tempInd = tempInd(:)';
        PEMap(portIdx) = [tempInd, tempInd + XPolJump];
    end
end