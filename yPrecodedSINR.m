function [sinr, F] = yPrecodedSINR(H, nVar, P)
    % Returns sinr, SINR of a MIMO channel for unit-power symbols and
    % F, collection of precoders for each subcarrier and symbol
    %
    % H - An array of size K-by-L-by-R-by-nLayer, where K is number of
    % subcarriers, L is the number of symbols, R is the number of receive
    % antennas, nLayer is the number of layers or transmission rank. In
    % other words, H is the channel matrix per subcarrier per symbol.
    %
    % nVar - Noise variance in linear scale.
    %
    % P - Precoder. When specified as 'MF' or 'ZF', the precoder is
    % calculated by matched filtering or zero-forcing respectively; when
    % specified as a matrix/vector, it is treated as the precoder itself.
    %
    % sinr - SINR of type scalar in linear scale
    
    % Construct function handles for calculating unnormalized precoder
    MFPrecoder = @(H) H';
    ZFPrecocer = @(H) H'/(H*H');
    DesigPrecoder = @(H) P;
    
    % Parse input
    if ischar(P)||isstring(P)
        switch P
            case 'MF' % Matched filtering
                Pcdr = MFPrecoder;
            case 'ZF' % Zero-forcing
                Pcdr = ZFPrecocer;
            case 'MMSE'
                % See equation 31 in Peel, C. B., Hochwald, B. M., & Swindlehurst, 
                % A. L. (2005). A vector-perturbation technique for near-capacity 
                % multiantenna multiuser communication-part I: Channel inversion 
                % and regularization. IEEE Transactions on Communications, 53(1), 
                % 195–202. https://doi.org/10.1109/TCOMM.2004.840638
                K = size(H,3);
                alpha_MMSE = K*nVar;
                Pcdr = @(H) H'/(H*H' + alpha_MMSE*eye(K));
            otherwise
                error('Not implement yet')
        end
    elseif isa(P, 'double')
        Pcdr = DesigPrecoder;
    else
        error('"P" must be of either char array/string or vector/matrix')
    end
    
    % Construct function handles for calculating effective sum sinr
    f = @(x) log2(1+x);
    finv = @(x) 2.^x - 1;
    effsum = @(X) finv(sum(f(X),'omitnan'));
    
    % Calculate the SINR values as per LMMSE estimation method and precoders per
    % subcarrier per symbol
    F = zeros(size(H,1,2,4,3)); % Initialize collection of precoders
    [NumSubcar, NumSym] = size(H,1,2); % Get the numbers of subcarriers and symbols
    sinrGrid = zeros(NumSubcar, NumSym);
    for k = 1:NumSubcar
        for l = 1:NumSym
            Htemp = H(k,l,:,:);
            Htemp = reshape(Htemp,size(H,3,4));
            W_ = Pcdr(Htemp); % Get precoder
            if ismatrix(W_)
                W = W_;
            elseif ndims(W_)==4 % if P is a collection of precoders for each subcarrier and each symbol
                W = squeeze(W_(k,l,:,:));
            else
                error('Invalid P.');
            end
            W = W/norm(W,'fro'); % Normalize the precoder
            F(k,l,:,:) = W;
            noise = nVar*eye(size(W,2)); % Noise variance multiplied by identity matrix
            den = noise/((W'*Htemp')*Htemp*W+noise);
            LayeredSinr = real(((1./diag(den))-1));
            sinrGrid(k,l) = effsum(LayeredSinr);
        end
    end
    
    % Function handles for calculating entropic mean
    entrmean = @(X) finv(mean(f(X),'omitnan'));
    
    sinr = entrmean(entrmean(sinrGrid));
    
end