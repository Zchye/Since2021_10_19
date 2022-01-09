function sinr = yPrecodedSINR(H,sigma,W)
    % Calculates effective SINR of a MIMO channel
    
    noise = sigma^2*eye(size(W,2)); % Noise variance per layer
    den = noise/( (W'*H')*H*W+noise );
    
    f = @(x) log2(1 + x); % Function handle for calculating entropy
    finv = @(x) 2^x - 1; % Inverse of f
    x = (1./diag(den))-1; % A vector whose elements are stream sinrs
    sinr = real(finv(sum(f(x)))); % Effective sinr
    
end