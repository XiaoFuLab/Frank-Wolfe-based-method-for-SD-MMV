function [out] = phiC(C, epsilon, varargin) 
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    L = size(C, 2);
    C = C/epsilon;
    % C = exp(C);
    H2 = max(C, [], 2);

    H1 = zeros(L, 1);
    for i=1:L
        H1(i) = sum(exp(C(i, :) - H2(i)));
    end
    % H = sum(C, 2);
    
    % H = H/L;

    out = epsilon * sum(log(H1)+H2 - log(L));
end


