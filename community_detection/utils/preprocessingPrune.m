function [S, d] = preprocessingPrune(X, r, q, epsilon, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    %% Prune function
    X = real(X);
    if nargin==1
        r = 10;
        q = 0.75;
        epsilon = 0.05;
    end

    vnorm = vecnorm(X');
    S0=find(vnorm>quantile(vnorm,q));

    if isempty(S0)
        S = [];
        d = 0;
    else
        [~,D]=knnsearch(X,X(S0,:),'K',r);
        d=mean(D,2);
        S=S0(d>quantile(d,1-epsilon));
    end

end
