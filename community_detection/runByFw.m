function [thetaHat, Tracking] = runByFw(X, N, varargin) 
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
    p.addOptional('maxIters', 1000);
    p.addOptional('lambda', -2);
    p.addOptional('isSparse', 0);
    p.addOptional('isSecondOrder', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    ops = p.Results;

    fprintf('Running Eigendecomposition ...\n')
    options = struct;
    options.maxIters = 120;
    L = size(X, 1);
    options.verbose = 1;
    options.N = N;
    options.certTol = -1;
    options.tol = -1;
    options.XSparse = false; 
    options.objTol = 0.00;
    options.estimatedN = 4*N;
    options.innerVerboseInterval = int64(5000);
    options.backend = 'mex';
    options.epsilon = 1e-5;
    options.numThreads = 8;
    options.iterStartAt = 1;
    % options.curvatureConst = 5000;
    options.dualGapThreshold = 0.5; 

    originalX = X;
    [originalV, ~] = eigs(X, N);

    fprintf('Before preprocessing size: %s \n', num2str(size(X, 1)));
    [V, D] = eigs(X, N);
    options.firstIterMethod = 'NONE';
    L = size(originalX, 1);

    fprintf('Original size: %s \n', num2str(size(V)));
    % colEnergy = vecnorm(V, 2, 2).^2;
    % sortColEnergy = sort(colEnergy, 'descend');
    % totalEnergy = sum(colEnergy);
    % sumEnergy = 0;
    % threshold = 0;
    % for i=1:L
    %     sumEnergy = sumEnergy + sortColEnergy(i);
    %     if sumEnergy > 0.99*totalEnergy
    %         threshold = sortColEnergy(i);
    %         break
    %     end
    % end
    % zeroFiltering = find(colEnergy < threshold);
    % V(zeroFiltering, :) = [];
    V = truncateData(V, 0.99);

    fprintf('After preprocessing size: %s \n', num2str(size(V)));
    X = V';
    
    options.lambda = 0; %1e-6;

    fprintf('Done\n')
    fprintf('Running FW ...\n')
    [hatC, Tracking] = fw(X, options, 'originalAdjacency', originalX);
    Tracking.fwX = X;
    fprintf('Done\n')
    [~, lambdaHat] = maxk(vecnorm(hatC, Inf, 2), N, 1);
    
    hatA = X(:, lambdaHat);
    disp('hatA: ');
    disp(hatA);
    hatS = originalV / hatA';
    tre = 1e-12;
    hatS(hatS<tre) = 0;
    hatS = normalize_row_l1_(hatS);
    thetaHat = hatS;
    Tracking.thetaHat = thetaHat;
    Tracking.lambdaHat = lambdaHat;
end

function [ theta_l1 ] = normalize_row_l1_( theta )

d_theta = sum(theta,2);
S = find(d_theta>0);
theta_l1 = theta;
theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));

S1 = find(d_theta==0);
theta_l1(S1, :) = ones(numel(S1), size(theta, 2)) .* mean(theta_l1(S, :));

end
