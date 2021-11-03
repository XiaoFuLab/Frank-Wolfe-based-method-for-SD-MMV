function [thetaHat, Tracking] = runByFg(X, N, isPruned, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    fprintf('Running Eigendecomposition ...\n')
    originalX = X;
    [V, D] = eigs(X, N);
    originalV = V;
    L = size(originalX, 1);

    if isPruned
        fprintf('Original size: %s \n', num2str(size(V)));
        colEnergy = vecnorm(V, 2, 2).^2;
        sortColEnergy = sort(colEnergy, 'descend');
        totalEnergy = sum(colEnergy);
        sumEnergy = 0;
        threshold = 0;
        for i=1:L
            sumEnergy = sumEnergy + sortColEnergy(i);
            if sumEnergy > 0.99*totalEnergy
                threshold = sortColEnergy(i);
                break
            end
        end
        zeroFiltering = find(colEnergy < threshold);

        V(zeroFiltering, :) = [];
        fprintf('After preprocessing size: %s \n', num2str(size(V)));
    end
    X = V';

    fprintf('Done\n')
    fprintf('Running FG ...\n')
    [~, K, ~, ~, ~, ~, Tracking] = fgnsr(X, N, 'maxiter', 20, ...
        'debug', 0, 'verbose', true, 'log_interval', 1);
    lambdaHat = K;
    
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
