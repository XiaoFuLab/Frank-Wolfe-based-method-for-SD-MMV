function [C, Tracking] = fwCoreMatlab(X, initC, lambda, iterStart, epsilon, maxIters,...
        innerVerboseInterval, objTol, verbose, dualGapThreshold, varargin) 
% ----------------------------------------------------------------------
% 
% FW-based algorithm for solving
% minimize_{C}   ||X-XC||_F^2 + \lambda*\varPhi_{\mu}(C)
% subject to     C \geq 0, 1^T C = 1^T
%
% This is Matlab-based implementation. This version is mainly for inspecting 
% method's behaviour. The implementation is also cleaner and easier to read. 
% It produces an optimal C and various tracking information during updating 
% process via variable `Tracking`.
%
% A faster version is provided in C++ via a mex file without tracking information.
%
% Papams:
% - X: data input
% - initC: initilization for C
% - lambda: lambda in the objective function above
% - iterStart: standard Frank-Wolfe start with iterStart=0, although a warm version can use iterStart>0
% - epsilon: hyperparameter of the regularization term
% - maxIters: 
% - innerVerboseInterval: just logging interval
% - objTol: a threshold in relative change of objective function 
% - verbose: logging or not
% - dualGapThreshold: another criterion for early stopping
% Output:
% - C: optimal C
% - Tracking: a structure with some information collecting during the process.

% Author: Tri Nguyen (nguyetr9@oregonstate.edu)
% Reference: 
% Nguyen, T., Fu, X. and Wu, R., 2021. Memory-efficient convex optimization for 
% self-dictionary separable nonnegative matrix factorization: A frank-wolfe approach. 
% arXiv preprint arXiv:2109.11135. (Submitted to TSP, accepted, 2022.)
% ----------------------------------------------------------------------

    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', false);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    Tracking = struct;
    Tracking.obj = [];
    Tracking.nnz = [];
    Tracking.lInfNorm = [];

    N = size(initC, 1);
    C = initC;
    XT = X';
    relative_change = 1e6;
    current_obj = 1e6;

    for iter=iterStart:maxIters
        CNew = C;
        stepSize = 2.0/double(iter + 1);
        
        sumCe1 = zeros(N, 1);
        sumCe2 = zeros(N, 1);
        for k=1:N
            muRowC = C(k, :) ./ epsilon;
            maxMuRowC = max(muRowC);
            muRowC = muRowC - maxMuRowC;
            sumCe1(k) = sum(full(exp(muRowC)));
            sumCe2(k) = maxMuRowC;
        end

        g1Mean = 0;
        g2Mean = 0;
        dualGap = 0;

        for k=1:N

            if verbose && mod(k, innerVerboseInterval) == 0
                fprintf('Inner iteration progress: %d/%d \n', k, N);
            end
            ck = C(:, k);

            g1 = 2*(XT*(X*ck - X(:, k)));
            g1Mean = g1Mean + mean(g1);
            
            g2 =  exp(ck./epsilon - sumCe2) ./ sumCe1;
            g2Mean = g2Mean + mean(g2);

            g = g1 + lambda*g2;
            [vvv, ind] = min(g);
            
            s = zeros(N, 1); s(ind) = 1.0;
            ckNew = ck + stepSize*(s - ck);
            CNew(:, k) = ckNew;
            dualGap = dualGap + (ck-s)'*g;
        end
        if options.debug
            c_max = max(C, [], 2);
            second_part = sum(c_max - epsilon*log(N) +  epsilon * log(sum(exp((C - c_max)/epsilon), 2)));
            Tracking.obj(end+1) = norm(X - X*C, 'fro')^2 + lambda*second_part;
        end
        Tracking.nnz(end+1) = nnz(C);
        Tracking.lInfNorm(end+1) = sum(vecnorm(C, Inf, 2));
        C = CNew;

        for k=1:N
            muRowC = C(k, :) ./ epsilon;
            maxMuRowC = max(muRowC);
            muRowC = muRowC - maxMuRowC;
            sumCe1(k) = sum(full(exp(muRowC)));
            sumCe2(k) = maxMuRowC;
        end
        fittingError = norm(X - X*C, 'fro')^2;
        regValue = epsilon*sum(log(sumCe1) + sumCe2 - log(N) );
        new_obj = fittingError + lambda*regValue;
        % new_obj = 0;
        relative_change = abs(new_obj - current_obj)/current_obj;
        current_obj = new_obj;

        if verbose
            fprintf('------------ Iter: %d ------------\n', iter);
            fprintf('\t stepSize: %f\n', stepSize);
            fprintf('\t obj rel change: %e \n', relative_change);
            fprintf('\t Fitting error: %e \n', fittingError);
            fprintf('\t Reg value: %e \n', regValue);
            fprintf('\t Current obj: %e \n', current_obj);
            fprintf('\t g1 mean: %e \n', g1Mean/N);
            fprintf('\t g2 mean: %e \n', g2Mean/N);
            fprintf('\t Dual gap: %e \n', dualGap);
            fprintf('\n');
        end

        if relative_change < objTol
            ending = 1;
            break
        end

        if iter>1 && dualGap < dualGapThreshold
            ending = 1;
            break
        end
    end
end
