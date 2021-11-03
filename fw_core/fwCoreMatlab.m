function [C, Tracking] = fwCoreMatlab(X, initC, lambda, iterStart, epsilon, maxIters,...
        innerVerboseInterval, objTol, verbose, dualGapThreshold, varargin) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

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
    fprintf('DEBUG: normX = %f', norm(X, 'fro'));

    L = size(initC, 1);
    C = initC;
    XT = X';
    relative_change = 1e6;
    current_obj = 1e6;
    % if iterStart ~= 1 && iterStart ~= 2
    %     error('iterStart must be either 1 or 2!');
    % end

    for iter=iterStart:maxIters
        CNew = C;
        stepSize = 2.0/double(iter + 1);
        
        sumCe1 = zeros(L, 1);
        sumCe2 = zeros(L, 1);
        for k=1:L
            muRowC = C(k, :) ./ epsilon;
            maxMuRowC = max(muRowC);
            muRowC = muRowC - maxMuRowC;
            sumCe1(k) = sum(full(exp(muRowC)));
            sumCe2(k) = maxMuRowC;
        end

        g1Mean = 0;
        g2Mean = 0;
        dualGap = 0;

        for k=1:L

            if verbose && mod(k, innerVerboseInterval) == 0
                fprintf('Inner iteration progress: %d/%d \n', k, L);
            end
            ck = C(:, k);

            g1 = 2*(XT*(X*ck - X(:, k)));
            g1Mean = g1Mean + mean(g1);
            
            g2 =  exp(ck./epsilon - sumCe2) ./ sumCe1;
            g2Mean = g2Mean + mean(g2);

            g = g1 + lambda*g2;
            [vvv, ind] = min(g);
            % ============ START DEBUG ============
            % ============  END DEBUG  ============
            
            s = zeros(L, 1); s(ind) = 1.0;
            ckNew = ck + stepSize*(s - ck);
            CNew(:, k) = ckNew;
            dualGap = dualGap + (ck-s)'*g;
        end
        if options.debug
            c_max = max(C, [], 2);
            second_part = sum(c_max - epsilon*log(L) +  epsilon * log(sum(exp((C - c_max)/epsilon), 2)));
            Tracking.obj(end+1) = norm(X - X*C, 'fro')^2 + lambda*second_part;
        end
        Tracking.nnz(end+1) = nnz(C);
        Tracking.lInfNorm(end+1) = sum(vecnorm(C, Inf, 2));
        C = CNew;

        for k=1:L
            muRowC = C(k, :) ./ epsilon;
            maxMuRowC = max(muRowC);
            muRowC = muRowC - maxMuRowC;
            sumCe1(k) = sum(full(exp(muRowC)));
            sumCe2(k) = maxMuRowC;
        end
        fittingError = norm(X - X*C, 'fro')^2;
        regValue = epsilon*sum(log(sumCe1) + sumCe2 - log(L) );
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
            fprintf('\t g1 mean: %e \n', g1Mean/L);
            fprintf('\t g2 mean: %e \n', g2Mean/L);
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
