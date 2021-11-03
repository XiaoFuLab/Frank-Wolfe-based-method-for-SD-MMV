function [C] = fwCoreMatlabTunning(X, initC, lambda, iterStart, epsilon, maxIters,...
        innerVerboseInterval, objTol, verbose, dualGapThreshold) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    fprintf('DEBUG: normX = %f \n', norm(X, 'fro'));

    L = size(initC, 1);
    C = initC;
    XT = X';
    relative_change = 1e6;
    current_obj = 1e6;

    XTX = XT*X;
    for iter=iterStart:maxIters
        CNew = C;
        stepSize = 2.0/double(iter + 1);
        
        dualGap = 0;


        
        % G1 = 2*XT*(X*C - X);
        G1 = 2*(XTX*C - XTX);

        maxRowC = max(C, [], 2);
        denominator = sum(exp((C-maxRowC)/epsilon), 2);
        nominator = exp((C-maxRowC)/epsilon);
        G2 = nominator ./ denominator;
        G = G1 + lambda*G2;
        [~, maxIndex] = min(G);
        % disp(G);
        S = sparse(maxIndex, 1:L, ones(L, 1), L, L);
        S = full(S);

        dualGap = dualGap + trace((C-S)'*G);
        C = C + stepSize*(S-C);

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
            fprintf('\t Dual gap: %e \n', dualGap);
            fprintf('\n');
        end

        if relative_change < objTol
            ending = 1;
            break
        end

        if dualGap < dualGapThreshold
            ending = 1;
            break
        end
    end
end
