function [C] = fwCoreMatlab__deprecated(X, initC, lambda, iterStart, epsilon, maxIters,...
        innerVerboseInterval, objTol, verbose) 
% ----------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    fprintf('DEBUG: normX = %f', norm(X, 'fro'));

    L = size(initC, 1);
    C = initC;
    XT = X';
    relative_change = 1e6;
    current_obj = 1e6;
    if iterStart ~= 1 && iterStart ~= 2
        error('iterStart must be either 1 or 2!');
    end

    for iter=iterStart:maxIters
        CNew = C;
        stepSize = 2.0/double(iter + 1);
        
        sumCe = zeros(L, 1);
        for k=1:L
            sumCe(k) = full(sum(exp(C(k, :)./epsilon)));
        end

        for k=1:L

            if verbose && mod(k, innerVerboseInterval) == 0
                fprintf('Inner iteration progress: %d/%d \n', k, L);
            end
            ck = C(:, k);

            g1 = 2*(XT*(X*ck - X(:, k)));
            g2 =  exp(ck./epsilon) ./ sumCe;
            g = g1 + lambda*g2;
            [vvv, ind] = min(g);
            s = zeros(L, 1); s(ind) = 1.0;
            tmp = ck + stepSize*(s - ck);
            CNew(:, k) = tmp;
        end
        C = CNew;
        fittingError = norm(X - X*C, 'fro')^2;
        new_obj = fittingError + lambda*epsilon*sum(log(sumCe/L));
        relative_change = abs(new_obj - current_obj)/current_obj;
        current_obj = new_obj;
        if relative_change < objTol
            ending = 1;
            break
        end

        if verbose
            fprintf('------------ Iter: %d ------------\n', iter);
            fprintf('\t stepSize: %f\n', stepSize);
            fprintf('\t obj rel change: %f \n', relative_change);
            fprintf('\t Current obj: %f \n', current_obj);
            fprintf('\n');
        end
    end
end
