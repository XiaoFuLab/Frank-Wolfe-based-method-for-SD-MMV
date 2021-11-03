function [C, Tracking] = fw(X, varargin) 
% -------------------------------------------------------------------------
%% frank_wolfe utilises Self-Dictionary to solve NMF under pure pixel 
%   assumption
%% Params:
%
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('objTol', 0);
    p.addOptional('maxIters', 500);
    p.addOptional('lambda', 1);
    p.addOptional('N', -1);
    p.addOptional('firstIterMethod', 'NONE');
    p.addOptional('innerVerboseInterval', 1000);
    p.addOptional('backend', '.');
    p.addOptional('numThreads', 8);
    p.addOptional('epsilon', 1e-8);
    p.addOptional('estimatedN', 30);
    p.addOptional('degreeOfBelieve', 1); % From >0 to 1, 1 means you absolutely believe that data is perfectly followed our model (ofcourse having noise)
    p.addOptional('iterStartAt', 1); 
    p.addOptional('dualGapThreshold', 1e-5); 
    p.addOptional('curvatureConst', 1200); 
    p.addOptional('originalAdjacency', nan); 
    p.addOptional('initC', nan); 

    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    [d, L] = size(X);

    Tracking = struct;
    Tracking.obj = zeros(options.maxIters, 1);
    Tracking.time = zeros(options.maxIters, 1);
    Tracking.stepSize = zeros(options.maxIters, 1);
    Tracking.endingCode = 0;
    Tracking.checkpoint = [];
    if options.debug
        Tracking.cert = -1*ones(L, options.maxIters);
    end

    iterDuration = 0;
    C = sparse(L, L);
    epsilon = 1e-2;

    nnzX = nnz(X)/prod(size(X));
    if options.verbose
        fprintf('X nnz ratio: %f \n', nnzX);
    end

    XT = X';
    lambdaInit = options.lambda;
    if options.verbose
        fprintf('Determining lambda ... \n');
    end
    [lambda, approxError, hatC_] = autoDetermineLambda(X, options.N, lambdaInit, epsilon, options.degreeOfBelieve, 'verbose', options.verbose);
    Tracking.lambda = lambda;
    if options.verbose
        fprintf('Lambda: %e - ApproxError: %e \n', ...
            lambda, approxError);
    end
    if options.dualGapThreshold == -1
        options.dualGapThreshold = approxError*5;
    end

    if ~strcmp(options.firstIterMethod, 'NONE')
        if strcmp(options.firstIterMethod, 'SPA')
            C = hatC_;
            % ============ START DEBUG ============
             [~, lambdaHat] = maxk(vecnorm(C, Inf, 2), options.N, 1);   
             disp('At initilization, FW picks:')
             disp(lambdaHat(:)');
            % ============  END DEBUG  ============
            
        elseif strcmp(options.firstIterMethod, 'SPACL_W1')
            [theta, tmp ]= runBySpacl2(options.originalAdjacency, options.N);
            originalL = size(theta, 1);
            ind = 1:originalL;
            ttt = zeros(1, originalL);
            ttt(tmp.lambdaHat) = 1;
            ttt(tmp.I_outlier) = [];
            spaLambdaHat = find(ttt);
            assert(numel(spaLambdaHat) == numel(tmp.lambdaHat));

            for i=1:L
                C(spaLambdaHat(mod(i-1, options.N)+1), i) = 1;
            end
        elseif strcmp(options.firstIterMethod, 'SPACL2')
            [theta, tmp] = runBySpacl2(options.originalAdjacency, ...
                options.N);
            theta = theta';
            C(tmp.lambdaHat, :) = theta;
        elseif strcmp(options.firstIterMethod, 'RANDOM')
            C = dirichletRand(ones(1, L), L);
        elseif strcmp(options.firstIterMethod, 'PREDEFINED')
            C = options.initC;
        else
            error('Mistyped');
        end
    end
    
    % ==============================
    % PICKING WARM START STEP SIZE
    % ==============================
    sumCe1 = zeros(L, 1);
    sumCe2 = zeros(L, 1);
    for k=1:L
        muRowC = C(k, :) ./ epsilon;
        maxMuRowC = max(muRowC);
        muRowC = muRowC - maxMuRowC;
        sumCe1(k) = sum(full(exp(muRowC)));
        sumCe2(k) = maxMuRowC;
    end
    dualityGap = 0;
    for k=1:L
        ck = C(:, k);

        g1 = 2*(XT*(X*ck - X(:, k)));
        g2 =  exp(ck./epsilon - sumCe2) ./ sumCe1;

        g = g1 + lambda*g2;
        [vvv, ind] = min(g);
        s = zeros(L, 1); s(ind) = 1.0;
        dualityGap = dualityGap + g'*(ck - s);
    end
    % ==============================
    % ENDING PICKING WARM START STEP SIZE
    % ==============================
    if options.verbose
        fprintf('Duality gap at initilization: %f \n', dualityGap);
    end

    if strcmp(options.iterStartAt, 'warm')
        iterStart = round(2*options.curvatureConst/dualityGap);
    elseif strcmp(options.iterStartAt, 'heuristic')
        iterStart = 1/sqrt(1/L * approxError^2);
    elseif strcmp(options.iterStartAt, 'heuristic2')
        iterStart = sqrt(options.N/(1/L * approxError^2));
    else
        iterStart = options.iterStartAt;
    end

    % ============ START DEBUG ============
     % maxIters = max(options.maxIters, iterStart+1);   
    maxIters = options.maxIters;
    % ============  END DEBUG  ============
    

    if maxIters>iterStart
        if strcmp(options.backend,  'matlab')
            [C, InnerTracking] = fwCoreMatlab(X, C, lambda, iterStart, options.epsilon, ...
                maxIters, options.innerVerboseInterval, ...
                options.objTol, options.verbose, options.dualGapThreshold, 'debug', options.debug);
            Tracking.InnerTracking = InnerTracking;
        elseif strcmp(options.backend,  'mex')
            estimatedN = int64(options.estimatedN);
            numThreads = options.numThreads;
            if issparse(X)

                if options.verbose
                    fprintf('Running mex with sparse version, since nnz(X)=%d \n', nnz(X));
                end
                C = fwCoreSparseMex(X, C, lambda, iterStart, ...
                    options.epsilon, maxIters, ...
                    options.innerVerboseInterval, ...
                    options.objTol, options.verbose, ...
                    estimatedN, numThreads, options.dualGapThreshold);
            else
                C = fwCoreMex(X, C, lambda, iterStart, options.epsilon, ...
                    maxIters, options.innerVerboseInterval, ...
                    options.objTol, options.verbose, ...
                    estimatedN, numThreads, options.dualGapThreshold);
            end
        elseif strcmp(options.backend, 'matlabTunning')
            C = fwCoreMatlabTunning(X, C, lambda, iterStart, ...
                options.epsilon, maxIters, options.innerVerboseInterval, ...
                options.objTol, options.verbose, options.dualGapThreshold);
        else
            error('backend is either mex or matlab');
        end
    end

end

function [lambda, approxError, hatC_] = autoDetermineLambda(X, N, lambdaInit, epsilon, degreeOfBelieve, varargin) 
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

    L = size(X, 2);
    if lambdaInit == -1
        error('Deprecated')
    elseif lambdaInit == -2
        error('Deprecated')
        % spaLambdaHat = SPAselect(X, N);
        % [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, ...
        %     'maxIters', 100);
        % approxError = norm(X - X(:, spaLambdaHat)*CSpaHat, 'fro')^2;
        % hatC_  = sparse(L, L);
        % hatC_(spaLambdaHat, :) = CSpaHat;
        % s = phiC(hatC_, epsilon);
        % lambda = approxError/s/10;
    elseif lambdaInit == -3
        spaLambdaHat = SPAselect(X, N);
        [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, ...
            'maxIters', 100);
        approxError = norm(X - X(:, spaLambdaHat)*CSpaHat, 'fro');
        L = size(X, 2);
        hatC_  = sparse(L, L);
        hatC_(spaLambdaHat, :) = CSpaHat;
        s = phiC(hatC_, epsilon);
        lambda = approxError/N;
    elseif lambdaInit == -4
        spaLambdaHat = SPAselect(X, N);
        [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, ...
            'maxIters', 100);
        approxError = norm(X - X(:, spaLambdaHat)*CSpaHat, 'fro');
        L = size(X, 2);
        hatC_  = sparse(L, L);
        hatC_(spaLambdaHat, :) = CSpaHat;
        s = phiC(hatC_, epsilon);
        lambda = 4*(approxError/s)/L;
    elseif lambdaInit == -5
        % spaLambdaHat = SPAselect(X, N);
        % [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, ...
        %     'maxIters', 100);
        % approxError = norm(X - X(:, spaLambdaHat)*CSpaHat, 'fro');
        % L = size(X, 2);
        % hatC_  = sparse(L, L);
        % hatC_(spaLambdaHat, :) = CSpaHat;
        % s = phiC(hatC_, epsilon);
        lambda = 5.6e-10*L;
        hatC_ = sparse(L, L);
        approxError = -1;
    elseif lambdaInit <0
        error('Invalid lambdaInit');
    else
        lambda = lambdaInit;
        spaLambdaHat = SPAselect(X, N);
        if options.verbose
            disp('SPA chooses: ');
            disp(spaLambdaHat);
        end
        [CSpaHat, ~] = computeApproxError(X, spaLambdaHat, ...
            'maxIters', 100);
        approxError = norm(X - X(:, spaLambdaHat)*CSpaHat, 'fro');
        L = size(X, 2);
        hatC_  = sparse(L, L);
        hatC_(spaLambdaHat, :) = CSpaHat;
    end
end

