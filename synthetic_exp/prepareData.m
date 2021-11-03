function [X, purePixelSet, S] = prepareData(typeA, typeS, M, N, L, varargin) 
% -------------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('epsilon', nan);
    p.addOptional('SNR', 20);
    p.addOptional('purePixelSet', []);
    p.addOptional('numOutliers', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    if strcmp(typeA, 'UNIFORM')
        % UNIFORM
        A = rand(M, N);
    elseif typeA == 2
        error('NotImplemented')
        % ILL-CONDITIONED
        A = rand(M, N);
        [U, S, V] = svd(A, 'econ');
        alpha = exp(-2*log(10)/(N-1));
        i = [1:N];
        S = diag(alpha.^i);
        A = U*S*V;
    else error('NotImplemented')
    end

    if options.verbose
        fprintf('Cond of A: %.2f\n', cond(A));
        fprintf('Eigen max of A^T*A: %.2f\n', eigs(A'*A, 1));
        fprintf('Eigen min of A^T*A: %f\n', eigs(A'*A, 1, 'sm'));
    end

    if strcmp(typeS, 'MIDDLE_POINT')
        %% MIDDLE POINT
        assert L == nchoosek(N, 2)+ N
        S = zeros(N, L);
        S(:, 1:N) = eye(N);
        t=N+1;
        for i=1:N
            for j=i+1:N
                S(:, t) = 0.5*(S(:, i) + S(:, j));
                t = t+1;
            end
        end

        Y = A*S;
        Noise = zeros(M, L);
        Noise(:, 1:N) = 0;
        center = mean(A, 2);
        Noise(:, N+1:end) = Y(:, N+1:end)-center;
        Noise = options.epsilon*Noise/norm(Noise, 'fro');
        X = Y + Noise;

        indices = randperm(L);
        X = X(:, indices);
        r_pure_pixel_set = [];
        pure_pixel_set = 1:N;
        for i=1:numel(pure_pixel_set)
            r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(i));
        end
        pure_pixel_set = r_pure_pixel_set;

    elseif typeS == 2
        error('Need to check later')
    elseif strcmp(typeS, 'DIRICHLET')
        S = zeros(N, L);
        S(:, 1:N) = eye(N);
        S(:, N+1:end) = dirichletRand(ones(1, N), L-N);
        Y = A*S;
        SNR = 10^(options.SNR/10);
        noise = randn(size(Y)); 
        sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / L / SNR;
        noise = sqrt(sigma2)*noise;
        X = Y + noise;
        if options.numOutliers>0
            X(:, end-options.numOutliers:end) = 50*X(:, end-options.numOutliers:end);
        end

        shuffledInd = randperm(L);
        X = X(:, shuffledInd);
        S = S(:, shuffledInd);

        purePixelSet = find(ismember(shuffledInd, 1:N));
    else error('NotImplemented')
    end
end
