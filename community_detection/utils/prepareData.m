function [X, groundTruth, G, Addition] = prepareData(source, ind, varargin) 
% ----------------------------------------------------------------------
% 
% Return:
%   - X: adjacency matrix
%   - groundTruth: true theta by normalizing over ...


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('isSecondOrder', false);
    p.addOptional('usingSubMatrix', true);
    p.addOptional('trialInd', 0);
    p.addOptional('removeDup', 1);
    p.addOptional('energyReserve', 1);
    
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    Addition = struct;
    root = [fileparts(which('prepareData.m')), '/..', ...
        '/data/coauthorship'];
    fileName = '';
    if strcmp(source, 'MAG')
        if options.trialInd == 0
            fileName = sprintf('MAG%d_adjacency', ind);
            path2File = sprintf('%s/%s.txt', root, fileName);
            path2GT = sprintf('%s/MAG%d_community.txt', root, ind);
            X = read_data(path2File);
            groundTruth = read_ground_truth(path2GT);
        else
            fileName = sprintf('MAG%d_trial_%d', ind, options.trialInd);
            Data = load(sprintf('%s/random-split2/%s.mat', root, fileName));
            X = Data.X;
            groundTruth = Data.GT;
        end
    elseif strcmp(source, 'DBLP')
        if ind <= 5
            if options.trialInd == 0
                path2File = sprintf('%s/DBLP%d_adjacency.txt', root, ind);
                path2GT = sprintf('%s/DBLP%d_community.txt', root, ind);
                X = read_data(path2File);
                groundTruth = read_ground_truth(path2GT);
            else
                fileName = sprintf('DBLP%d_trial_%d', ind, options.trialInd);
                Data = load(sprintf('%s/random-split2/%s.mat', root, fileName));
                X = Data.X;
                groundTruth = Data.GT;
            end
        elseif ind > 5 && ind <= 10
            error('Not use anymore !');
            X = [];
            groundTruth = [];
            ind = ind - 5;
            for s = [1:ind-1,ind+1:5]
                path2File = sprintf('%s/DBLP%d_adjacency.txt', root, s);
                path2GT = sprintf('%s/DBLP%d_community.txt', root, s);
                X = blkdiag(X, read_data(path2File));
                groundTruth = blkdiag(...
                    groundTruth, read_ground_truth(path2GT));
            end
        elseif ind>10 && ind <= 14
            error('Not use anymore !');
            candidates = [5 10 15 19];
            tmp = candidates(ind-10);
            Data = load(sprintf('%s/random-split/%s/%d/data.mat', ...
                root, num2str(tmp), options.trialInd));
            X = Data.Xr;
            groundTruth = Data.GTr;
        else
            error('Check it')
        end
    elseif strcmp(source, 'EMAIL')
        root = [fileparts(which('prepareData.m')), '/..', '/data/email-eu-core/'];
        if ind>0
            Data = load(sprintf('%s/random-split/%s/%d/data.mat', ...
                root, num2str(ind), options.trialInd));
            X = Data.Xr;
            groundTruth = Data.GTr;
        else
            [X, theta] = read_email_dataset(root);
            numCommunities = max(theta);
            numNodes = size(X, 1);
            groundTruth = sparse(1:numNodes, theta, 1, numNodes, numCommunities);
        end
    elseif strcmp(source, 'YOUTUBE')
        root = [fileparts(which('prepareData.m')), '/..', '/data/youtube/'];
        Data = load(sprintf('%s/random-split/%s/%d/data.mat', ...
            root, num2str(ind), options.trialInd));
        X = Data.Xr;
        groundTruth = Data.GTr;
    else
        error('NotImplemented')
    end

    G = nan;
    if options.isSecondOrder
        numGroups = size(groundTruth, 2);
        Addition.originalX = X;
        Addition.originalGt = groundTruth;
        [X, G] = getSecondOrder(X, numGroups);
        groundTruth = groundTruth(G, :);

        if ~options.usingSubMatrix
            X = Addition.originalX;
        end
    end

    
    % ==============================
    % New
    % ==============================
    fprintf('Before cut-off: %s \n', num2str(size(X, 1)));
    if options.energyReserve < 0 || options.energyReserve > 1
        error('energyReserve should be within 0-1');
    end

    colEnergy = vecnorm(X, 2, 1).^2;
    sortColEnergy = sort(colEnergy, 'descend');
    totalEnergy = sum(colEnergy);
    sumEnergy = 0;
    threshold = 0;
    L = size(X, 1);
    for i=1:L
        sumEnergy = sumEnergy + sortColEnergy(i);
        if sumEnergy > options.energyReserve*totalEnergy
            threshold = sortColEnergy(i);
            break
        end
    end
    zeroFiltering = find(colEnergy < threshold);

    X(zeroFiltering, :) = [];
    X(:, zeroFiltering) = [];
    groundTruth(zeroFiltering, :) = [];
    fprintf('After cut-off: %s \n', num2str(size(X, 1)));
    

end

function [X, GT] = random_select(X, GT, numCommunities, varargin) 
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

    N = size(GT, 2);
    L = size(X, 1);
    communities = randsample(N, numCommunities);
    ind = [];
    for i=1:L
        if sum(GT(i, communities)) > 0
            ind(end+1) = i;
        end
    end
    X = X(ind, ind);
    GT = GT(ind, communities);
end


function [A, theta] = read_email_dataset(root, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    path2File = sprintf('%s/email-Eu-core.txt', root);
    X = readmatrix(path2File);
    X = X+1;
    num_nodes = max(X(:));
    A = [X; [X(:, 2) X(:, 1)]];
    A = sparse(A(:, 1), A(:, 2), 1, num_nodes, num_nodes);
    for i=1:num_nodes
        for j=1:i
            if A(i, j)>0
                A(i, j)=1;
                A(j, i)=1;
            end
        end
    end

    % Ground truth
    path2File = sprintf('%s/email-Eu-core-department-labels.txt', root);
    theta = readmatrix(path2File);
    theta = theta + 1;
    assert(num_nodes == max(theta(:, 1)));
    num_communities = max(theta(:, 2));
    theta = theta(:, 2);
end

function [out, G] = getSecondOrder(X, numGroups, varargin) 
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

    n = 5000;
    N = size(X, 1);
    G = randperm(N, n);
    a = sum(X);
    a(G) = 0;
    k = numGroups;
    [~,H] = maxk(a,k);

    R = 1:N; R(G) = []; R(H) = [];
    out = X(R,H)'*X(R,G);
end

function [A] = read_data(path_to_file, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    X = readmatrix(path_to_file);
    num_nodes = max(X(:));
    A = [X; [X(:, 2) X(:, 1)]];
    A = sparse(A(:, 1), A(:, 2), 1, num_nodes, num_nodes);
end


function [theta] = read_ground_truth(path_to_file, varargin) 
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
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    theta = readmatrix(path_to_file);
    num_nodes = max(theta(:, 1));
    num_communities = max(theta(:, 2));

    theta = sparse(theta(:, 1), theta(:, 2), theta(:, 3));
    theta = full(theta);
end
