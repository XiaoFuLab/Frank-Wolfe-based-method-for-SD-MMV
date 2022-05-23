function [X, groundTruth] = prepareData(source, ind, varargin) 
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    root = [fileparts(which('prepareData.m')), '/', ...
        '/data/coauthorship'];
    fileName = '';
    if strcmp(source, 'MAG')
        fileName = sprintf('MAG%d_adjacency', ind);
        path2File = sprintf('%s/%s.txt', root, fileName);
        path2GT = sprintf('%s/MAG%d_community.txt', root, ind);
        X = read_data(path2File);
        groundTruth = read_ground_truth(path2GT);
    elseif strcmp(source, 'DBLP')
        path2File = sprintf('%s/DBLP%d_adjacency.txt', root, ind);
        path2GT = sprintf('%s/DBLP%d_community.txt', root, ind);
        X = read_data(path2File);
        groundTruth = read_ground_truth(path2GT);
    end
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
