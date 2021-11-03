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
