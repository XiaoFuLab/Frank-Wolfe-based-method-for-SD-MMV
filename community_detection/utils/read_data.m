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
