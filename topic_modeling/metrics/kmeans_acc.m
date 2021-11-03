function [Acc, rand_index, match] = kmeans_acc(ground_truth, pred, N, varargin) 
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

    topics = unique(ground_truth);
    for i=1:N
        ground_truth(ground_truth == topics(i)) = i;
    end
    assert(min(pred) == 1);
    assert(max(pred) == N);
    [Acc,rand_index,match] = AccMeasure(ground_truth, pred);
end
