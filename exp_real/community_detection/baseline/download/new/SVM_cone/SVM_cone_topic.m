function [C, pure] = SVM_cone_topic(A,k,beta)
%
% SVM_cone of topic model inference
% Input: A: Word-document matrix
%        k: number of topics
%        beta: regularization parameter that does not use the points whose
%        norm is almong the buttom beta quantile among all the points,
%        default as 0.1
% Output: theta: node-community matrix, theta_{ij} is the probability node i is in community j
%         B: community-community matrix, B_{ij} is the probability there is
%         an edge between a node in community i and a node in community j

% Author: Xueyu Mao 
% Email: maoxueyu@gmail.com
% Last Update: Jan 15, 2019
%
% Reference: Overlapping Clustering Models, and One (class) SVM to Bind
% Them All, in Neural Information Processing Systems, 2018. ArXiv: 1806.06945

if nargin < 3
    beta = 0.1;
end

H = normalize_row_l1_s( A' );

H = H';

[v,~]=svds(H,k);

[M, pure] = SVM_cone(v,k,0,1,beta);

C = normalize_row_l1_( M' );
C = sparse(C');

end

% normalize rows of theta by l1 norm
function [ theta_l1, d_theta ] = normalize_row_l1_s( theta )
    d_theta = sum(theta,2);
    S = find(d_theta>0);
    theta_l1 = theta;
    theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));
end



