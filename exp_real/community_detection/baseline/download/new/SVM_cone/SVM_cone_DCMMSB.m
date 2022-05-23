function [theta, B] = SVM_cone_DCMMSB(A,k,beta)
%
% SVM_cone of DCMMSB inference (for OCCAM inference, change line as the comment instructs)
% Input: A: Adjacency matrix
%        K: number of clusters
%        beta: regularization parameter that does not use the points whose
%        norm is almong the buttom beta quantile among all the points,
%        default as 0.2
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
    beta = 0.2;
end

n = size(A,1);
[v,e]=eigs(A,k);
deg = sum(A,2);

if mean(deg)>=50
    beta = 0;
end

[M, pure] = SVM_cone(v,k,deg,0,beta);

vp = v(pure,:);
np = diag(1./vecnorm(vp'));

gamma_B_square = diag(np* vp * e * vp'* np);
gamma_B_square(gamma_B_square<=0)=1e-12;
gamma_p = diag(sqrt(gamma_B_square/max(gamma_B_square)));

gamma_theta = M*gamma_p;
tre = 1e-12;

gamma_theta(gamma_theta<tre) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DCMMSB inference
[ theta, gamma ] = normalize_row_l1_s( gamma_theta );
% for OCCAM inference, use the following code instead
% gamma = vecnorm(gamma_theta');
% theta = normalize_row_l2_( gamma_theta );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta(isnan(theta))=0;

rho_tmp = sum(gamma)/(n);
gamma = gamma/rho_tmp;

B = diag(1./gamma(pure))*vp*e*vp'*diag(1./gamma(pure));
B(B<0)=0;

B = B / max(max(B));
end

% normalize rows of theta by l1 norm
function [ theta_l1, d_theta ] = normalize_row_l1_s( theta )
    d_theta = sum(theta,2);
    S = find(d_theta>0);
    theta_l1 = theta;
    theta_l1(S,:) = bsxfun(@times, theta(S,:), 1./(sum(theta(S,:), 2)));
end

function [ theta_l2 ] = normalize_row_l2_( theta )
    d_theta = vecnorm(theta');
    S = find(d_theta>0);
    theta_l2 = theta;
    theta_l2(S,:) = bsxfun(@times, theta(S,:), 1./(vecnorm(theta(S,:)'))');
end



