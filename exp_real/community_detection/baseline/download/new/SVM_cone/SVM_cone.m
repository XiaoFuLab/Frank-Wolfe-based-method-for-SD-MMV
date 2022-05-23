function [M, pure] = SVM_cone(Z,k,deg,topic_flag,beta)
%
% SVM_cone: conrer finding and regression algorithm of data that lies on a
% cone
% Input: Z: data matrix that lies on a cone (n by k)
%        k: dimension of the data point
%        beta: regularization parameter that does not use the points whose
%        norm is almong the buttom beta quantile among all the points,
%        default as 0
%        deg: degree of observed data matrix (if observation is network
%        adjacency matrix), 0 for topic models
%        topic_flag: 1 if SVM_cone is used for topic model inference
% Output: M: weight matrix of data points to the corners (Z=MY_p, please see our paper for details)
%         pure: estimated corner indices

% Author: Xueyu Mao 
% Email: maoxueyu@gmail.com
% Last Update: Jan 15, 2019
%
% Reference: Overlapping Clustering Models, and One (class) SVM to Bind
% Them All, in Neural Information Processing Systems, 2018. ArXiv: 1806.06945

if nargin < 5
    beta = 0;
end
    
n = size(Z,1);

row_norm = vecnorm(Z');
mean_row_norm = quantile(row_norm,beta);

if  mean_row_norm < eps
    beta = beta+0.1;
    mean_row_norm = quantile(row_norm,beta);
end

low_norm_idx = find(row_norm<mean_row_norm);
uA = normalize_row_l2_( Z );

% throw out nodes which have very small row norms in Z
uA(low_norm_idx,:)=10*ones(length(low_norm_idx),1)*mean(uA);

if topic_flag
    cond_tre = 150;
elseif mean(deg)>50
    cond_tre = 1.5;
elseif mean(deg)<15
    cond_tre = k*6;
else
    cond_tre = k*3;
end

flag = 1;
iii=0;

% One-class SVM (https://github.com/cjlin1/libsvm)
if k>10
    nu=2*k/n;
else
    nu=k/n;
end

model_A = svmtrain(ones(n,1), uA, strcat(['-s 2 -t 0 -n ', num2str(nu)]));

w = model_A.sv_coef'*model_A.SVs;
b = model_A.rho;

b_y = uA*w';
   
SVs=full(model_A.SVs);
SV_indices = model_A.sv_indices(sum(abs(SVs),2)>0);

epsilon_0 = 0;
epsilon_1 = 0.02;

while flag
    iii=iii+1;
    flag = 0;    
    
    pures_idx = find(b_y<=b*(1+epsilon_0));

    new_SV_indices = unique([SV_indices'  pures_idx' ]');
    new_SVs = uA(new_SV_indices,:);
    
    % K-means clustring to cluster the pure nodes set 
    sumd_min = Inf;
    kmeans_iter = 10;

    for ii = 1:kmeans_iter
        [idx_tmp,~,sumd] = kmeans(new_SVs,k);
        
        if sum(sumd) < sumd_min
            sumd_min = sum(sumd);
            idx = idx_tmp;
        end
    end
    
    pure = [];
    for i=1:k
        tmp=new_SV_indices(idx==i); 
        dist_b = vecnorm((uA(tmp,:)-mean(uA(tmp,:)))');
        [~,id] = sort(dist_b);
        i_p = 1;
        while sum(abs(Z(id(i_p),:)))<=1e-12 && i_p <= length(id)-1
            i_p = i_p + 1;
        end
        pure=[pure tmp(id(i_p))];
    end
    
    vp = Z(pure,:);
    np = diag(1./vecnorm(vp'));
    yp = np*vp;

    if iii >50
        break;
    end
    
    if cond(yp) > cond_tre
        epsilon_0 = epsilon_0 + epsilon_1;
        flag = 1;
    end
end

vp = Z(pure,:);
np = diag(1./vecnorm(vp'));
yp = np*vp;

M = Z/yp;
tre = 1e-12;
M(M<tre) = 0;
M = sparse(M);
end

% normalize rows of theta by l2 norm
function [ theta_l2 ] = normalize_row_l2_( theta )
    d_theta = vecnorm(theta');
    S = find(d_theta>0);
    theta_l2 = theta;
    theta_l2(S,:) = bsxfun(@times, theta(S,:), 1./(vecnorm(theta(S,:)'))');
end




