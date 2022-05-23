%% Clear all things
clc; clear; close all; path(pathdef);

addpath('./../../utils/')
addpath('~/code/matlab/common/tensorlab')
addpath('./../SPOC/')

load networks.mat

src = zeros(4,7);
src2 = zeros(4,7);
ttt = zeros(4,7);
err = zeros(4,7);
Src = zeros(4,7);
nz = zeros(1, 7);


%% DBLP
for t = 1:5
    t
    
    A = []; C = [];
    for s = [1:t-1,t+1:5]
        A = blkdiag(A,adjacency{s});
        C = blkdiag(C,community{s});
    end
    
    k = size(C,2);
    
    A = A + A';
    
    N = size(A,1);
    n = 1000;
    
    tic
    [M1,~] = GeoNMF(A,k);
    ttt(1,t) = toc;
    
    tic
    [M2,~] = SPOC(A,k);
    ttt(2,t) = toc;
    
    G = randperm(N,n);
    a = full(sum(A));
    a(G)  = 0;
    [~,H] = maxk(a,k-1);
    % a(H) = 0;
    % [~,G] = maxk(a,n);


    tic
    M3 = anima_tensor( A, G, k );
    ttt(3,t) = toc;


    % Huang&Fu
    tic

    R = 1:N; R(G) = []; R(H) = [];
    Y = full(A(R,H)'*A(R,G));
    [M4,~] = CDMVS(Y);
    ttt(4,t) = toc;

    tmp = abs(Y);
    nz(t) = sum(sum(tmp, 1) == 0);

    % evaluate
    M0 = C(G,:)';

    M1 = M1(G,:)';
    M2 = M2(:,G);
    
% ============ START DEBUG ============

    M0 = C';
    M1 = M1';
    M2 = M2;
% ============  END DEBUG  ============



    R = corr(M0',M1','type','Spearman'); p = Hungarian(-R);
    src(1,t) = sum(diag(p'*R));
    src2(1,t) = mean(diag(p'*R));

    R = corr(M0',M2','type','Spearman'); p = Hungarian(-R);
    src(2,t) = sum(diag(p'*R));
    src2(2,t) = mean(diag(p'*R));

    R = corr(M0',M3','type','Spearman'); p = Hungarian(-R);
    src(3,t) = sum(diag(p'*R));
    src2(3,t) = mean(diag(p'*R));

    R = corr(M0',M4','type','Spearman'); p = Hungarian(-R);
    src(4,t) = sum(diag(p'*R));
    src2(4,t) = mean(diag(p'*R));
    
    P = Hungarian(-M1*M0'); M1 = P'*M1;
    err(1,t) = norm(M1(:)-M0(:),1);
    P = Hungarian(-M2*M0'); M2 = P'*M2;
    err(2,t) = norm(M2(:)-M0(:),1);

    P = Hungarian(-M3*M0'); M3 = P'*M3;
    err(3,t) = norm(M3(:)-M0(:),1);
    P = Hungarian(-M4*M0'); M4 = P'*M4;
    err(4,t) = norm(M4(:)-M0(:),1);
    %
    
    Src(1, t) = computeSRC(M0', M1');
    Src(2, t) = computeSRC(M0', M2');
    Src(3, t) = computeSRC(M0', M3');
    Src(4, t) = computeSRC(M0', M4');
end


%% MAG

% clear

% for t = 6:7
%     
%     A = adjacency{t};
%     C = community{t};
%     
%     k = size(C,2);
%     
%     A = A + A';
%     
%     N = size(A,1);
%     n = 1000;
%     
%     tic
%     [M1,~] = GeoNMF(A,k);
%     ttt(1,t) = toc;
%     
%     tic
%     [M2,~] = SPOC(A,k);
%     ttt(2,t) = toc;
%     
%     
%     a = full(sum(A));
%     [~,H] = maxk(a,k-1);
%     a(H) = 0;
%     [~,G] = maxk(a,n);
%     
%     
%     tic
%     M3 = anima_tensor( A, G, k );
%     ttt(3,t) = toc;
%     
%     
%     % Huang&Fu
%     tic
%     
%     R = 1:N; R(G) = []; R(H) = [];
%     Y = full(A(R,H)'*A(R,G));
%     [M4,~] = CDMVS(Y);
%     ttt(4,t) = toc;
%     
%     
%     % evaluate
%     M0 = C(G,:)';
%     
%     M1 = M1(G,:)';
%     M2 = M2(:,G);
%     
%     R = corr(M0',M1','type','Spearman'); p = Hungarian(-R);
%     src(1,t) = sum(diag(p'*R));
%     R = corr(M0',M2','type','Spearman'); p = Hungarian(-R);
%     src(2,t) = sum(diag(p'*R));
%     R = corr(M0',M3','type','Spearman'); p = Hungarian(-R);
%     src(3,t) = sum(diag(p'*R));
%     R = corr(M0',M4','type','Spearman'); p = Hungarian(-R);
%     src(4,t) = sum(diag(p'*R));
%     
%     P = Hungarian(-M1*M0'); M1 = P'*M1;
%     err(1,t) = norm(M1(:)-M0(:),1);
%     P = Hungarian(-M2*M0'); M2 = P'*M2;
%     err(2,t) = norm(M2(:)-M0(:),1);
%     P = Hungarian(-M3*M0'); M3 = P'*M3;
%     err(3,t) = norm(M3(:)-M0(:),1);
%     P = Hungarian(-M4*M0'); M4 = P'*M4;
%     err(4,t) = norm(M4(:)-M0(:),1);
%     
%     
% end
%
ttt
src
err
