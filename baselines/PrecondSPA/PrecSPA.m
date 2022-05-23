% SDP-based Preconditioning 
% 
% *** Description ***
% First, the semidefinite programming based preconditioning for 
% near-separable matrices is computed. Then, the successive projection 
% algorithm (SPA) is run on the preconditioned matrix (see FastSepNMF.m). 
%
% See Algorithm 2 in  N. Gillis and S.A. Vavasis, `Semidefinite Programming 
% Based Preconditioning for More Robust Near-Separable Nonnegative Matrix
% Factorization', arXiv, 2013. 
% 
% [K,Q,u] = PrecSPA(M,r,normalize,dimred,affi)
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
%              (default = 0.)
% dimred     : parameter allowing to choose the dimensionality reduction
%              technique used (for m > r): 
%              dimred=1: use the truncated SVD, that is, svds(M,r). 
%              dimred=2: use the truncated SVD, but computation based on MM'.
%              dimred~=1,2: use a random projection (u = randn(m,r))
%              (default = 1.)
% affi       : affi=1 displays the differents steps with running times  
%              affi=0 no displays 
%              (default = 1.)
% 
% ****** Output ******
% K        : index set of the extracted columns. 
% Q        : the preconditioning (M <- Q*M).  
% u        : the linear dimensionality reduction (M <- u'*M). 

function [K,Q,u] = PrecSPA(M,r,normalize,dimred,affi)

[m,n] = size(M); 
if r > min(m,n)
    error('r has to be smaller than m and n'); 
end
if nargin <= 2, normalize = 0; end
if normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags((sum(M).^(-1))', 0, n, n); M = M*D; 
end
if nargin <= 3, 
    dimred = 1; 
    % If m << n, use SVD of MM^T. It is computationally much cheaper although
    % it is numerically much less stable (condition number squared). 
    if m^2 <= n, dimred = 2; end
    % If m,n,r are too large, use a random projection
    if m*n*r*min(m,n) >= 1e12, dimred = 3; end
end
if nargin <= 4, affi = 1; end

% Linear dimensionality rediction
etim = cputime; 
if affi == 1
    fprintf('Linear dimensionality reduction started...\n')
end
[Mr, u] = lindimred(M,r,dimred); 
if affi == 1
    fprintf(' Done. \n')
    fprintf('The linear dimensionality reduction took %2.2f seconds. \n',cputime-etim); 
end

% Preconditioning based on minimum volume ellipsoid 
Q = minvolell(Mr,1e-6,r,affi); 

% Preconditioned SPA
K = FastSepNMF(Q*Mr,r); 