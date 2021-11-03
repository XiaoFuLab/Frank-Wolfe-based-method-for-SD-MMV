% *** Description ***
% 
% Select indices of columns of M using a heursitic based SPA. These columns
% should be `good' candidates to be on the border of the minimum volume 
% ellipsoid containing all columns of M. 
% 
% See Algorithm 3 in N. Gillis and S.A. Vavasis, `Semidefinite Programming 
% Based Preconditioning for More Robust Near-Separable Nonnegative Matrix
% Factorization', arXiv, 2013. 
% 
% ind = SPAselect(M,Nbr) 
%  
% ****** Input ******
% M    : an m-by-n matrix
% Nbr  : number of columns to be extracted. 
% 
% ****** Output ******
% inf  : index set of the extracted columns. 

function ind = SPAselect(M,Nbr); 

[r,n] = size(M); 

ind = []; 

if Nbr > n
    error('Nbr has to be larger than the number of columns'); 
end

while length(ind) < Nbr && max(M(:)) > 0
    t = min(r, Nbr-length(ind)); 
	ind = [ind FastSepNMF(M,t)]; 
    M(:,ind) = 0; 
end