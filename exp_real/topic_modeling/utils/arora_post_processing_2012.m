function [C, R, W] = arora_post_processing_2012(P, J, r, D, varargin) 
% ----------------------------------------------------------------------
%% INPUT
%
% P: word by word matrix
% J: set of pure pixel indices
% r: rank
% D: vocab by document, where D is used in model: D = CW
%   Entry of D can be any measurement, such as tf, tf-idf

%% OUTPUT:
%
% C: vocab by topic
% R: WW'
% W: topic by document

% For more details: Learning Topic Models — Going beyond SVD
%
% Author: Tri Nguyen (nguyetr9@oregonstate.edu)
% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    % ==============================
    % ARORA part of identifying A and R
    % ==============================
    
    [nWords,~]=size(P);
    P_perm(1:r,:)= P(J,:);
    row_perm = [1:nWords];
    row_perm(J)=[];
    P_perm(r+1:nWords,:)= P(row_perm,:);
    
    P_perm_perm(:,1:r) = P_perm(:,J);
    col_perm = row_perm;
    P_perm_perm(:,r+1:nWords)= P_perm(:,col_perm);
    
    P_new = P_perm_perm;
    
    Qss = P_new(1:r,1:r);
    Qs = P_new(1:r,:);
    ps = Qs*ones(nWords,1);
    z = (Qss'*Qss)\Qss'*ps;
    At = (Qss*diag(z))\Qs;
    R = diag(z)*Qss*diag(z);
    A = At';

   
    C = A; % Convert to our notation
    % ==============================
    % Given topic model: D=CW, find W as
    % min_{W>0} ||D - CW||_F^2
    % ==============================
    W = nnlsHALSupdt(D, A, D\A);

end
