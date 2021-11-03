function [ A,R,J ] = arora( P,r,ver_opt )
% r: rank of the model, or number of clusters
% P: the correlation matrix, word by word
% ver_opt = 2012 or 2013, two versions of finding A


% Arora's high-level algorithm looks like this


% Let's run Gillis' SPA
% addpath sunsal


[nWords,~]=size(P);
normalize = 1;

M = P;

%==================================================================
% stage 1: find the anchor words - any separable NMF algorithm
%==================================================================


if nWords<1000
    disp('not enough dimension for projection')
    [J,normM,U] = FastSepNMF(M,r,normalize); % This is SPA
else
    
%     ------------ or 'FastAnchorWord by Arora himself'--------
    
%     generate a number

    m_proj = 1000;
    Uproj = randn(m_proj,nWords);
    
    D = spdiags(1./(sum(P)+eps)', 0, nWords, nWords);
    
    P_normalize = P*D;
    
    M = Uproj*P_normalize;
    
    
    
    [J,normM,U] = FastSepNMF(M,r,0);
    
    JJ = J;
    for j=1:r
        A_j = P_normalize(:,J);
        A_j(:,j)=[];
        [~,i] = max(sum( (P_normalize-A_j*inv(A_j'*A_j)*(A_j'*P_normalize)).^2 ));
        JJ(j)= i;
    end
    
     
    
    
end


%========================================
% recover A and R (notations in Arora'13)
%==========================================

% Original Recover (closed-form)

% permute the rows and columns of P;

if ver_opt == 2012
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
elseif ver_opt == 2013
    
    D = spdiags(((sum(P')+1e-16).^(-1))', 0, nWords, nWords); P_normalize = D*P; 

    Q_bar = P_normalize; % taking transpose since we are asking for the row normalization
    Q_S = Q_bar(J,:);% taking out the basis,
   
    theta = (Q_S*Q_S')'\(Q_S*Q_bar');
    U = zeros(r,nWords) ;
    %% LS
    
%     [Ct]= ADMM_Theta(Q_S',Q_bar',theta,U,r,ones(1,nWords));
%     C = Ct';
    
     [Ct,e] = FGMfcnls(Q_bar',Q_S',theta,5000);
%      Ct = nnlsHALSupdt(Q_bar',Q_S',theta,5000);
      C = Ct';
    
    
%     cvx_begin
%     variable Ct(J,nWords)
%     minimize norm(Q_bar'-Q_S'*Ct,'fro')
%     subject to 
%     ones(1,J)*Ct==ones(1,nWords);
%     Ct>=0;
%     cvx_end
%     
    
    %

    
    
    pw = ((sum(P')+1e-16))';% P*ones(nWords,1);
    Dpw = spdiags(pw,0,nWords,nWords);
    A_prime = Dpw*C;
    DD= spdiags(((sum(A_prime)+1e-16).^(-1))', 0, r, r); A = A_prime*DD; 
    Ap = (A'*A)\A';
    R = Ap*P*Ap';
    
end

end

function [theta,U]= ADMM_Theta(B,Y,theta,U,K,TrR)
rho = 1/eigs(B'*B,1);
Big = [2*B'*B + rho*eye(K), ones(K,1);ones(1,K),0];
inverseBig = inv(Big);

[K,M]=size(theta);
X = theta;
U = X;
Z = 0;
Maxiter = 2000;
ob=[];
for i=1:Maxiter
    XX =inverseBig*[2*B'*Y+rho*(Z-U);TrR];
    X = XX(1:K,:);
    Z = max(X + U,0);
    U = U + X - Z;
    if sum(sqrt(sum((X-Z).^2)))<=1e-10;
        break;
    end
    
    
end


theta = X;
end
