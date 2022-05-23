function [ M ] = anima_tensor( A, sub, F )%, a )

R = 1:size(A,1);
R(sub) = [];

N = length(sub);
I = sub(1:round(N/3));
J = sub(round(N/3)+1:round(2*N/3));
K = sub(round(2*N/3)+1:end);

% mu_A = full(sum(A(I,R),2));
% mu_B = full(sum(A(J,R),2));
% mu_C = full(sum(A(K,R),2));

% T = 2*a^2*length(R)*cpdgen({mu_A,mu_B,mu_C});
T = zeros(length(I),length(J),length(K));
flag = 0; i = 0;
while flag==0
    if (i+1)*1000>length(R)
        flag = 1;
        n = R(i*1000+1:length(R));
    else
        n = R(i*1000+1:(i+1)*1000);
    end
    T = T + cpdgen({full(A(I,n)),full(A(J,n)),full(A(K,n))});
%     T = T + (a+1)*(a+2)*cpdgen({full(A(I,n)),full(A(J,n)),full(A(K,n))});
%     T = T - a*(a+1)*cpdgen({full(A(I,n)),full(A(J,n)),repmat(mu_C,1,length(n))});
%     T = T - a*(a+1)*cpdgen({full(A(I,n)),repmat(mu_B,1,length(n)),full(A(K,n))});
%     T = T - a*(a+1)*cpdgen({repmat(mu_A,1,length(n)),full(A(J,n)),full(A(K,n))});
    i = i+1;
end



Mc = cpd(T,F);

M = [ Mc{1}', Mc{2}', Mc{3}' ];
M = max(M,eps);
M = bsxfun(@times,M,1./sum(M));