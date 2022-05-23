function [ theta, B ] = SPOC( A, K )

[U,L] = eigs(A,K,'la');

Up = U';
F  = zeros(K);
for k=1:K
    [~,i] = max(sum(Up.^2));
    F(:,k) = U(i,:)';
    Q = orth(F(:,1:k)); P = eye(K) - Q*Q';
    Up = P*Up;
end

B = F'*L*F;
theta = F\U';

theta = max(eps,theta);
theta = bsxfun( @times, theta, 1./sum(theta) );