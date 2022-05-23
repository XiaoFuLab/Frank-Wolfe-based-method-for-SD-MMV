function [ theta, B , Tracking] = SPOC_pruned( A, K )

Tracking = struct;
Tracking.lambdaHat = [];
[U,L] = eigs(A,K,'la');


% ==============================
% Modified by Tri Nguyen
% ==============================
U_original = U;
fprintf('Before truncating: %d \n', size(U, 1));
U = truncateData(U, 0.99);
fprintf('After truncating: %d \n', size(U, 1));
% ==============================
% End modifying
% ==============================

Up = U';
F  = zeros(K);
for k=1:K
    [~,i] = max(sum(Up.^2));
    Tracking.lambdaHat(end+1) = i;
    F(:,k) = U(i,:)';
    Q = orth(F(:,1:k)); P = eye(K) - Q*Q';
    Up = P*Up;
end

B = F'*L*F;

% ==============================
% Modified by Tri Nguyen
% ==============================
% theta = F\U';
theta = F\U_original';

% ==============================
% End modifying
% ==============================



theta = max(eps,theta);
theta = bsxfun( @times, theta, 1./sum(theta) );
theta = theta';
