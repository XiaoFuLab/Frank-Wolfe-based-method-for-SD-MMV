%% Clear all things
clc; clear; close all; path(pathdef);

% d = 224; L = 5000; k = 10;
d = 10; L = 50; k = 3;
A = rand(d, k);
S = drchrnd(ones(1, k), L);
S(1:k, 1:k) = eye(k);
S(k+1:end, 1:k) = 0;
X = A*S;
SNR = 40;
noise = randn(size(X)); 
sigma2 = sum(vecnorm(A*S, 2, 1).^2) / d / L / SNR;
noise = sqrt(sigma2) *noise;
X = X + noise;
noise_eps = max(vecnorm(noise, 2, 1));
% noise_eps = noise_eps + (rand(1)/6.7 - .15)*noise_eps;
delta = 2*noise_eps;

[lambda, C] = convex_relaxation(X, noise_eps);

figure();
bar(vecnorm(C, 2, 2));
xlabel('Index n')
ylabel('||C(n,:)||_2')
saveas(gcf, 'results/convex-relaxation-demo.eps', 'epsc')
function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
    r = r';
end

