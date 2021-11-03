%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common/')

L = 1000;
S = dirichletRand(10*ones(3, 1), L);
S(:, 1:3) = [
0.9 0.2 0.1;
0.1 0.8 0.1;
0   0.  0.8;
];
C = [S; zeros(L-3, L)];


result = [];
params = [1 0.5 0.1 0.05 0.01 0.005 0.001 1e-4 1e-5 1e-6 1e-7 1e-8 1e-10 1e-12 1e-15];
for epsilon=params
    x1 = phiC(C, epsilon);
    x2 = sum(vecnorm(C, Inf, 2));

    result(end+1) = abs(x1-x2);
end

figure();
loglog(params, result, 'o-');
xlabel('epsilon')
ylabel('$abs( \phi_{\mu}(C) - ||C||_{\infty, 1})$')
exportgraphics(gcf, 'phiC.png', 'resolution', 300);

