% Example of Preconditioned SPA on the 2-by-3 ill-conditioned matrix
% 
% from Section 5.1. in  N. Gillis and W.-K. Ma, `Enhancing Pure-Pixel 
% Identification Performance via Preconditioning', arXiv, 2014. 

clear all; clc; close all;

disp('****Don''t forget to run cvx_setup****');

% Dimensions 
r = 2;m = r;
% Parameter k  
k = 10; 
% You can also try
% k = 100;  % or k = 1000;  

% Noise levels
error = logspace(-2*(log10(k)+1),0,10);

% Results
T = length(error);
nalgo = 6;
results = zeros(nalgo,T);
timings = zeros(1,nalgo);
fprintf('Total number of tested noise levels: %2.0f \n', T);

for i = 1 : T
        W = eye(r) + k*ones(r);
        H = [eye(r) nchoose2(r)/2];
        
        [r,n] = size(H);
        M = W*H;
        [m,n] = size(M);
        
        normmatM = norm(M,'fro');
        
        Noise = [-W M(:,r+1:end)];
        Noise = error(i)*Noise;
        
        % Input matrix
        M = W*H + Noise;
        
        for algoix = 1 : nalgo
            e = cputime; 
            if algoix == 1 % SPA
                K = FastSepNMF(M,r); 
            elseif algoix == 2 % SDP-SPA
                [K,Q,u] = PrecSPA(M,r,0,1,0); 
            elseif algoix == 3 % PW-SPA
                [uh,sh,vh] = svds(M,r); 
                K = FastSepNMF(vh',r); 
            elseif algoix == 4 % SPA-SPA
                K = FastSepNMF(M,r); 
                [uh,sh,vh] = svds(M(:,K),r);  
                K = FastSepNMF( diag(diag(sh).^(-1))*uh'*M,r);
            elseif algoix == 5 % VCA 
                [A_est, K] = VCA(M,'Endmembers',r,'verbose','off'); 
            elseif algoix == 6 % XRAY
                K = FastConicalHull(M,r); K = K'; 
            end
            timings(algoix) = timings(algoix) + cputime-e;
            rcur = sum(unique(K) <= r);
            
            results(algoix,i) = results(algoix,i) + rcur/r;
            
        end
    
    fprintf('%1.0f...',i);
    if mod(i,10) == 0, fprintf('\n'); end
end

% Display performance
scriptfiglog;

% Display total time and robustness
for i = 1 : nalgo
    [a,b] = min(results(i,:) >= 1-1e-16);
    if b == 1, rob(i) = 0; else rob(i) = error(b-1);
    end
end
fprintf('\n');

fprintf('  Algorithm    | Total time (s.) | Robustness \n')
fprintf('--------------------------------------------------- \n')
fprintf(' SPA           |       %3.2f      |   %2.2f  \n', timings(1), rob(1));
fprintf(' SDP-SPA       |       %3.2f      |   %2.2f  \n', timings(2), rob(2));
fprintf(' PW-SPA        |       %3.2f      |   %2.2f  \n', timings(3), rob(3));
fprintf(' SPA-SPA       |       %3.2f      |   %2.2f  \n', timings(4), rob(4));
fprintf(' VCA           |       %3.2f      |   %2.2f  \n', timings(5), rob(5));
fprintf(' XRAY          |       %3.2f      |   %2.2f  \n', timings(6), rob(6));