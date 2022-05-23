% Example of Preconditioned SPA on the Middle Points Experiment
% 
% from Section 5.2. in N. Gillis and W.-K. Ma, `Enhancing Pure-Pixel 
% Identification Performance via Preconditioning', arXiv, 2014. 

clear all; clc; close all; 
m = 40; r = 20; 

disp('****Don''t forget to run cvx_setup****'); 

% Noise level 
numnoiselev = 10;
maxnoiselevel = 0.6; 
error = 0:maxnoiselevel/numnoiselev:maxnoiselevel; 

% Number of matrices generated for each noise level
iter = 1;  

% Results 
T = length(error); 
nalgo = 6; 
results = zeros(nalgo,T);  
timings = zeros(1,nalgo); 
fprintf('Total number of tested noise levels: %2.0f \n', T); 
for i = 1 : T
    for k = 1 : iter
        % Generate M randomly : 
        % - W is randomly generated with rand
        W = rand(m,r); 
        
        % - H is generated so that the first r columns of M are the columns 
        % of W, and the remaining ones are on the middle point of any 
        % combination of two columns of W
        H = [eye(r) nchoose2(r)/2]; 
        [r,n] = size(H); 
        
        % The noise moves the data points toward the outside of conv(W)
        M = W*H; 
        Noise = error(i)*[zeros(m,r) M(:,r+1:end)-repmat(mean(W,2),1,n-r)]; 

        % Input matrix 
        M = W*H + Noise;
        % Run all algorithms 
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
            results(algoix,i) = results(algoix,i) + rcur/iter/r; 
        end
    end
    fprintf('%1.0f...',i);
    if mod(i,10) == 0, fprintf('\n'); end
end

% Display performance
scriptfig; 

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