function coh_score = Coherence(C,N,D);
% input: 
% C: F columns, F topics;
% N: the first N words
% D: the original data corpus, doc x word

% extract indices
epsilon = 0.01;
F = size(C,2);
for f=1:F
   [~,idx]=sort(C(:,f),'descend');
   idxC(:,f)=idx(1:N);
end

% creating the NxN coherence score matrix
for f=1:F

    Coh_Mat = zeros(N,N);
    for n=1:N
        pw1 = nnz(D(:,idxC(n,f)));
        for m=n+1:N
            d = D(:,idxC(n,f))& D(:,idxC(m,f));
            pw1w2 = nnz(d); 
            Coh_Mat(n,m) =  log((pw1w2+epsilon)/pw1);
            
        end
    end
    Coh{f}= Coh_Mat + Coh_Mat';
    
    coh_score_individual(f)=sum(sum(Coh{f}));

end

coh_score= mean(coh_score_individual);