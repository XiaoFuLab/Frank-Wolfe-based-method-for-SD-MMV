function [ Count ] = UniequeWordsCount( C,N )
% input:
% C: F columns, F topics;
% N: the first N words
% D: the original data corpus, doc x word

F = size(C,2);
for f=1:F
    [~,idx]=sort(C(:,f),'descend');
    idxC(:,f)=idx(1:N);
end

Count = zeros(F,F);
for f=1:F
    for k=f+1:F
        Count(f,k) = length(intersect(idxC(:,f),idxC(:,k))); 
    end 
end
Count = Count + Count';

