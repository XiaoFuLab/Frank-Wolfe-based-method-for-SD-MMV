function [src,M1] = getSRC(M1,M)
    k=size(M,1);
    R = corr(M',M1','type','Spearman');
    p = hungarian(-R);
    src = sum(diag(p*R))/k;
    M1  = p*M1;
end
