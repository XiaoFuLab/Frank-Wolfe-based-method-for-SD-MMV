function r = dirichlet_rnd(a, n)
    a = a(:)';
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
    r = r';
end

