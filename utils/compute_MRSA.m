function [out, M] = compute_MRSA(W_hat, W) 
    
    % In Gillis' paper
    
    K = size(W, 2);
    W_hat_bar = mean(W_hat, 2);
    W_bar = mean(W, 2);
    for i=1:K
        for j=1:K
            tmp1 = W_hat(:, i) - W_hat_bar;
            tmp2 = W(:, j) - W_bar;
            tmp = tmp1'*tmp2/norm(tmp1, 'fro')/norm(tmp2, 'fro');
            tmp = min(tmp, 1);
            tmp = max(tmp, -1);
            R(i,j) = acos(tmp);
        end
    end
    M = matchpairs(R, 1000000);
    ind = sub2ind(size(R), M(:, 1), M(:, 2));
    out = mean(R(ind));
end

