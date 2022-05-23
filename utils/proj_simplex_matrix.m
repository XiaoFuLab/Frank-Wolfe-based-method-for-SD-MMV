function X = proj_simplex_matrix(X)
% project a matrix X to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

    for i=1:size(X, 2)
        X(:, i) = proj_simplex(X(:, i));
    end

end
