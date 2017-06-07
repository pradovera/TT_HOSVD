function [W, Z] = truncate_randomized_rank1(A, r, p)
%%% Algorithm 1
    n = size(A, 2);
    Omega0 = randn(ceil(n^.5), r+p);
    Omega = kr(Omega0,Omega0);
    Omega = Omega(1:n, :);
    [W, ~] = qr(A * Omega, 0);
    Z = A' * W;
    [U, S, V] = svd(Z, 'econ');
    W = W * V(:, 1:r);
    Z = U(:, 1:r) * S(1:r, 1:r);
end

