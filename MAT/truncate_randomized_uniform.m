function [W, Z] = truncate_randomized_uniform(A, r, p)
%%% Algorithm 1
    Omega = rand(size(A, 2), r+p);
    [W, ~] = qr(A * Omega, 0);
    Z = A' * W;
    [U, S, V] = svd(Z, 'econ');
    W = W * V(:, 1:r);
    Z = U(:, 1:r) * S(1:r, 1:r);
end

