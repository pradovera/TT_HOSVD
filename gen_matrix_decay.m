function [Z, U, S, V] = gen_matrix_decay(m, n, sigmas)
    [U, ~] = qr(randn(m, length(sigmas)), 0);
    [V, ~] = qr(randn(n, length(sigmas)), 0);
    S = diag(sigmas);
    Z = bsxfun(@times, U, sigmas) * V';
end

