function [ W, Z ] = truncate_randomized_rank1( A, r, p, exponent)
%%% Algorithm 1
    if(nargin < 4)
        exponent = 0;
    end
    Omega0 = randn(ceil(size(A, 2)^.5), r+p);
    Omega = kr(Omega0,Omega0);
    Omega = Omega(1:size(A, 2), :);
    Y = A * Omega;
    if exponent > 0
        G = A * A';
    end
    for i = 1:exponent
        Y = G * Y;
    end
    [W, ~] = qr(Y, 0);
    Z = A' * W;
    [U, S, V] = svd(Z, 'econ');
    W = W * V(:, 1:r);
    Z = U(:, 1:r) * S(1:r, 1:r);
end

