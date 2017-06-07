function [ W, Z ] = truncate_hada_randomized(W1, Z1, W2, Z2, r, p)
    m = size(W1, 1);
    n = size(Z1, 1);
    if m ~= size(W2, 1) || n ~= size(Z2, 1)
        error('The size of the two matrices must coincide')
    end
    
    if size(W1, 2) < size(W2, 2)
        t = W1; W1 = W2; W2 = t;
        t = Z1; Z1 = Z2; Z2 = t;
        clear t;
    end
    
    Y = zeros(m, r + p);
    for ii = 1:r + p
        Y1 = Z1' * bsxfun(@times, Z2, randn(n, 1));
        Y(:, ii) = sum((W1 * Y1) .* W2, 2);
    end
    [W, ~] = qr(Y, 0);
    
    Z = zeros(n, r + p);
    for ii = 1:r + p
        Y1 = W1' * bsxfun(@times, W2, W(:, ii));
        Z(:, ii) = sum((Z1 * Y1) .* Z2, 2);
    end
    
    [U, S, V] = svd(Z, 'econ');
    W = W * V(:, 1:r);
    Z = U(:, 1:r) * S(1:r, 1:r);
end
