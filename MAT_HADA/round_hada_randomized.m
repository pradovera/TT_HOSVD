function [ W, Z ] = round_hada_randomized(W1, Z1, W2, Z2, epsilon, q)
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
    
    epEff = epsilon^2;
    if epEff < 2 * eps
        warning('Effective precision smaller than machine epsilon');
        epEff = 2 * eps;
    end
    
    Y = zeros(m, q);
    for ii = 1:q
        Y1 = Z1' * bsxfun(@times, Z2, randn(n, 1));
        Y(:, ii) = sum((W1 * Y1) .* W2, 2);
    end
    [W, ~] = qr(Y, 0);
    
    Z = zeros(n, q);
    for ii = 1:q
        Y1 = W1' * bsxfun(@times, W2, W(:, ii));
        Z(:, ii) = sum((Z1 * Y1) .* Z2, 2);
    end
    
    normfro2 = cumsum(sum(Z.^2, 1))';
    
    j = 0;
    while (j == 0) || ((normfro2(j + q) - normfro2(j)) > epEff * normfro2(j + q))
        j = j + 1;
        Y1 = Z1' * bsxfun(@times, Z2, randn(n, 1));
        w = sum((W1 * Y1) .* W2, 2);
        w = w / norm(w);
        for jj = 1:2
            w = w - W * (W' * w);
        end
        w = w / norm(w);
        W = [W, w];
        
        Y1 = W1' * bsxfun(@times, W2, w);
        z = sum((Z1 * Y1) .* Z2, 2);
        Z = [Z, z];
        
        normfro2 = [normfro2; normfro2(end) + sum(z.^2)];
    end
end